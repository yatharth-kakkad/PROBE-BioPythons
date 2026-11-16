"""Panel selection and the conditional Gaussian predictor used by both PROBE models."""
from dataclasses import dataclass

import numpy as np


@dataclass
class Panel:
    input_mean: np.ndarray
    input_scale: np.ndarray
    output_mean: np.ndarray
    covariance: np.ndarray        # input-input covariance (standardized inputs)
    cross_covariance: np.ndarray  # input-output covariance
    noise: float                  # ridge regularization

    def __post_init__(self):
        n, t = self.input_mean.size, self.output_mean.size
        for value, shape in [(self.input_mean, (n,)), (self.input_scale, (n,)),
                             (self.output_mean, (t,)), (self.covariance, (n, n)),
                             (self.cross_covariance, (n, t)), (np.asarray(self.noise), ())]:
            if value.shape != shape or not np.isfinite(value).all():
                raise ValueError('Invalid model dimensions or nonfinite parameters')
        if not n or not t or (self.input_scale <= 0).any() or self.noise < 0:
            raise ValueError('Model requires nonempty genes, positive scales and nonnegative regularization')

    def predict(self, values):
        """Rows of input values with NaN for unmeasured inputs. Conditions on the measured subset."""
        values = np.asarray(values, float)
        if values.ndim != 2 or values.shape[1] != len(self.input_mean) or np.isinf(values).any():
            raise ValueError('Expected one column per input gene, NaN for missing')
        out = np.tile(self.output_mean, (len(values), 1))
        for row, x in zip(out, values):
            ids = np.flatnonzero(np.isfinite(x))
            if len(ids):
                a = self.covariance[np.ix_(ids, ids)] + self.noise * np.eye(len(ids))
                row += ((x[ids] - self.input_mean[ids]) / self.input_scale[ids]) @ np.linalg.solve(a, self.cross_covariance[ids])
        return out


def _standardize(x, y, noise):
    if (x.ndim != 2 or y.ndim != 2 or len(x) != len(y) or len(x) < 2 or
            not x.shape[1] or not y.shape[1] or np.isinf(x).any() or np.isinf(y).any() or
            not np.isfinite(x).any(axis=0).all() or not np.isfinite(y).any(axis=0).all()):
        raise ValueError('Training requires aligned matrices with at least two rows and observations in every column')
    if not np.isfinite(noise) or noise < 0:
        raise ValueError('Regularization must be finite and nonnegative')
    mean, scale = np.nanmean(x, axis=0), np.nanstd(x, axis=0)
    if (scale < 1e-8).any():
        raise ValueError('Constant input column')
    a = np.where(np.isfinite(x), (x - mean) / scale, 0)
    y_mean = np.nanmean(y, axis=0)
    b = np.where(np.isfinite(y), y - y_mean, 0)
    return a, b, mean, scale, y_mean


def fit(x, y, noise):
    """Fit the predictor from training inputs x (models × candidate inputs) and outputs y."""
    a, b, mean, scale, y_mean = _standardize(np.asarray(x, float), np.asarray(y, float), noise)
    n = len(a) - 1
    return Panel(mean, scale, y_mean, a.T @ a / n, a.T @ b / n, float(noise))


def select(x, y, budget, noise):
    """Greedy A-optimal selection: each step adds the input that most reduces the summed
    remaining variance of the (variance-normalized) outputs. Returns column indices in order."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    a, b, *_ = _standardize(x, y, noise)
    if not isinstance(budget, (int, np.integer)) or not 0 < budget <= x.shape[1] or noise == 0:
        raise ValueError('Selection requires an integer budget within the input count and positive regularization')
    variance = np.nanvar(y, axis=0)
    if (variance <= 0).any():
        raise ValueError('Panel selection requires varying response columns')
    n = len(a) - 1
    cross = a.T @ (b / np.sqrt(variance)) / n
    diagonal = np.sum(a * a, axis=0) / n
    updates = np.empty((a.shape[1], budget))
    chosen, allowed = [], np.ones(a.shape[1], bool)
    for step in range(budget):
        gain = np.sum(cross * cross, axis=1) / (diagonal + noise)
        j = int(np.argmax(np.where(allowed, gain, -np.inf)))
        chosen.append(j); allowed[j] = False
        residual = a.T @ a[:, j] / n - updates[:, :step] @ updates[j, :step]
        direction = residual / np.sqrt(diagonal[j] + noise)
        cross -= direction[:, None] * (cross[j] / np.sqrt(diagonal[j] + noise))
        diagonal = np.maximum(diagonal - direction ** 2, 0)
        updates[:, step] = direction
    return np.asarray(chosen)
