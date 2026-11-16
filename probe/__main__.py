"""Predict knockout responses from a CSV of measured inputs, or write a blank input template.

    python -m probe --model models/broad/model.json --input models/broad/example_inputs.csv --output predictions.csv
    python -m probe --model models/broad/model.json --template --output my_inputs.csv
"""
import argparse
import csv
from pathlib import Path

from .model import Model, read_measurements


def main():
    parser = argparse.ArgumentParser(prog='python -m probe', description='Predict gene-knockout responses with a frozen PROBE model.')
    parser.add_argument('--model', type=Path, default=Path('models/broad/model.json'))
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--input', type=Path, help='CSV with gene_id, gene, gene_effect')
    mode.add_argument('--template', action='store_true', help='write a blank input CSV for this model')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        model = Model.load(args.model)
        if args.template:
            rows = [dict(gene_id=g['id'], gene=g['symbol'], gene_effect='') for g in model.probes]
        else:
            rows = sorted(model.predict(read_measurements(args.input.read_text(), model.probes)),
                          key=lambda r: r['predicted_effect'])
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f'Wrote {len(rows)} rows to {args.output}')
    except (ValueError, OSError, KeyError, csv.Error) as error:
        parser.exit(2, f'probe: {error}\n')


if __name__ == '__main__':
    main()
