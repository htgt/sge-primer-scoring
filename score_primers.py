#!/usr/bin/env python3

import argparse
from os import path

from src.scoring import Scoring


def non_empty_file(arg):
    if not path.isfile(arg):
        raise argparse.ArgumentTypeError(f'File does not exist: {arg}')
    if path.getsize(arg) == 0:
        raise argparse.ArgumentTypeError(f'File is empty: {arg}')
    return arg


def positive_int(arg):
    if int(arg) < 0:
        raise argparse.ArgumentTypeError('Mismatch number cannot be negative')
    return int(arg)


def new_file_path(arg):
    if arg.endswith('/') or path.isdir(arg):
        raise argparse.ArgumentTypeError(
            f'Directory provided rather than file path: {arg}')
    if path.isfile(arg):
        raise argparse.ArgumentTypeError(f'File already exists: {arg}')
    return arg


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Tool to score primer pairs using output from Exonerate iPCRess'),
        epilog=(
            './score_primers.py examples/example_input.txt 4 '
            'examples/example_output.tsv'))
    parser.add_argument(
        'ipcress_file',
        help='File containing output from Exonerate iPCRess',
        type=non_empty_file)
    parser.add_argument(
        'mismatch',
        help='Mismatch number used for Exonerate iPCRess',
        type=positive_int)
    parser.add_argument(
        'output_tsv',
        help='Path for output TSV file',
        type=new_file_path)
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 0.0')
    return parser.parse_args()


def main():
    args = parse_arguments()
    scoring = Scoring(args.ipcress_file, args.mismatch)
    scoring.add_scores_to_df()
    scoring.save_mismatches(args.output_tsv)


if __name__ == '__main__':
    main()
