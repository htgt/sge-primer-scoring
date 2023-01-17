#!/usr/bin/env python3

# Copyright (c) 2022, 2023 Genome Research Ltd.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


import argparse
from os import path

from src.scoring import Scoring


def non_empty_file(arg):
    if not path.isfile(arg):
        raise argparse.ArgumentTypeError(f"File does not exist: '{arg}'")
    if path.getsize(arg) == 0:
        raise argparse.ArgumentTypeError(f"File is empty: '{arg}'")
    return arg


def positive_int(arg):
    if int(arg) < 0:
        raise argparse.ArgumentTypeError('Mismatch number cannot be negative')
    return int(arg)


def new_file_path(arg):
    if arg.endswith('/') or path.isdir(arg):
        raise argparse.ArgumentTypeError(
            f"Directory provided rather than file path: '{arg}'"
        )
    if path.isfile(arg):
        raise argparse.ArgumentTypeError(f"File already exists: '{arg}'")
    return arg


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Tool to score primer pairs using output from Exonerate iPCRess'
        ),
        epilog=(
            './score_primers.py examples/example_ipcress_file.txt'
            ' 4 example_output.tsv'
        ))
    parser.add_argument(
        'ipcress_file',
        help='File containing output from Exonerate iPCRess',
        type=non_empty_file
    )
    parser.add_argument(
        'mismatch',
        help='Mismatch number used for Exonerate iPCRess',
        type=positive_int
    )
    parser.add_argument(
        'output_tsv',
        help='Path for output TSV file',
        type=new_file_path
    )
    parser.add_argument(
        '--targeton_csv',
        help=(
            'CSV of primer pairs and corresponding targetons'
            '- adds targeton column to output'
        ),
        type=non_empty_file
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    scoring = Scoring(args.ipcress_file, args.mismatch, args.targeton_csv)
    scoring.add_scores_to_df()
    scoring.save_mismatches(args.output_tsv)
    print(f"Scoring complete! File saved to '{args.output_tsv}'")


if __name__ == '__main__':
    main()
