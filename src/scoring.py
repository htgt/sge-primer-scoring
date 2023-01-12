from collections import defaultdict
import re
from pathlib import Path

import pandas as pd
import numpy as np


class ScoringError(Exception):
    pass


class Scoring:
    def __init__(self, ipcress_file, mismatches, targeton_csv=None):
        self._mismatch_df = self.mismatches_to_df(
            ipcress_file, mismatches, targeton_csv)
        self._csv = targeton_csv

    @staticmethod
    def mismatches_to_df(ipcress_file, mismatches, targeton_csv=None):
        mismatch_counts = defaultdict(
            lambda: {str(i): 0 for i in range(2 * mismatches + 1)})

        with open(ipcress_file) as ipcress_fh:
            for line in ipcress_fh:
                if line == '-- completed ipcress analysis\n':
                    break
                exp_id, primer_5, mismatch_5, primer_3, mismatch_3 = (
                    Scoring._parse_ipcress_line(line, ipcress_file))
                total_mismatches = str(int(mismatch_5) + int(mismatch_3))
                try:
                    mismatch_counts[(exp_id, primer_5)][mismatch_5] += 1
                    mismatch_counts[(exp_id, primer_3)][mismatch_3] += 1
                    mismatch_counts[(exp_id, 'Total')][total_mismatches] += 1
                except KeyError:
                    raise ScoringError(
                        f'Mismatch number too low: {mismatches}')

        df = pd.DataFrame.from_dict(mismatch_counts, orient='index')
        if df.empty:
            raise ScoringError(f'{ipcress_file}: No data in ipcress file')
        df.index.set_names(['Primer pair', 'A/B/Total'], inplace=True)
        if targeton_csv:
            Scoring._add_targeton_column(df, targeton_csv)
        df.sort_index(inplace=True)  # order A, B, Total
        df['WGE format'] = df.apply(lambda row: row.to_dict(), axis=1)
        return df

    @staticmethod
    def _parse_ipcress_line(line, ipcress_file):
        regex = (
            r'ipcress: \S+ '  # ipcress: 11:filter(unmasked)
            r'(\S+) \d+ '  # SMARCA4_exon24_1 204
            r'([A|B]) \d+ (\d+) '  # A 3231378 4
            r'([A|B]) \d+ (\d+) '  # A 3231564 4
            r'[a-zAB_]+\n'  # single_A
        )
        valid_line = re.fullmatch(regex, line)
        if not valid_line:
            raise ScoringError(f'{ipcress_file}: Invalid ipcress file')
        return valid_line.groups()

    @staticmethod
    def _add_targeton_column(df, targeton_csv):
        targetons = defaultdict(str)
        with open(targeton_csv) as fh:
            for line in fh:
                valid_line = re.fullmatch(r'(\S+),(\S+)\n', line)
                if not valid_line:
                    raise ScoringError(f'{targeton_csv}: Invalid targeton csv')
                primer_pair, targeton = valid_line.groups()
                if (primer_pair in targetons) and (
                        targetons[primer_pair] != targeton):
                    raise ScoringError(
                        f'{targeton_csv}: Conflicting entries'
                        f' in targeton csv for {primer_pair}')
                targetons[primer_pair] = targeton
        df['Targeton'] = df.apply(lambda row: targetons[row.name[0]], axis=1)
        df.set_index('Targeton', append=True, inplace=True)
        df.index = df.index.reorder_levels(['Targeton', 'Primer pair', 'A/B/Total'])

    @property
    def mismatch_df(self):
        return self._mismatch_df

    def add_scores_to_df(self):
        df = self.mismatch_df
        df['Score'] = df.apply(self.score_mismatches, axis=1)
        df['Sum'] = df.groupby('Primer pair')['Score'].transform('sum')
        if self._csv:
            df.sort_values(
                ['Targeton', 'Sum', 'Primer pair', 'A/B/Total'], inplace=True)
        else:
            df.sort_values(['Sum', 'Primer pair', 'A/B/Total'], inplace=True)
        df.drop('Sum', axis=1, inplace=True)

    @staticmethod
    def score_mismatches(row):
        if row.name[-1] != 'Total':
            return np.nan
        weights = {str(i): 10 ** (8 - i) for i in range(2, 9)}
        weights['0'] = 10 ** 10  # fail
        weights['1'] = 10 ** 10  # fail
        score = 0
        for col, val in row.items():
            if col not in weights.keys():
                continue
            if col == '0':
                if val == 0:
                    raise ScoringError(
                        f'No on-target hit found for {row.name[-2]}')
                val -= 1  # take away on-target hit
            score += val * weights[col]
        return score

    def save_mismatches(self, output_file):
        Path(output_file).parent.mkdir(exist_ok=True, parents=True)
        self.mismatch_df.to_csv(output_file, sep='\t')
