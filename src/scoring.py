import csv
from collections import defaultdict

import pandas as pd


class Scoring:
    def __init__(self, ipcress_file, mismatches):
        self._mismatch_df = self.mismatches_to_df(ipcress_file, mismatches)

    @staticmethod
    def mismatches_to_df(ipcress_file, mismatches):
        ipcress_fields = [
            'ipcress', 'seq', 'id', 'length',
            'primer_5', 'pos_5', 'mismatch_5',
            'primer_3', 'pos_3', 'mismatch_3',
            'description'
        ]
        mismatch_counts = defaultdict(
            lambda: {str(i): 0 for i in range(2 * mismatches + 1)})

        with open(ipcress_file, newline='') as ipcress_fh:
            reader = csv.DictReader(
                ipcress_fh, fieldnames=ipcress_fields, delimiter=' ')
            for row in reader:
                if row['ipcress'] == '--':
                    break
                mismatch_counts[(
                    row['id'], row['primer_5'])][row['mismatch_5']] += 1
                mismatch_counts[(
                    row['id'], row['primer_3'])][row['mismatch_3']] += 1
                total_mismatches = str(
                    int(row['mismatch_5']) + int(row['mismatch_3']))
                mismatch_counts[(row['id'], 'Total')][total_mismatches] += 1

        df = pd.DataFrame.from_dict(mismatch_counts, orient='index')
        df.sort_index(inplace=True)  # order A, B, Total
        df['WGE format'] = df.apply(lambda row: row.to_dict(), axis=1)
        return df

    @property
    def mismatch_df(self):
        return self._mismatch_df

    def save_mismatches(self, output_file):
        self.mismatch_df.to_csv(output_file, sep='\t')
