from unittest.mock import patch
from io import StringIO

import pandas as pd
import numpy as np
from pyfakefs.fake_filesystem_unittest import TestCase

from scoring import Scoring


class TestScoring(TestCase):
    def setUp(self):
        self.setUpPyfakefs()
        file_contents = (
            'ipcress: 10:filter(unmasked) SMARCA4_exon24_1 '
            '300 A 48790792 1 A 48791074 2 single_A\n'
            'ipcress: 12:filter(unmasked) SMARCA4_exon24_1 '
            '225 B 132750362 2 A 132750569 2 revcomp\n'
            'ipcress: 19:filter(unmasked) SMARCA4_exon24_1 '
            '252 A 11027747 0 B 11027978 0 forward\n'
            'ipcress: 1:filter(unmasked) SMARCA4_exon24_3 '
            '272 A 171950952 2 B 171951204 1 forward\n'
            'ipcress: 19:filter(unmasked) SMARCA4_exon24_3 '
            '278 A 11027755 0 B 11028013 0 forward\n'
            '-- completed ipcress analysis'
        )
        self.fs.create_file('/test_input.txt', contents=file_contents)

        index = pd.MultiIndex.from_tuples([
            ('SMARCA4_exon24_1', 'A'),
            ('SMARCA4_exon24_1', 'B'),
            ('SMARCA4_exon24_1', 'Total'),
            ('SMARCA4_exon24_3', 'A'),
            ('SMARCA4_exon24_3', 'B'),
            ('SMARCA4_exon24_3', 'Total'),
        ], names=['Primer pair', None])
        self.df = pd.DataFrame({
            '0': [1, 1, 1, 1, 1, 1],
            '1': [1, 0, 0, 0, 1, 0],
            '2': [2, 1, 0, 1, 0, 0],
            '3': [0, 0, 1, 0, 0, 1],
            '4': [0, 0, 1, 0, 0, 0],
            'WGE format': [
                {'0': 1, '1': 1, '2': 2, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 1, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 0},
            ]
        }, index=index)

        index = pd.MultiIndex.from_tuples([
            ('SMARCA4_exon24_3', 'A'),
            ('SMARCA4_exon24_3', 'B'),
            ('SMARCA4_exon24_3', 'Total'),
            ('SMARCA4_exon24_1', 'A'),
            ('SMARCA4_exon24_1', 'B'),
            ('SMARCA4_exon24_1', 'Total'),
        ], names=['Primer pair', None])
        self.score_df = pd.DataFrame({
            '0': [1, 1, 1, 1, 1, 1],
            '1': [0, 1, 0, 1, 0, 0],
            '2': [1, 0, 0, 2, 1, 0],
            '3': [0, 0, 1, 0, 0, 1],
            '4': [0, 0, 0, 0, 0, 1],
            'WGE format': [
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 1, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 0},
                {'0': 1, '1': 1, '2': 2, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
            ], 'Score': [np.nan, np.nan, 100000, np.nan, np.nan, 110000]
        }, index=index)

    def test_mismatches_to_df_success(self):
        # arrange
        expected = self.df

        # act
        actual = Scoring.mismatches_to_df('/test_input.txt', 2)

        # assert
        pd.testing.assert_frame_equal(actual, expected)

    @patch('scoring.Scoring.mismatches_to_df')
    def check_score_df(self, mock_mismatches_to_df, check_like):
        # arrange
        mock_mismatches_to_df.return_value = self.df
        scores = [np.nan, np.nan, 110000, np.nan, np.nan, 100000]
        expected = self.score_df

        # act
        scoring = Scoring('/test_input.txt', 2)
        with patch.object(pd.DataFrame, 'apply', return_value=scores):
            scoring.add_scores_to_df()
        actual = scoring.mismatch_df

        # assert
        pd.testing.assert_frame_equal(actual, expected, check_like=check_like)

    def test_add_scores_to_df_adds_scores(self):
        self.check_score_df(check_like=True)

    def test_add_scores_to_df_orders_by_score(self):
        self.check_score_df(check_like=False)

    def test_score_mismatches_returns_score_for_total_row(self):
        # arrange
        mismatches = pd.Series({
            '0': 1, '1': 0, '2': 0, '3': 1, '4': 1,
            'WGE format': {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
        }, name=('SMARCA4_exon24_1', 'Total'))
        expected = 110000

        # act
        actual = Scoring.score_mismatches(mismatches)

        # assert
        self.assertEqual(actual, expected)

    def test_score_mismatches_returns_nan_for_non_total_row(self):
        # arrange
        mismatches = pd.Series({
            '0': 1, '1': 1, '2': 2, '3': 0, '4': 0,
            'WGE format': {'0': 0, '1': 1, '2': 2, '3': 0, '4': 0},
        }, name=('SMARCA4_exon24_1', 'A'))
        expected = np.nan

        # act
        actual = Scoring.score_mismatches(mismatches)

        # assert
        np.testing.assert_equal(actual, expected)

    @patch('builtins.print')
    def test_score_mismatches_prints_warning_if_no_match(self, mock_print):
        # arrange
        mismatches = pd.Series({
            '0': 0, '1': 0, '2': 0, '3': 1, '4': 1,
            'WGE format': {'0': 0, '1': 0, '2': 0, '3': 1, '4': 1},
        }, name=('SMARCA4_exon24_1', 'Total'))

        # act
        actual = Scoring.score_mismatches(mismatches)

        # assert
        mock_print.assert_called_with('No match found for SMARCA4_exon24_1')

    @patch('scoring.Scoring.mismatches_to_df')
    def test_save_mismatches_success(self, mock_mismatches_to_df):
        # arrange
        mock_mismatches_to_df.return_value = self.df
        fake_output_file = StringIO()
        expected = (
            "Primer pair\t\t0\t1\t2\t3\t4\tWGE format\n"
            "SMARCA4_exon24_1\tA\t1\t1\t2\t0\t0\t"
            "{'0': 1, '1': 1, '2': 2, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_1\tB\t1\t0\t1\t0\t0\t"
            "{'0': 1, '1': 0, '2': 1, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_1\tTotal\t1\t0\t0\t1\t1\t"
            "{'0': 1, '1': 0, '2': 0, '3': 1, '4': 1}\n"
            "SMARCA4_exon24_3\tA\t1\t0\t1\t0\t0\t"
            "{'0': 1, '1': 0, '2': 1, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_3\tB\t1\t1\t0\t0\t0\t"
            "{'0': 1, '1': 1, '2': 0, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_3\tTotal\t1\t0\t0\t1\t0\t"
            "{'0': 1, '1': 0, '2': 0, '3': 1, '4': 0}\n"
        )

        # act
        Scoring('/test_input.txt', 2).save_mismatches(fake_output_file)
        fake_output_file.seek(0)
        actual = fake_output_file.read()

        # assert
        self.assertEqual(actual, expected)
