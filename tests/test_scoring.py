from unittest.mock import patch
from os import path

import pandas as pd
import numpy as np
from pyfakefs.fake_filesystem_unittest import TestCase

from scoring import Scoring, ScoringError


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
            'ipcress: 13:filter(unmasked) BRCA1_exon1_1 '
            '207 A 32315485 0 B 32315669 0 forward\n'
            '-- completed ipcress analysis\n'
        )
        self.fs.create_file('/ipcress.txt', contents=file_contents)
        file_contents = (
            'SMARCA4_exon24_1,Targeton_1\n'
            'SMARCA4_exon24_3,Targeton_1\n'
            'BRCA1_exon1_1,Targeton_2'
        )
        self.fs.create_file('/targetons.csv', contents=file_contents)

        index = pd.MultiIndex.from_tuples([
            ('BRCA1_exon1_1', 'A'),
            ('BRCA1_exon1_1', 'B'),
            ('BRCA1_exon1_1', 'Total'),
            ('SMARCA4_exon24_1', 'A'),
            ('SMARCA4_exon24_1', 'B'),
            ('SMARCA4_exon24_1', 'Total'),
            ('SMARCA4_exon24_3', 'A'),
            ('SMARCA4_exon24_3', 'B'),
            ('SMARCA4_exon24_3', 'Total'),
        ], names=['Primer pair', 'A/B/Total'])
        self.df = pd.DataFrame({
            '0': [1, 1, 1, 1, 1, 1, 1, 1, 1],
            '1': [0, 0, 0, 1, 0, 0, 0, 1, 0],
            '2': [0, 0, 0, 2, 1, 0, 1, 0, 0],
            '3': [0, 0, 0, 0, 0, 1, 0, 0, 1],
            '4': [0, 0, 0, 0, 0, 1, 0, 0, 0],
            'WGE format': [
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 1, '2': 2, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 1, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 0},
            ]
        }, index=index)

        index = pd.MultiIndex.from_tuples([
            ('Targeton_1', 'SMARCA4_exon24_1', 'A'),
            ('Targeton_1', 'SMARCA4_exon24_1', 'B'),
            ('Targeton_1', 'SMARCA4_exon24_1', 'Total'),
            ('Targeton_1', 'SMARCA4_exon24_3', 'A'),
            ('Targeton_1', 'SMARCA4_exon24_3', 'B'),
            ('Targeton_1', 'SMARCA4_exon24_3', 'Total'),
            ('Targeton_2', 'BRCA1_exon1_1', 'A'),
            ('Targeton_2', 'BRCA1_exon1_1', 'B'),
            ('Targeton_2', 'BRCA1_exon1_1', 'Total'),
        ], names=['Targeton', 'Primer pair', 'A/B/Total'])
        self.targeton_df = pd.DataFrame({
            '0': [1, 1, 1, 1, 1, 1, 1, 1, 1],
            '1': [1, 0, 0, 0, 1, 0, 0, 0, 0],
            '2': [2, 1, 0, 1, 0, 0, 0, 0, 0],
            '3': [0, 0, 1, 0, 0, 1, 0, 0, 0],
            '4': [0, 0, 1, 0, 0, 0, 0, 0, 0],
            'WGE format': [
                {'0': 1, '1': 1, '2': 2, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 1, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
            ]
        }, index=index)

    def test_mismatches_to_df_no_targeton_csv_success(self):
        # arrange
        expected = self.df

        # act
        actual = Scoring.mismatches_to_df('/ipcress.txt', 2)

        # assert
        pd.testing.assert_frame_equal(actual, expected)

    def test_mismatches_to_df_targeton_csv_success(self):
        # arrange
        expected = self.targeton_df

        # act
        actual = Scoring.mismatches_to_df('/ipcress.txt', 2, '/targetons.csv')

        # assert
        pd.testing.assert_frame_equal(actual, expected)

    def test_mismatches_to_df_invalid_ipcress_file_fail(self):
        # arrange
        self.fs.create_file('/invalid_input.txt', contents='invalid')
        expected = "Invalid ipcress file: '/invalid_input.txt'"

        # act
        with self.assertRaises(ScoringError) as cm:
            Scoring.mismatches_to_df('/invalid_input.txt', 2)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_mismatches_to_df_low_mismatch_fail(self):
        # arrange
        expected = "Mismatch number too low for ipcress file: '1'"

        # act
        with self.assertRaises(ScoringError) as cm:
            Scoring.mismatches_to_df('/ipcress.txt', 1)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_mismatches_to_df_no_ipcress_data_fail(self):
        # arrange
        file_contents = '-- completed ipcress analysis\n'
        self.fs.create_file('/empty_input.txt', contents=file_contents)
        expected = "No data in ipcress file: '/empty_input.txt'"

        # act
        with self.assertRaises(ScoringError) as cm:
            Scoring.mismatches_to_df('/empty_input.txt', 2)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_mismatches_to_df_invalid_targeton_csv_fail(self):
        # arrange
        self.fs.create_file('/invalid.csv', contents='invalid')
        expected = "Invalid targeton csv: '/invalid.csv'"

        # act
        with self.assertRaises(ScoringError) as cm:
            Scoring.mismatches_to_df('/ipcress.txt', 2, '/invalid.csv')

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_mismatches_to_df_conflicting_targeton_csv_entry_fail(self):
        # arrange
        file_contents = (
            'SMARCA4_exon24_1,Targeton_1\n'
            'SMARCA4_exon24_1,Targeton_2\n'
        )
        self.fs.create_file('/duplicates.csv', contents=file_contents)
        expected = ("Conflicting entries in targeton csv "
                    "for SMARCA4_exon24_1: '/duplicates.csv'")

        # act
        with self.assertRaises(ScoringError) as cm:
            Scoring.mismatches_to_df('/ipcress.txt', 2, '/duplicates.csv')

        # assert
        self.assertEqual(str(cm.exception), expected)

    @patch('scoring.Scoring.mismatches_to_df')
    def check_score_df(self, mock_mismatches_to_df, check_like):
        # arrange
        mock_mismatches_to_df.return_value = self.df
        apply_output = [
            np.nan, np.nan, 0, np.nan, np.nan,
            110000, np.nan, np.nan, 100000
        ]
        index = pd.MultiIndex.from_tuples([
            ('BRCA1_exon1_1', 'A'),
            ('BRCA1_exon1_1', 'B'),
            ('BRCA1_exon1_1', 'Total'),
            ('SMARCA4_exon24_3', 'A'),
            ('SMARCA4_exon24_3', 'B'),
            ('SMARCA4_exon24_3', 'Total'),
            ('SMARCA4_exon24_1', 'A'),
            ('SMARCA4_exon24_1', 'B'),
            ('SMARCA4_exon24_1', 'Total'),
        ], names=['Primer pair', 'A/B/Total'])
        scores = [
            np.nan, np.nan, 0, np.nan, np.nan,
            100000, np.nan, np.nan, 110000
        ]
        expected = pd.DataFrame({
            '0': [1, 1, 1, 1, 1, 1, 1, 1, 1],
            '1': [0, 0, 0, 0, 1, 0, 1, 0, 0],
            '2': [0, 0, 0, 1, 0, 0, 2, 1, 0],
            '3': [0, 0, 0, 0, 0, 1, 0, 0, 1],
            '4': [0, 0, 0, 0, 0, 0, 0, 0, 1],
            'WGE format': [
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 1, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 0},
                {'0': 1, '1': 1, '2': 2, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
            ], 'Score': scores}, index=index)

        # act
        scoring = Scoring('/ipcress.txt', 2)
        with patch.object(pd.DataFrame, 'apply', return_value=apply_output):
            scoring.add_scores_to_df()
        actual = scoring.mismatch_df

        # assert
        pd.testing.assert_frame_equal(actual, expected, check_like=check_like)

    def test_add_scores_to_df_adds_scores(self):
        self.check_score_df(check_like=True)

    def test_add_scores_to_df_orders_by_score(self):
        self.check_score_df(check_like=False)

    @patch('scoring.Scoring.mismatches_to_df')
    def test_add_scores_to_df_sorts_by_targeton(self, mock_mismatches_to_df):
        # arrange
        mock_mismatches_to_df.return_value = self.targeton_df
        apply_output = [
            np.nan, np.nan, 110000, np.nan,
            np.nan, 100000, np.nan, np.nan, 0
        ]
        index = pd.MultiIndex.from_tuples([
            ('Targeton_1', 'SMARCA4_exon24_3', 'A'),
            ('Targeton_1', 'SMARCA4_exon24_3', 'B'),
            ('Targeton_1', 'SMARCA4_exon24_3', 'Total'),
            ('Targeton_1', 'SMARCA4_exon24_1', 'A'),
            ('Targeton_1', 'SMARCA4_exon24_1', 'B'),
            ('Targeton_1', 'SMARCA4_exon24_1', 'Total'),
            ('Targeton_2', 'BRCA1_exon1_1', 'A'),
            ('Targeton_2', 'BRCA1_exon1_1', 'B'),
            ('Targeton_2', 'BRCA1_exon1_1', 'Total'),
        ], names=['Targeton', 'Primer pair', 'A/B/Total'])
        scores = [
            np.nan, np.nan, 100000, np.nan,
            np.nan, 110000, np.nan, np.nan, 0
        ]
        expected = pd.DataFrame({
            '0': [1, 1, 1, 1, 1, 1, 1, 1, 1],
            '1': [0, 1, 0, 1, 0, 0, 0, 0, 0],
            '2': [1, 0, 0, 2, 1, 0, 0, 0, 0],
            '3': [0, 0, 1, 0, 0, 1, 0, 0, 0],
            '4': [0, 0, 0, 0, 0, 1, 0, 0, 0],
            'WGE format': [
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 1, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 0},
                {'0': 1, '1': 1, '2': 2, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
                {'0': 1, '1': 0, '2': 0, '3': 0, '4': 0},
            ], 'Score': scores}, index=index)

        # act
        scoring = Scoring('/ipcress.txt', 2, '/targetons.csv')
        with patch.object(pd.DataFrame, 'apply', return_value=apply_output):
            scoring.add_scores_to_df()
        actual = scoring.mismatch_df

        # assert
        pd.testing.assert_frame_equal(actual, expected)

    def test_score_mismatches_returns_score_for_total_row_no_targeton(self):
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

    def test_score_mismatches_returns_score_for_total_row_with_targeton(self):
        # arrange
        mismatches = pd.Series({
            '0': 1, '1': 0, '2': 0, '3': 1, '4': 1,
            'WGE format': {'0': 1, '1': 0, '2': 0, '3': 1, '4': 1},
        }, name=('Targeton_1', 'SMARCA4_exon24_1', 'Total'))
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

    def test_score_mismatches_no_on_target_hit_fail_no_targeton(self):
        # arrange
        mismatches = pd.Series({
            '0': 0, '1': 0, '2': 0, '3': 1, '4': 1,
            'WGE format': {'0': 0, '1': 0, '2': 0, '3': 1, '4': 1},
        }, name=('SMARCA4_exon24_1', 'Total'))
        expected = 'No on-target hit found for SMARCA4_exon24_1'

        # act
        with self.assertRaises(ScoringError) as cm:
            Scoring.score_mismatches(mismatches)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_score_mismatches_no_on_target_hit_fail_with_targeton(self):
        # arrange
        mismatches = pd.Series({
            '0': 0, '1': 0, '2': 0, '3': 1, '4': 1,
            'WGE format': {'0': 0, '1': 0, '2': 0, '3': 1, '4': 1},
        }, name=('Targeton_1', 'SMARCA4_exon24_1', 'Total'))
        expected = 'No on-target hit found for SMARCA4_exon24_1'

        # act
        with self.assertRaises(ScoringError) as cm:
            Scoring.score_mismatches(mismatches)

        # assert
        self.assertEqual(str(cm.exception), expected)

    @patch('scoring.Scoring.mismatches_to_df')
    def test_save_mismatches_creates_file(self, mock_mismatches_to_df):
        # arrange
        mock_mismatches_to_df.return_value = self.df

        # act
        Scoring('/ipcress.txt', 2).save_mismatches('output.tsv')

        # assert
        self.assertTrue(path.exists('/output.tsv'))

    @patch('scoring.Scoring.mismatches_to_df')
    def test_save_mismatches_creates_parent_dir(self, mock_mismatches_to_df):
        # arrange
        mock_mismatches_to_df.return_value = self.df

        # act
        Scoring('/ipcress.txt', 2).save_mismatches('/test/output.tsv')

        # assert
        self.assertTrue(path.exists('/test/output.tsv'))

    @patch('scoring.Scoring.mismatches_to_df')
    def test_save_mismatches_correct_file_content(self, mock_mismatches_to_df):
        # arrange
        mock_mismatches_to_df.return_value = self.df
        expected = (
            "Primer pair\tA/B/Total\t0\t1\t2\t3\t4\tWGE format\n"
            "BRCA1_exon1_1\tA\t1\t0\t0\t0\t0\t"
            "{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0}\n"
            "BRCA1_exon1_1\tB\t1\t0\t0\t0\t0\t"
            "{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0}\n"
            "BRCA1_exon1_1\tTotal\t1\t0\t0\t0\t0\t"
            "{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0}\n"
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
        Scoring('/ipcress.txt', 2).save_mismatches('output.tsv')
        with open('/output.tsv') as f:
            actual = f.read()

        # assert
        self.assertEqual(actual, expected)
