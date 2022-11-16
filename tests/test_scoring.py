import unittest
from unittest.mock import patch
from io import StringIO

import pandas as pd
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
            'ipcress: 1:filter(unmasked) SMARCA4_exon24_3 '
            '272 A 171950952 2 B 171951204 1 forward\n'
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
        ])
        self.df = pd.DataFrame({
            '0': [0, 0, 0, 0, 0, 0],
            '1': [1, 0, 0, 0, 1, 0],
            '2': [2, 1, 0, 1, 0, 0],
            '3': [0, 0, 1, 0, 0, 1],
            '4': [0, 0, 1, 0, 0, 0],
            'WGE format': [
                {'0': 0, '1': 1, '2': 2, '3': 0, '4': 0},
                {'0': 0, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 0, '1': 0, '2': 0, '3': 1, '4': 1},
                {'0': 0, '1': 0, '2': 1, '3': 0, '4': 0},
                {'0': 0, '1': 1, '2': 0, '3': 0, '4': 0},
                {'0': 0, '1': 0, '2': 0, '3': 1, '4': 0},
            ]
        }, index=index)

    def test_mismatches_to_df_success(self):
        # arrange
        expected = self.df

        # act
        actual = Scoring.mismatches_to_df('/test_input.txt', 2)

        # assert
        pd.testing.assert_frame_equal(actual, expected)

    @patch('scoring.Scoring.mismatches_to_df')
    def test_save_mismatches_success(self, mock_mismatches_to_df):
        # arrange
        mock_mismatches_to_df.return_value = self.df
        fake_output_file = StringIO()
        expected = (
            "\t\t0\t1\t2\t3\t4\tWGE format\n"
            "SMARCA4_exon24_1\tA\t0\t1\t2\t0\t0\t"
            "{'0': 0, '1': 1, '2': 2, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_1\tB\t0\t0\t1\t0\t0\t"
            "{'0': 0, '1': 0, '2': 1, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_1\tTotal\t0\t0\t0\t1\t1\t"
            "{'0': 0, '1': 0, '2': 0, '3': 1, '4': 1}\n"
            "SMARCA4_exon24_3\tA\t0\t0\t1\t0\t0\t"
            "{'0': 0, '1': 0, '2': 1, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_3\tB\t0\t1\t0\t0\t0\t"
            "{'0': 0, '1': 1, '2': 0, '3': 0, '4': 0}\n"
            "SMARCA4_exon24_3\tTotal\t0\t0\t0\t1\t0\t"
            "{'0': 0, '1': 0, '2': 0, '3': 1, '4': 0}\n"
        )

        # act
        Scoring('/test_input.txt', 2).save_mismatches(fake_output_file)
        fake_output_file.seek(0)
        actual = fake_output_file.read()

        # assert
        self.assertEqual(actual, expected)
