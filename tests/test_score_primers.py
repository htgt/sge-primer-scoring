import argparse

from pyfakefs.fake_filesystem_unittest import TestCase

from score_primers import positive_int, non_empty_file, new_file_path


class TestScorePrimers(TestCase):
    def setUp(self):
        self.setUpPyfakefs()
        self.fs.create_file('non_empty_file.txt', contents='test')
        self.fs.create_dir('existing_dir')
        self.fs.create_file('empty_file.txt')

    def test_non_empty_file_non_empty_file_success(self):
        # arrange
        test_arg = 'non_empty_file.txt'
        expected = 'non_empty_file.txt'

        # act
        actual = non_empty_file(test_arg)

        # assert
        self.assertEqual(actual, expected)

    def test_non_empty_file_non_existing_file_fail(self):
        # arrange
        test_arg = 'non_existing_file.txt'
        expected = "File does not exist: 'non_existing_file.txt'"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as cm:
            non_empty_file(test_arg)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_non_empty_file_dir_fail(self):
        # arrange
        test_arg = 'existing_dir'
        expected = "File does not exist: 'existing_dir'"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as cm:
            non_empty_file(test_arg)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_non_empty_file_empty_file_fail(self):
        # arrange
        test_arg = 'empty_file.txt'
        expected = "File is empty: 'empty_file.txt'"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as cm:
            non_empty_file(test_arg)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_positive_int_positive_arg_success(self):
        # arrange
        test_arg = '1'
        expected = 1

        # act
        actual = positive_int(test_arg)

        # assert
        self.assertEqual(actual, expected)

    def test_positive_int_negative_arg_fail(self):
        # arrange
        test_arg = '-1'
        expected = 'Mismatch number cannot be negative'

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as cm:
            positive_int(test_arg)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_new_file_path_new_file_path_success(self):
        # arrange
        test_arg = 'new_file.txt'
        expected = 'new_file.txt'

        # act
        actual = new_file_path(test_arg)

        # assert
        self.assertEqual(actual, expected)

    def test_new_file_path_new_dir_path_fail(self):
        # arrange
        test_arg = 'new_dir/'
        expected = "Directory provided rather than file path: 'new_dir/'"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as cm:
            new_file_path(test_arg)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_new_file_path_existing_dir_path_fail(self):
        # arrange
        test_arg = 'existing_dir'
        expected = "Directory provided rather than file path: 'existing_dir'"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as cm:
            new_file_path(test_arg)

        # assert
        self.assertEqual(str(cm.exception), expected)

    def test_new_file_path_existing_file_path_fail(self):
        # arrange
        test_arg = 'empty_file.txt'
        expected = "File already exists: 'empty_file.txt'"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as cm:
            new_file_path(test_arg)

        # assert
        self.assertEqual(str(cm.exception), expected)
