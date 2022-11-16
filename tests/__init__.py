import unittest
import os
import sys

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import src


def suite():
    return unittest.TestLoader().discover("satlib_app.tests", pattern="*.py")
