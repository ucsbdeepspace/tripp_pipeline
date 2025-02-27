"""
Author: Benjamin Fogiel - bfogiel@ucsb.edu
Tests the read and write functions
"""

import os
import unittest
from click.testing import CliRunner
from astropy.io import fits
import tripp
import glob
import tempfile
import getpass
from datetime import datetime

class TestRead(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/science")
    nopath = os.path.join(os.path.dirname(__file__), "fixtures/gobbledygook")

    def setUp(self):
        self.output = [s for s in tripp.read(TestRead.path)]
        self.noutput = [s for s in tripp.read(TestRead.nopath)]

    def tearDown(self):
        for hdul in self.output:
            hdul.close()

    def test_length(self):
        self.assertEqual(len(self.output), 10, "Did not read the correct number of HDULs in directroy, fixtures/science")
        self.assertEqual(len(self.noutput), 0, "Did not read zero HDULs in from empty directory.")

    def test_type(self):
        for o in self.output:
            self.assertIsInstance(o, fits.HDUList, "Did not read type fits.HDUList")

    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(tripp._read_cmd, ['-d', 'filePath'])
        assert result.exit_code == 0

class TestWrite(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        ## Creates a unique dir to write to
        self.testDir = getpass.getuser() + datetime.now().strftime("%H%M%S")
        self.tmpdir = tempfile.TemporaryDirectory()
        self.dirName = os.path.join(self.tmpdir.name, self.testDir)
        os.mkdir(self.dirName)
        self.h = tripp.read(os.path.join(os.path.dirname(__file__), "fixtures/science"))
        self.output = tripp.write(self.h, os.path.join(os.path.dirname(__file__), self.dirName), "{number}.fits")
        self.paths = glob.glob("{}/*.fits*".format(os.path.join(os.path.dirname(__file__), self.dirName)))

    @classmethod
    def tearDownClass(self):
        for hdul in self.h:
            hdul.close(hdul)
        self.tmpdir.cleanup()

    def test_type(self):
        for o in self.output:
            self.assertIsInstance(o, fits.HDUList, "Did not write type fits.HDUList")

    def test_dir(self):
        self.assertEqual(len(self.paths), 10,f"Did not write the correct number of HDULs to directory {self.paths}")

    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(tripp._write_cmd, ['-d', 'filePath', '-f', 'format'])
        assert result.exit_code == 0

if __name__ == "__main__":
    unittest.main()
