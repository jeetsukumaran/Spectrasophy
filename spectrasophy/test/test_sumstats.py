#! /usr/bin/env python

import os
import unittest
from spectrasophy import utility
from spectrasophy.test import TESTS_DATA_DIR


class SfsFixture(object):

    def __init__(self, filename_prefix):
        self.filename_prefix = filename_prefix
        self.dirpath = os.path.join(TESTS_DATA_DIR, "genomic")
        self.read_source_data()
        print(self.dna)
        print(self.snps)

    def read_source_data(self):
        self.dna = []
        self.snps = []
        for idx in range(0, 1):
            self.dna.append([])
            src_path = os.path.join(self.dirpath, "{}.data.dna.deme{}.txt".format(self.filename_prefix, idx))
            with utility.universal_open(src_path) as src:
                for row in src:
                    self.dna[idx].append(row.strip())
            self.snps.append([])
            src_path = os.path.join(self.dirpath, "{}.data.snps.deme{}.txt".format(self.filename_prefix, idx))
            with utility.universal_open(src_path) as src:
                for row in src:
                    self.snps[idx].append(row.strip())

_test_fixtures = [
        SfsFixture("testdataset01"),
        ]

# class SiteFrequencySpectrumCalcTestCase(unittest.TestCase):

#     def tes


if __name__ == "__main__":
    unittest.main()


