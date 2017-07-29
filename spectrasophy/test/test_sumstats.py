#! /usr/bin/env python

import os
import unittest
from spectrasophy import utility
from spectrasophy import sumstats
from spectrasophy.test import TESTS_DATA_DIR


class SfsTestDataSet(object):

    def __init__(self, filename_prefix):
        self.filename_prefix = filename_prefix
        self.dirpath = os.path.join(TESTS_DATA_DIR, "genomic")
        self.read_source_data()

    def read_source_data(self):
        self.dna = []
        self.snps = []
        self.dna_fasta_filepaths = []
        self.snps_fasta_filepaths = []
        for idx in range(2):
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
            self.dna_fasta_filepaths.append(os.path.join(self.dirpath, "testdataset01.data.dna.deme{}.fasta".format(idx)))
            self.snps_fasta_filepaths.append(os.path.join(self.dirpath, "testdataset01.data.snps.deme{}.fasta".format(idx)))

_TEST_DATASETS = [
        SfsTestDataSet("testdataset01"),
        ]

class DataReadingTestCases(unittest.TestCase):

    def test_read_dna(self):
        sc = sumstats.SpectrasophySummaryStatsCalculator()
        for test_ds in _TEST_DATASETS:
            for idx in range(2):
                data = sc.read_data(
                        filepath=test_ds.dna_fasta_filepaths[idx],
                        datatype="dna",
                        schema="fasta")
                obs_sequences = [s.symbols_as_string() for s in data.sequences()]
                self.assertEqual(obs_sequences, test_ds.dna[idx])

    def test_read_snps(self):
        sc = sumstats.SpectrasophySummaryStatsCalculator()
        for test_ds in _TEST_DATASETS:
            for idx in range(2):
                data = sc.read_data(
                        filepath=test_ds.snps_fasta_filepaths[idx],
                        datatype="standard",
                        schema="fasta")
                obs_sequences = [s.symbols_as_string() for s in data.sequences()]
                self.assertEqual(obs_sequences, test_ds.snps[idx])

class SiteFrequencySpectrumCalcTestCase(unittest.TestCase):

    def test_unfolded_sfs(self):
        print("ok")
        sc = sumstats.SpectrasophySummaryStatsCalculator()
        for test_ds in _TEST_DATASETS:
            for datatype, source_filepaths in (
                    ("dna", test_ds.dna_fasta_filepaths,),
                    ("standard", test_ds.snps_fasta_filepaths),
                    ):
                for idx in range(2):
                    data = sc.read_data(
                            filepath=source_filepaths[idx],
                            datatype=datatype,
                            schema="fasta")
                    print("> {}: {}: {}".format(datatype, idx, data.folded_site_frequency_spectrum()))
                    print("--")

if __name__ == "__main__":
    unittest.main()


