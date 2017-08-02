#! /usr/bin/env python

##############################################################################
## Copyright (c) 2017 Jeet Sukumaran.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
## IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
## THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
## AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import dendropy
import collections
from spectrasophy import model

class SpectrasophySummaryStatsCalculator(object):

    def __init__(self, **kwargs):
        self.output_prefix = kwargs.pop("output_prefix", "spectrasophy")
        self.is_unfolded_site_frequency_spectrum = kwargs.pop("is_unfolded_site_frequency_spectrum", False)
        self.is_calculate_single_population_sfs = kwargs.pop("is_calculate_single_population_sfs", False)
        self.is_calculate_joint_population_sfs = kwargs.pop("is_calculate_joint_population_sfs", True)
        self.stat_label_prefix = kwargs.pop("stat_label_prefix", "stat")
        self.supplemental_labels = kwargs.pop("supplemental_labels", None)
        locus_info = kwargs.pop("locus_info", None)
        if locus_info:
            self.model = model.SpectrasophyModel(params_d=None, locus_info=locus_info,)
        else:
            self.model = None
        if kwargs:
            raise Exception("Unrecognized configuration entries: {}".format(kwargs))

    def read_data(self, filepath, datatype, schema):
        if datatype == "dna":
            data = dendropy.DnaCharacterMatrix.get(path=filepath, schema=schema)
        elif datatype == "standard" or datatype == "snp":
            data = dendropy.StandardCharacterMatrix.get(path=filepath, schema=schema)
        return data

    def folded_joint_site_frequency_spectrum(self, d0, d1):
        # weirdly, FastsimCoal2 puts first deme second axis, i.e. columns,
        # while second deme gets put on rows
        demes = (d0, d1)
        jsfs = [[0 for i in range(len(d1)+1)] for j in range(len(d0)+1)]
        num_demes = 2
        nsites = None
        deme_site_columns = []
        for deme_idx in range(num_demes):
            deme_sites = zip(*(s.symbols_as_string() for s in demes[deme_idx].sequences()))
            if nsites is None:
                nsites = len(deme_sites)
            else:
                assert len(deme_sites) == nsites
            deme_site_columns.append(deme_sites)
        for site_idx in range(len(deme_site_columns[0])):
            deme_counters = []
            pooled_counter = collections.Counter()
            for deme_idx in range(num_demes):
                deme_counter = collections.Counter(deme_site_columns[deme_idx][site_idx])
                deme_counters.append(deme_counter)
                pooled_counter.update(deme_counter)
            if len(pooled_counter) == 1:
                jsfs[0][0] += 1
                continue
            majority_allele = pooled_counter.most_common(1)[0][0]
            for deme_idx in range(num_demes):
                del deme_counters[deme_idx][majority_allele]
            jsfs[sum(deme_counters[1].values())][sum(deme_counters[0].values())] += 1
        return jsfs
        # site_columns_list = [zip(*self.sequences())]
        # num_sites = len(site_columns[0])
        # for arg in args:
        #     c = zip(*arg.sequences())
        #     assert len(c) == num_sites#, "{} != {}".format(len(c), num_sites)
        #     site_columns_list.append(c)

        # for site_idx in range(num_sites):
        #     pooled_counter = collections.Counter()
        #     deme_site_counters = []
        #     for site_columns in site_columns_list:
        #         deme_site_counter = collections.Counter(site_columns[site_idx])
        #         pooled_counter.update(deme_site_counter)
        #         deme_site_counters.append(deme_site_counter)

