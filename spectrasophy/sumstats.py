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
        elif datatype == "snp":
            data = dendropy.StandardCharacterMatrix.get(path=filepath, schema=schema)
        return data

