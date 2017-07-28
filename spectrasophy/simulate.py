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

import collections
import random
import sys
import os
import time
try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing
import traceback

import spectrasophy
from spectrasophy import utility
from spectrasophy import model
from spectrasophy import fsc2

class SimulationWorker(multiprocessing.Process):

    def __init__(self,
            name,
            model,
            work_queue,
            results_queue,
            fsc2_path,
            working_directory,
            run_logger,
            logging_frequency,
            messenger_lock,
            random_seed,
            is_calculate_single_population_sfs,
            is_calculate_joint_population_sfs,
            is_unfolded_site_frequency_spectrum,
            stat_label_prefix,
            is_include_model_id_field,
            supplemental_labels,
            debug_mode,
            ):
        multiprocessing.Process.__init__(self, name=name)
        self.fsc2_handler = fsc2.Fsc2Handler(
                name=name,
                fsc2_path=fsc2_path,
                working_directory=working_directory,
                is_calculate_single_population_sfs=is_calculate_single_population_sfs,
                is_calculate_joint_population_sfs=is_calculate_joint_population_sfs,
                is_unfolded_site_frequency_spectrum=is_unfolded_site_frequency_spectrum)
        self.model = model
        self.rng = random.Random(random_seed)
        self.work_queue = work_queue
        self.results_queue = results_queue
        self.run_logger = run_logger
        self.logging_frequency = logging_frequency
        self.messenger_lock = messenger_lock
        self.is_unfolded_site_frequency_spectrum = is_unfolded_site_frequency_spectrum
        self.stat_label_prefix = stat_label_prefix
        self.is_include_model_id_field = is_include_model_id_field
        self.supplemental_labels = supplemental_labels
        self.is_debug_mode = debug_mode
        self.kill_received = False
        self.num_tasks_received = 0
        self.num_tasks_completed = 0

    def send_worker_message(self, msg, level):
        if self.run_logger is None:
            return
        # if self.run_logger.messaging_level > level or self.messenger.silent:
        #     return
        msg = "{}: {}".format(self.name, msg)
        self.messenger_lock.acquire()
        try:
            self.run_logger.log(msg, level=level)
        finally:
            self.messenger_lock.release()

    def send_worker_critical(self, msg):
        self.send_worker_message(msg, utility.RunLogger.CRITICAL_MESSAGING_LEVEL)

    def send_worker_debug(self, msg):
        self.send_worker_message(msg, utility.RunLogger.DEBUG_MESSAGING_LEVEL)

    def send_worker_info(self, msg):
        self.send_worker_message(msg, utility.RunLogger.INFO_MESSAGING_LEVEL)

    def send_worker_warning(self, msg):
        self.send_worker_message(msg, utility.RunLogger.WARNING_MESSAGING_LEVEL)

    def send_worker_error(self, msg):
        self.send_worker_message(msg, utility.RunLogger.ERROR_MESSAGING_LEVEL)

    def run(self):
        result = None
        while not self.kill_received:
            try:
                rep_idx = self.work_queue.get_nowait()
            except queue.Empty:
                break
            self.num_tasks_received += 1
            # self.send_worker_critical("Received task: '{task_name}'".format(
            #     task_count=self.num_tasks_received,
            #     task_name=rep_idx))
            # rng = random.Random(random_seed)
            try:
                result = self.simulate()
            except (KeyboardInterrupt, Exception) as e:
                # traceback.print_exc()
                e.worker_name = self.name
                e.traceback_exc = traceback.format_exc()
                self.results_queue.put(e)
                break
            if self.kill_received:
                break
            self.results_queue.put(result)
            self.num_tasks_completed += 1
            # self.send_info("Completed task {task_count}: '{task_name}'".format(
            if rep_idx and self.logging_frequency and rep_idx % self.logging_frequency == 0:
                self.run_logger.info("Completed replicate {task_name}".format(
                    task_count=self.num_tasks_received,
                    task_name=rep_idx))
        if self.kill_received:
            self.send_worker_warning("Terminating in response to kill request")

    def simulate(self):
        results_d = collections.OrderedDict()
        if self.is_include_model_id_field:
            results_d["model.id"] = None
        if self.supplemental_labels:
            for key in self.supplemental_labels:
                results_d[key] = self.supplemental_labels[key]
        params, fsc2_run_configurations = self.model.sample_parameter_values_from_prior(rng=self.rng)
        results_d.update(params)
        for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
            for locus_definition in lineage_pair.locus_definitions:
                self.fsc2_handler.run(
                        field_name_prefix="{}.{}.{}".format(
                                self.stat_label_prefix,
                                lineage_pair.label,
                                locus_definition.locus_label),
                        fsc2_config_d=fsc2_run_configurations[locus_definition],
                        random_seed=self.rng.randint(1, 1E6),
                        results_d=results_d,
                        )
        if self.is_include_model_id_field:
            results_d["model.id"] = results_d["param.divTimeModel"]
        return results_d

class SpectrasophySimulator(object):

    def __init__(self,
            config_d,
            num_processes=None,
            logging_frequency=1000,
            package_id=None,
            is_verbose_setup=True):
        # configure
        if package_id is None:
            self.package_id = spectrasophy.package_id()
        else:
            self.package_id = package_id
        self.elapsed_time = 0.0 # need to be here for logging
        config_d = dict(config_d) # make copy so we can pop items
        self.is_verbose_setup = is_verbose_setup
        self.configure_simulator(config_d)
        self.num_cpus = multiprocessing.cpu_count()
        if num_processes is None or num_processes <= 0:
            self.num_processes = num_cpus
        elif num_processes == 1 and self.num_cpus > 1 and self.is_verbose_setup:
            self.run_logger.info(
                    ("Multiple processors ({num_cpus}) available:"
                    " consider using the '-M' or '-m' options to"
                    " parallelize processing of trees"
                    ).format(num_cpus=self.num_cpus))
            self.num_processes = 1
        else:
            self.num_processes = num_processes
        if self.is_verbose_setup:
            self.run_logger.info("Will run up to {} processes in parallel".format(self.num_processes))
            self.run_logger.info("{} lineage pairs in analysis:".format(self.model.num_lineage_pairs))
            for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
                self.run_logger.info("  - '{}': {:>2d} loci (Samples: {})".format(
                        lineage_pair.label,
                        len(lineage_pair.locus_definitions),
                        ", ".join("{}/{}".format(locus.num_genes_deme0, locus.num_genes_deme1) for locus in lineage_pair.locus_definitions),
                        ))
        self.worker_class = SimulationWorker

    def configure_simulator(self, config_d, verbose=True):
        self.title = config_d.pop("title", "spectrasophy-{}-{}".format(time.strftime("%Y%m%d%H%M%S"), id(self)))
        self.output_prefix = config_d.pop("output_prefix", self.title)
        self.working_directory = config_d.pop("working_directory", self.title)
        self.run_logger = config_d.pop("run_logger", None)
        if self.run_logger is None:
            self.run_logger = utility.RunLogger(
                    name="spectrasophy-simulate",
                    stderr_logging_level=config_d.pop("standard_error_logging_level", "info"),
                    log_to_file=config_d.pop("log_to_file", True),
                    log_to_stderr=config_d.pop("log_to_stderr", True),
                    log_path=self.output_prefix + ".log",
                    file_logging_level=config_d.pop("file_logging_level", "info"),
                    )
        self.run_logger.system = self
        self.logging_frequency = config_d.pop("logging_frequency", 1000)
        if self.is_verbose_setup:
            self.run_logger.info("Running: {}".format(self.package_id))
            self.run_logger.info("Configuring simulation: '{}'".format(self.title))
            self.run_logger.info("Output directory: '{}'".format(os.path.dirname(os.path.abspath(self.output_prefix))))
            self.run_logger.info("Output filename prefix: '{}'".format(os.path.basename(self.output_prefix)))
            self.run_logger.info("Working directory: '{}'".format(self.working_directory))
        self.fsc2_path = config_d.pop("fsc2_path", "fsc25")
        if self.is_verbose_setup:
            self.run_logger.info("FastSimCoal2 path: '{}'".format(self.fsc2_path))
        self.rng = config_d.pop("rng", None)
        if self.rng is None:
            self.random_seed = config_d.pop("random_seed", None)
            if self.random_seed is None:
                self.random_seed = random.randint(0, sys.maxsize)
            if self.is_verbose_setup:
                self.run_logger.info("Initializing with random seed: {}".format(self.random_seed))
            self.rng = random.Random(self.random_seed)
        else:
            if "random_seed" in config_d:
                raise TypeError("Cannot specify both 'rng' and 'random_seed'")
            if self.is_verbose_setup:
                self.run_logger.info("Using existing random number generator")
        self.is_debug_mode = config_d.pop("debug_mode", False)
        if self.is_verbose_setup and self.is_debug_mode:
            self.run_logger.info("Running in DEBUG mode")
        self.site_frequency_spectrum_type = config_d.pop("site_frequency_spectrum_type", "unfolded").lower()
        self.is_unfolded_site_frequency_spectrum = config_d.pop("is_unfolded_site_frequency_spectrum", False)
        self.is_calculate_single_population_sfs = config_d.pop("is_calculate_single_population_sfs", False)
        self.is_calculate_joint_population_sfs = config_d.pop("is_calculate_joint_population_sfs", True)
        if not self.is_calculate_single_population_sfs and not self.is_calculate_joint_population_sfs:
            raise ValueError("Neither single-population nor joint site frequency spectrum will be calculated!")
        self.stat_label_prefix = config_d.pop("stat_label_prefix", "stat")
        self.supplemental_labels = config_d.pop("supplemental_labels", None)
        self.is_include_model_id_field = config_d.pop("is_include_model_id_field", False)
        if "params" not in config_d:
            raise ValueError("Missing 'params' entry in configuration")
        params_d = config_d.pop("params")
        if "locus_info" not in config_d:
            raise ValueError("Missing 'locus_info' entry in configuration")
        locus_info = config_d.pop("locus_info")
        self.model = model.SpectrasophyModel(params_d=params_d, locus_info=locus_info,)
        if config_d:
            raise Exception("Unrecognized configuration entries: {}".format(config_d))

    def execute(self,
            nreps,
            results_csv_writer=None,
            results_store=None,
            is_write_header=True,
            ):
        # load up queue
        self.run_logger.info("Creating work queue")
        work_queue = multiprocessing.Queue()
        for rep_idx in range(nreps):
            work_queue.put( rep_idx )
        time.sleep(0.1) # to avoid: 'IOError: [Errno 32] Broken pipe'; https://stackoverflow.com/questions/36359528/broken-pipe-error-with-multiprocessing-queue
        self.run_logger.info("Launching {} worker processes".format(self.num_processes))
        results_queue = multiprocessing.Queue()
        messenger_lock = multiprocessing.Lock()
        workers = []
        for pidx in range(self.num_processes):
            worker = self.worker_class(
                    # name=str(pidx+1),
                    name="{}-{}".format(self.title, pidx+1),
                    model=self.model,
                    work_queue=work_queue,
                    results_queue=results_queue,
                    fsc2_path=self.fsc2_path,
                    working_directory=self.working_directory,
                    run_logger=self.run_logger,
                    logging_frequency=self.logging_frequency,
                    messenger_lock=messenger_lock,
                    random_seed=self.rng.randint(1, sys.maxsize),
                    is_calculate_single_population_sfs=self.is_calculate_single_population_sfs,
                    is_calculate_joint_population_sfs=self.is_calculate_joint_population_sfs,
                    is_unfolded_site_frequency_spectrum=self.is_unfolded_site_frequency_spectrum,
                    stat_label_prefix=self.stat_label_prefix,
                    is_include_model_id_field=self.is_include_model_id_field,
                    supplemental_labels=self.supplemental_labels,
                    debug_mode=self.is_debug_mode,
                    )
            worker.start()
            workers.append(worker)

        # collate results
        result_count = 0
        try:
            while result_count < nreps:
                result = results_queue.get()
                if isinstance(result, KeyboardInterrupt):
                    raise result
                elif isinstance(result, Exception):
                    self.run_logger.error("Exception raised in worker process '{}'"
                                          "\n>>>\n{}<<<\n".format(
                                              result.worker_name,
                                              result.traceback_exc))
                    raise result
                if results_store is not None:
                    results_store.append(result)
                if results_csv_writer is not None:
                    if result_count == 0 and is_write_header:
                        results_csv_writer.fieldnames = result.keys()
                        results_csv_writer.writeheader()
                    results_csv_writer.writerow(result)
                # self.run_logger.info("Recovered results from worker process '{}'".format(result.worker_name))
                result_count += 1
                # self.info_message("Recovered results from {} of {} worker processes".format(result_count, self.num_processes))
        except (Exception, KeyboardInterrupt) as e:
            for worker in workers:
                worker.terminate()
            raise
        self.run_logger.info("All {} worker processes terminated".format(self.num_processes))
        return results_store

