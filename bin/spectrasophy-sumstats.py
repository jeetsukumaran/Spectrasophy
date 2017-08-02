#! /usr/bin/env python

import os
import sys
import argparse
import traceback
import time
import spectrasophy
from spectrasophy import sumstats
from spectrasophy import utility

def main():
    parser = argparse.ArgumentParser()
    package_id = spectrasophy.package_id()
    parser.add_argument("--version", action="version", version=package_id)

    simulator_options = parser.add_argument_group("Configuration")
    simulator_options.add_argument("configuration_filepath",
            metavar="CONFIGURATION-FILE",
            help="Path to the configuration file listing the data.")
    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument('-o', '--output-name-prefix',
            action='store',
            dest='output_name_prefix',
            type=str,
            default=None,
            metavar='NAME-PREFIX',
            help="Prefix for output filenames (default: same as configuration filename stem).")
    output_options.add_argument('-O', '--output-directory',
            action='store',
            dest='output_directory',
            type=str,
            default=None,
            metavar='DIRECTORY',
            help="Directory for output files (default: current working directory).")
    output_options.add_argument(
            "-U",
            "--unfolded-site-frequency-spectrum",
            "--derived-site-frequency-spectrum",
            action="store_true",
            default=False,
            help="Calculate the unfolded or derived site frequency spectrum."
            " Otherwise, defaults to the folded or minor site frequency"
            " spectrum."
            )
    output_options.add_argument(
            "--calculate-single-population-site-frequency-spectrum",
            action="store_true",
            default=False,
            help="Calculate the single (within) population site frequency"
            " spectrum in addition to the joint."
            )
    output_options.add_argument("-l", "--labels",
            action="append",
            help="Addition field/value pairs to add to the output (in format <FIELD-NAME>:value;)")
    output_options.add_argument('--field-delimiter',
            type=str,
            default='\t',
            help="Delimiter string separating fields in output (default: <TAB>').")
    output_options.add_argument('--summary-stats-label-prefix',
            type=str,
            default='stat',
            metavar='PREFIX',
            help="Prefix for summar statistic field labels (default: '%(default)s').")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append instead of overwriting output file(s).")
    output_options.add_argument( "--no-write-header",
            action="store_true",
            default=False,
            help="Do not writer header row.")

    args = parser.parse_args()

    config_d = {}
    utility.parse_legacy_configuration(
            filepath=args.configuration_filepath,
            config_d=config_d)
    config_d["output_prefix"] = utility.output_prefix(
            primary_source_filepath=args.configuration_filepath,
            output_name_prefix=args.output_name_prefix,
            output_directory=args.output_directory)
    config_d["is_unfolded_site_frequency_spectrum"] = args.unfolded_site_frequency_spectrum
    config_d["is_calculate_single_population_sfs"] = args.calculate_single_population_site_frequency_spectrum
    config_d["is_calculate_joint_population_sfs"] = True
    config_d["stat_label_prefix"] = args.summary_stats_label_prefix
    config_d["supplemental_labels"] = utility.parse_fieldname_and_value(args.labels)

    sscalc = sumstats.SpectrasophySummaryStatsCalculator(**config_d)
    filepath = config_d["output_prefix"] + ".obs.sumstats.tsv"
    dest = utility.open_destput_file_for_csv_writer(
            filepath=filepath,
            is_append=args.append)
    if args.append or args.no_write_header:
        is_write_header = False
    else:
        is_write_header = True
    with dest:
        writer = utility.get_csv_writer(
                dest=dest,
                delimiter=args.field_delimiter)
        try:
            results = sscalc.write_summary_stats(
                    results_csv_writer=writer,
                    results_store=None,
                    is_write_header=is_write_header)
        except Exception as e:
            sys.stderr.write("Traceback (most recent call last):\n  {}{}\n".format(
                "  ".join(traceback.format_tb(sys.exc_info()[2])),
                e))
            sys.exit(1)

if __name__ == "__main__":
    main()

