#!/usr/bin/env python

import os
import sys
import errno
import argparse
from os.path import exists

def parse_args(input_string, args=None):
    Description = "check samplesheet file contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument("run_bam_subsetting", help="Input string argument", default=input_string)
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

# TODO print warning if there is no purity ploidy for exome? saying can only detect allelic imbalance 
def print_warning(warning, context="Line", context_str=""):
    warning_str = "WARNING: {}".format(warning)
    if context != "" and context_str != "":
        warning_str = "WARNING: {}\n{}: '{}'".format(
            warning, context.strip(), context_str.strip()
        )
    print(warning_str)

def find_missing_patients(samples_list):
    patients = [sample[0] for sample in samples_list]
    wxs_normal_samples = [sample[0] for sample in samples_list if sample[4] == 'wxs' and sample[2] == 'normal']
    missing_patients = list(set(patients) - set(wxs_normal_samples))
    return missing_patients

def find_missing_normals(samples_list):
    normal_sample_ids = [sample[1] for sample in samples_list if sample[2] == 'normal']
    assigned_normal_ids = [sample[7] for sample in samples_list if sample[2] == 'tumour' and sample[7] != '']
    missing_normal_ids = list(set(assigned_normal_ids) - set(normal_sample_ids))
    return missing_normal_ids

def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    ["patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_id"]
    """
    sample_mapping_dict = {}
    samples_list = []
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 5
        HEADER = ["patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_id"] 
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )

            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            patient, sample_name, sample_type, bam_path, sequencing_type, purity, ploidy, normal_sample_id = lspl[: len(HEADER)] 
            patient = patient.replace(" ", "_")
            if not patient:
                print_error("patient entry has not been specified!", "Line", line)

            ## Check bam file extension
            for bam in [bam_path]:
                if bam:
                    if bam.find(" ") != -1:
                        print_error("BAM file contains spaces!", "Line", line)
                    if not bam.endswith(".bam"):
                        print_error(
                            "BAM file does not have extension '.bam'!",
                            "Line",
                            line,
                        )
            
            ## Check seqtype is only 'rnaseq' or 'wxs'
            if sequencing_type not in ["wxs","rnaseq","rna"]:
                print_error("Sequencing type input is not 'rna', 'rnaseq' or 'wxs'",
                    "Line", line
                )
            
            ## Update naming to common output 
            if sequencing_type == "rna":
                # update rna to rnaseq
                sequencing_type = "rnaseq"
            
            ## Check wxs samples have purity and ploidy inputs
            if sequencing_type == "wxs" and sample_type == "tumour" and purity == "" or sequencing_type == "wxs" and sample_type == "tumour" and ploidy == "":
                print_warning("purity or ploidy value missing for wxs tumour sample! MHC_HAMMER will only be able to detect allelic imbalance!",
                    "Line", line
                )         
            
            ## Check purity input is between 0 and 1
            if sequencing_type == "wxs":
                if purity != "" and not 0 < float(purity) <= 1:
                    print_error("purity is not between 0 and 1",
                        "Line", line
                    )
            
            ## make mapping dictionary
            sample_info = [sample_name, sample_type, bam_path, sequencing_type, purity, ploidy, normal_sample_id]  

            ## Create sample mapping dictionary 
            if patient not in sample_mapping_dict:
                sample_mapping_dict[patient] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[patient]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[patient].append(sample_info)

            ## Add sample to samples_list
            samples_list.append([patient] + sample_info)

    ## Find missing patients
    missing_patients = find_missing_patients(samples_list)
    if missing_patients:
        print("Error: The following patients are missing a WXS normal sample:")
        for patient in missing_patients:
            print(f"Patient: {patient}")
        sys.exit(1)  

    ## Find missing normal ids
    missing_normals = find_missing_normals(samples_list)
    if missing_normals:
        print("Error: The following sample id is not present in the sample_name column where sample_type == normal, but has been assigned in the normal_sample_id column:")
        for sample in missing_normals:
            print(f"sample_id: {sample}")
        sys.exit(1)     

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_id"]) + "\n")
            for patient in sorted(sample_mapping_dict.keys()):
                for idx, val in enumerate(sample_mapping_dict[patient]):
                    fout.write(",".join(["{}".format(patient)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def check_samplesheet_subset_bam_input(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    ["patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_id", "flagstat_path"]
    """
    sample_mapping_dict = {}
    samples_list = []
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 6
        HEADER = ["patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_id", "flagstat_path"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            patient, sample_name, sample_type, bam_path, sequencing_type, purity, ploidy, normal_sample_id, flagstat_path = lspl[: len(HEADER)] 
            patient = patient.replace(" ", "_")
            if not patient:
                print_error("patient entry has not been specified!", "Line", line)

            sample_name = sample_name.replace(" ", "_")
            if not sample_name:
                print_error("sample_name entry has not been specified!", "Line", line)

            ## Check bam file extension
            for bam in [bam_path]:
                if bam:
                    if bam.find(" ") != -1:
                        print_error("BAM file contains spaces!", "Line", line)
                    if not bam.endswith(".bam"):
                        print_error(
                            "BAM file does not have extension '.bam'!",
                            "Line",
                            line,
                        )

            
            ## Check flagstat file extension
            for flagstat in [flagstat_path]:
                if flagstat:
                    if flagstat.find(" ") != -1:
                        print_error("flagstat file contains spaces!", "Line", line)
                    if not flagstat.endswith(".csv"):
                        print_error(
                            "flagstat file does not have extension '.csv'!",
                            "Line",
                            line,
                        )
            
            ## Check seqtype is only 'rnaseq' or 'wxs'
            if sequencing_type not in ["wxs","rnaseq","rna"]:
                print_error("Sequencing type input is not 'rna', 'rnaseq' or 'wxs'",
                    "Line", line
                )
            
            ## Update naming to common output 
            if sequencing_type == "rna":
                # update rna to rnaseq
                sequencing_type = "rnaseq"

            ## Check wxs samples have purity and ploidy inputs
            if sequencing_type == "wxs" and purity == "" or sequencing_type == "wxs" and ploidy == "":
                print_warning("purity or ploidy value missing for wxs tumour sample! MHC_HAMMER will only be able to detect allelic imbalance!",
                    "Line", line
                )         
            
            ## Check purity input is between 0 and 1
            if sequencing_type == "wxs":
                if purity != "" and not 0 < float(purity) <= 1:
                    print_error("purity is not between 0 and 1",
                        "Line", line
                    )
            
            ## make mapping dictionary
            sample_info = [sample_name, sample_type, bam_path, sequencing_type, purity, ploidy, normal_sample_id, flagstat_path]  

            ## Create sample mapping dictionary 
            if patient not in sample_mapping_dict:
                sample_mapping_dict[patient] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[patient]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[patient].append(sample_info)

            ## Add sample to samples_list
            samples_list.append([patient] + sample_info)

    ## Find missing patients
    missing_patients = find_missing_patients(samples_list)
    if missing_patients:
        print("Error: The following patients are missing a WXS normal sample:")
        for patient in missing_patients:
            print(f"Patient: {patient}")
        sys.exit(1)  

    ## Find missing normal ids
    missing_normals = find_missing_normals(samples_list)
    if missing_normals:
        print("Error: The following sample id is not present in the sample_name column where sample_type == normal, but has been assigned in the normal_sample_id column:")
        for sample in missing_normals:
            print(f"sample_id: {sample}")
        sys.exit(1) 

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_id", "flagstat_path"]) + "\n")
            for patient in sorted(sample_mapping_dict.keys()):
                for idx, val in enumerate(sample_mapping_dict[patient]):
                    fout.write(",".join(["{}".format(patient)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    print(args)
    if args.run_bam_subsetting == "true":
        check_samplesheet(args.FILE_IN, args.FILE_OUT)
    else:
        check_samplesheet_subset_bam_input(args.FILE_IN, args.FILE_OUT)

if __name__ == "__main__":
    sys.exit(main())

