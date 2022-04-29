#!/usr/bin/env python

import sys
from SigProfilerExtractor import sigpro as sig

run_dir = str(sys.argv[1])
study = str(sys.argv[2])
mut_type = str(sys.argv[3])
sig_type = str(sys.argv[4])
results_dir = str(sys.argv[5])
threads = int(sys.argv[6])

def main_function(run_directory, study_name, mutation_type, signature_type, results_directory, ncpu):
    if mutation_type.startswith( "CN" ):
        curr_input_data = run_directory+"/result/SigProfiler/"+study_name+"_"+mutation_type+"_signature_analysis/output/"+mutation_type+"/"+signature_type+".matrix.tsv"
    else :
        curr_input_data = run_directory+"/result/SigProfiler/"+study_name+"_mutation_signature_analysis/output/"+mutation_type+"/"+study_name+"_signature_analysis."+signature_type+".all"
    
    sig.sigProfilerExtractor(
        input_type="matrix",
        output=results_directory,
        input_data = curr_input_data,
        opportunity_genome="GRCh38",
        context_type=signature_type,
        cpu=ncpu,
        minimum_signatures=1, maximum_signatures=10, nmf_replicates=500,
        cosmic_version=3.2)
if __name__== "__main__":
    main_function(run_dir, study, mut_type, sig_type, results_dir, threads)