#!/usr/bin/env python

import sys
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna

run_dir = str(sys.argv[1])
study = str(sys.argv[2])
ref_genome = str(sys.argv[3])
mut_type = str(sys.argv[4])

def main_function(run_directory, study_name, reference_genome, mutation_type):
    
    if mutation_type.startswith( "CN" ):
        matrices = scna.generateCNVMatrix(
            file_type = "ASCAT",
            input_file = run_directory+"/result/SigProfiler/"+study_name+"_"+mutation_type+"_signature_analysis/"+study_name+"_merged_"+mutation_type+".txt",
            project = mutation_type,
            output_path = run_directory+"/result/SigProfiler/"+study_name+"_"+mutation_type+"_signature_analysis/output/")
    else:
        matrices = matGen.SigProfilerMatrixGeneratorFunc(
            project=study_name+"_signature_analysis",
            genome=reference_genome,
            vcfFiles=run_directory+"/result/SigProfiler/"+study_name+"_"+mutation_type+"_signature_analysis",
            exome=False, bed_file=None, chrom_based=False, plot=True, tsb_stat=True, seqInfo=True, cushion=100)
if __name__== "__main__":
    main_function(run_dir, study, ref_genome, mut_type)