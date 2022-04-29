"""Add lines to VCF header.
Does not verify the format of header lines or modify the following data in any way.
This is a solution for adding the PEDIGREE information for the cancer-mode of snpEff.
Example: python add_to_vcf_header.py original.vcf '##PEDIGREE=<Derived=TUMOR,Original=NORMAL>'
"""

import sys

def process(infile, add_to_header):
    """Do the work of looping over the input file to separate header from data.
    """

    old_header = ""
    line = None

    # First, find the true header lines
    # That start with two hashes
    for line in infile:
        
        if line.startswith("##"):
            old_header += line;
        else:
            break

    old_header += "%s\n" % add_to_header 

    # Write out the new header
    sys.stdout.write(old_header)

    # Write the line just after the header
    # Due to the 'break' in the loop above
    sys.stdout.write(line)

    # Write the rest of the file
    for line in infile:
        sys.stdout.write(line)


def main():
    
    fn = sys.argv[1]
    add_header = sys.argv[2]

    process(file(fn), add_header)

if __name__ == "__main__":
    main()