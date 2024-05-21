# usage: run_trimmomatic_master.py [basedir] [input_files_list.txt] [output_directory]
# input file format:
# sample_name	fastq1	fastq2
import sys, glob, fileinput, re, math, os
from subprocess import call


# ---- parse arguments ---- #
if len(sys.argv[1:]) == 3:
	vals = sys.argv[1:]
	print( "arguments successfully loaded!" )
	print( "input file: {0}; output_directory: {1}".format( vals[0], vals[1] ))
else:
	print( "\nUsage: <run_trimmomatic_master_1.py> [input_files_list] [output_directory]" )

infile=str(vals[1])
basedir=str(vals[0])
outdir=str(vals[2])

# ---- get paired fastq file names {sample: [fq1, fq2]} ---- #
my_fastqs={}

with open(infile) as my_dat:
	for cnt, line in enumerate(my_dat):
		line = line.rstrip('\n')
		line = re.split(r'\t', line)
		my_fastqs[line[0]]=[basedir+line[1],basedir+line[2]]


# ---- submit jobs ---- #
for key in my_fastqs:
	os.system( 'sbatch run_trimmomatic_1.sh {0} {1} {2} {3}'.format(key, my_fastqs[key][0], my_fastqs[key][1], outdir) )

