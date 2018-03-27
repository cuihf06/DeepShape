import argparse

parser = argparse.ArgumentParser(usage='''python TrueRibo_RecursiveShape.py [-r PATH_REF] [-a PATH_ALIGNMENT] [-s PATH_START_TPM] [-o PATH_OUT_DIR]''')
parser.add_argument('-r', dest='path_ref', help='The destination path of reference.')
parser.add_argument('-a', dest='path_alignment', help='The destination path of STAR mapping result.')
parser.add_argument('-s', dest='path_start_TPM', help='The destination path of start TPM.')
#parser.add_argument('-f', dest='TPM_fixed', help='Set this flag to 1 if you need to fix the given TPM as unchangeable.')
parser.add_argument('-o', dest='path_out_dir', help='The directory to save results.')
parser.add_argument('-l', dest='numloops')
args = parser.parse_args()

import DeepShape_functions as DF

# Flag
# TPM_fixed = args.TPM_fixed

# Main
path_ref = args.path_ref
[transcript_TPMs, transcript_CDS_len, transcript_CDS_range, transcript_codon_IDs, transcript_shape] = DF.initialization(path_ref)

# Read in start TPM
path_start_TPM = args.path_start_TPM
file_start_TPM = open(path_start_TPM)
for line in file_start_TPM:
	transID = line.strip().split()[0]
	this_TPM = float(line.strip().split()[1])
	transcript_TPMs[transID] = this_TPM

path_out_dir = args.path_out_dir
path_alignment = args.path_alignment
numloops = int(args.numloops)+1
for loopnum in range(1, numloops):
	reads_distributions = DF.get_reads_distributions(path_alignment, transcript_TPMs, transcript_shape)
#	if TPM_fixed == 0:
#		transcript_TPMs = DF.update_TPMs(reads_distributions, transcript_CDS_len)
	transcript_TPMs = DF.update_TPMs(reads_distributions, transcript_CDS_len)
	[transcript_shape, model] = DF.update_Ribosome_Shape_Model(transcript_codon_IDs, reads_distributions)
	DF.write_runlog(reads_distributions, transcript_TPMs, transcript_shape, model, path_out_dir, loopnum)
