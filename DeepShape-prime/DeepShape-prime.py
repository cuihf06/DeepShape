import sys
import string
import numpy

# Initialize abundance of all transcript to 1/#transcripts
def transcript_initialize(path_ref):
	f_ref = open(path_ref)
	transcript_TPMs = {}
	transcript_CDS_len = {}
	for line in f_ref:
		if line.startswith('>'):
			line = line.lstrip('>').rstrip().rstrip('|')
			linearray = line.split('|')
			transID = linearray[0]
			transcript_TPMs[transID] = 1
			for a in linearray:
				if a.startswith("CDS:"):
					start, stop = a.lstrip("CDS:").split('-')
					transcript_CDS_len[transID] = string.atoi(stop) - string.atoi(start) + 1
	transNum = len(transcript_TPMs)
	for key, value in transcript_TPMs.items():
		transcript_TPMs[key] = 1000000.0/transNum
	f_ref.close()
	return [transcript_TPMs, transcript_CDS_len]

def update_transcript_readcount(transmap, transcript_readcount, transcript_TPMs):
	total_TPM = sum([transcript_TPMs[transmap[i]] for i in range(0, len(transmap))])
	for i in range(0, len(transmap)):
		transID = transmap[i]
		TPM = transcript_TPMs[transmap[i]]
		transcript_readcount[transID] += float(TPM)/total_TPM

def readcount2TPM(transcript_readcount, transcript_CDS_len):
	transcript_TPMs = {}
	total_reads = 0
	for key, value in transcript_readcount.items():
		transcript_TPMs[key] = transcript_readcount[key]/transcript_CDS_len[key]*1000
		total_reads += transcript_TPMs[key]
	for key, value in transcript_TPMs.items():
		transcript_TPMs[key] = transcript_TPMs[key]/total_reads*1000000
	return transcript_TPMs

# E step
def E_step(path_alignment, transcript_TPMs):
	f_align = open(path_alignment)
	transcript_readcount = {}
	for key, value in transcript_TPMs.items():
		transcript_readcount[key] = 0
	transmap = []
	readID = ''
	for line in f_align:
		linearray = line.strip().split('\t')
		if linearray[0] == readID:
			transmap.append(linearray[5])
			continue
		if readID <> '':
			update_transcript_readcount(transmap, transcript_readcount, transcript_TPMs)
			transmap = [linearray[5]]
		readID = linearray[0]
	update_transcript_readcount(transmap, transcript_readcount, transcript_TPMs)
	f_align.close()
	return transcript_readcount

# M step
def M_step(transcript_readcount, transcript_CDS_len):
	transcript_TPMs_new = readcount2TPM(transcript_readcount, transcript_CDS_len)
	TPM_diff = {}
	for key, value in transcript_TPMs.items():
		TPM_diff[key] = (transcript_TPMs_new[key] - value)
	return [transcript_TPMs_new, TPM_diff]

# Output TPM_diff and TPM_diff_relative of each loop
def write_runlog(f_out, loopnum, transcript_TPMs, TPM_diff):
	TPM_list = [value for key, value in transcript_TPMs.items()]
	f_out.write('Loop ' + str(loopnum) + ' TPM:\t' + '\t'.join(map(str, TPM_list)) + '\n')
	f_out.write('Loop ' + str(loopnum) + ' TPM_diff:\t' + '\t'.join(map(str, TPM_diff)) + '\n')
	TPM_diff_square = sum([numpy.power(TPM_diff[i],2) for i in range(0, len(TPM_diff))])
	TPM_diff_mean_square = TPM_diff_square/len(TPM_diff)
	f_out.write('Loop ' + str(loopnum) + ' quadratic sum & mean square:\t' + str(TPM_diff_square) + '\t' + str(TPM_diff_mean_square) + '\n')

def write_runlog_1(f_out, f_summary_out, loopnum, transcript_TPMs, TPM_diff):
	for key, value in transcript_TPMs.items():
		f_out.write(key + '\t' + str(value) + '\t' + str(TPM_diff[key]) + '\n')
	TPM_diff_square = sum([numpy.power(value,2) for key, value in TPM_diff.items()])
	TPM_diff_mean_square = TPM_diff_square/len(TPM_diff)
	f_out.write('Loop ' + str(loopnum) + ' quadratic sum & mean square:\t' + str(TPM_diff_square) + '\t' + str(TPM_diff_mean_square) + '\n')


# Main
path_alignment = sys.argv[1]
path_ref = sys.argv[2]
[transcript_TPMs, transcript_CDS_len] = transcript_initialize(path_ref)

path_out_dir = sys.argv[3]
f_summary_out = open(path_out_dir+"/summary.txt", 'w')

numloops = int(sys.argv[4])

for loopnum in range(0, numloops):
	# E step: allot reads to transcripts according to transcript TPM distribution.
	transcript_readcount = E_step(path_alignment, transcript_TPMs)
	# M step: update transcript TPM distribution according to reads number.
	[transcript_TPMs, TPM_diff] = M_step(transcript_readcount, transcript_CDS_len)
	#write_runlog(f_out, loopnum, transcript_TPMs, TPM_diff)
	f_out = open(path_out_dir + "/runlog_"+str(loopnum)+".txt", 'w')
	write_runlog_1(f_out, f_summary_out, loopnum, transcript_TPMs, TPM_diff)
	f_out.close()
	
	
	

f_out.close()
