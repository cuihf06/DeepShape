import sys
import numpy as np

ifname_quant = sys.argv[1]
ofname_TPM = sys.argv[2]
ofname_TTPM = sys.argv[4]
ofname_TE = sys.argv[3]


trans2TE = {}
trans2TPM = {}
trans2TTPM = {}

def get_TPM(ifname_quant, trans2TPM):
	file_quant = open(ifname_quant)
	line = file_quant.readline()
	for line in file_quant:
		array = line.strip().split("\t")
		trans2TPM[array[0]] = int(float(array[3]))
	file_quant.close()
	total = sum(trans2TPM.values())
	diff = 1000000 - total
	translist = trans2TPM.keys()
	for i in range(0, diff):
		whichtransID = translist[np.random.randint(0,len(translist))]
		trans2TPM[whichtransID] = trans2TPM[whichtransID] + 1

def get_TTPM(trans2TPM, trans2TTPM, trans2TE):
	for key, value in trans2TPM.items():
		TE_before_normalization = 10**np.random.normal(0, 0.3)
		trans2TTPM[key] = value*TE_before_normalization
	total = sum(trans2TTPM.values())
	for key, value in trans2TPM.items():
		if value > 0:
			trans2TTPM[key] = trans2TTPM[key] / total * 1e6
			trans2TE[key] = trans2TTPM[key] / value
		else:
			trans2TTPM[key] = 0
			trans2TE[key] = -1

def write_files(trans2TPM, trans2TTPM, trans2TE, ifname_quant, ofname_TPM, ofname_TTPM, ofname_TE):
	file_quant = open(ifname_quant)
	file_TPM = open(ofname_TPM, 'w+')
	file_TTPM = open(ofname_TTPM, 'w+')
	file_TE = open(ofname_TE, 'w+')
	file_TPM.write("Name\tTPM\n")
	file_TTPM.write("Name\tTTPM\n")
	file_TE.write("Name\tTE\n")
	line = file_quant.readline()
	for line in file_quant:
		array = line.strip().split("\t")
		transID = array[0]
		file_TPM.write(transID + "\t" + str(trans2TPM[transID]) + "\n")
		file_TTPM.write(transID + "\t" + str(trans2TTPM[transID]) + "\n")
		file_TE.write(transID + "\t" + str(trans2TE[transID]) + "\n")
	file_quant.close()
	file_TPM.close()
	file_TE.close()
	file_TTPM.close()

# main
get_TPM(ifname_quant, trans2TPM)
get_TTPM(trans2TPM, trans2TTPM, trans2TE)
write_files(trans2TPM, trans2TTPM, trans2TE, ifname_quant, ofname_TPM, ofname_TTPM, ofname_TE)
