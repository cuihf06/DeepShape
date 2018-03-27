import numpy as np
import string
from Bio import SeqIO
import random
import pandas as pd
import DeepShape_model
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import backend as K
import os

# Initialization: TPMs and ribosome shapes.
# TPMs: 1M/#transcripts
# Ribosome shapes: 1/#codon
def initialization(path_ref):
	f_ref = open(path_ref)
	transcript_TPMs = {}
	transcript_CDS_len = {}
	transcript_shape = {}
	transcript_CDS_range = {}
	transcript_codon_IDs = {}
	for rec in SeqIO.parse(f_ref, "fasta"):
		namearray = rec.name.split('|')
		transID = namearray[0]
		# Get CDS length and CDS range. Initialize transcript shape vector.
		for a in namearray:
			if a.startswith("CDS:"):
				start, stop = map(string.atoi,a.lstrip("CDS:").split('-'))
				transcript_CDS_len[transID] = stop - start + 1
				transcript_CDS_range[transID] = [start, stop]
				transcript_codon_len = (stop - start + 1) / 3
				transcript_shape[transID] = np.ones(transcript_codon_len)
		# Get transcript codon sequences (in number).
		# The CDS range annotation is started from 1. But in python the sequence location is started from 0.
		seq = rec.seq
		transcript_codon_IDs[transID] = np.array([codon2ID(seq[(start+i*3-1):(start+i*3+3-1)]) for i in range(0, transcript_codon_len)])
	transNum = len(transcript_CDS_len)
	for key, value in transcript_CDS_len.items():
		transcript_TPMs[key] = 1000000.0/transNum
	f_ref.close()
	return [transcript_TPMs, transcript_CDS_len, transcript_CDS_range, transcript_codon_IDs, transcript_shape]

# Transfer codon to numbers: a subfunction of initialization().
def codon2ID(codon):
	tempdict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
	return tempdict[codon[0]]*16 + tempdict[codon[1]]*4 + tempdict[codon[2]]

# Get reads distributions: input alignment result file, parse it, and get the reads distributions along all transcripts.
def get_reads_distributions(path_alignment, transcript_TPMs, transcript_shape):
	f_align = open(path_alignment)
	reads_distributions = {}
	for key, value in transcript_shape.items():
		reads_distributions[key] = np.zeros(len(value))
	transmap = []
	readID = ''
	for line in f_align:
		linearray = line.strip().split('\t')
		if linearray[0] == readID:
			transmap.append([linearray[5], string.atoi(linearray[4])])
			continue
		if readID <> '':
			reads_distributions = update_reads_distributions(transmap, reads_distributions, transcript_TPMs, transcript_shape)
			transmap = [[linearray[5], string.atoi(linearray[4])]]
		readID = linearray[0]
	reads_distributions = update_reads_distributions(transmap, reads_distributions, transcript_TPMs, transcript_shape)
	f_align.close()
	return reads_distributions

# Update reads distributions: a subfunction of get_reads_distributions().
def update_reads_distributions(transmap, reads_distributions, transcript_TPMs, transcript_shape):
	total_proportion = sum([transcript_TPMs[transmap[i][0]]*transcript_shape[transmap[i][0]][int((transmap[i][1]-1)/3)] for i in range(0, len(transmap))])
	if total_proportion != 0:
		for i in range(0, len(transmap)):
			transID = transmap[i][0]
			codon_position = int((transmap[i][1]-1)/3)
			this_proportion = transcript_TPMs[transmap[i][0]]*transcript_shape[transmap[i][0]][codon_position]
			reads_distributions[transID][codon_position] += float(this_proportion)/total_proportion
	return reads_distributions

# Update TPMs: give reads distributions and return TPMs.
def update_TPMs(reads_distributions, transcript_CDS_len):
	transcript_TPMs = {}
	total_reads = 0
	for key, value in reads_distributions.items():
		transcript_TPMs[key] = sum(value)/transcript_CDS_len[key]*1000
		total_reads += transcript_TPMs[key]
	for key, value in transcript_TPMs.items():
		transcript_TPMs[key] = transcript_TPMs[key]/total_reads*1000000
	return transcript_TPMs

# Update ribosome shape model: give reads distributions, generate train data, and train the shape predict model.
def update_Ribosome_Shape_Model(transcript_codon_IDs, reads_distributions):
	transcript_shape_new = {}
	[X, Y, X_val, Y_val] = generate_train_val_data(transcript_codon_IDs, reads_distributions)
	[X_withoutZeros, Y_withoutZeros] = removeZeros(X, Y)
	[X_withoutOutliers, Y_withoutOutliers, Y_mean, Y_std] = removeOutliers(X_withoutZeros, Y_withoutZeros)
	[data_train_X, data_train_Y] = resampleData(Y_withoutOutliers, X_withoutOutliers, 10, 50000, Y_mean, Y_std)
	[data_val_X, data_val_Y] = [X_val, Y_val]
	model = TrueRibo_model.setmodel(121)
	myearlystopping = EarlyStopping(monitor='val_loss', patience=1, verbose=1)
	model.compile(loss=user_Klogmse, optimizer='nadam')
	model.fit(data_train_X, data_train_Y, validation_data=(data_val_X, data_val_Y), epochs=200, callbacks=[myearlystopping])
	for key, value in transcript_codon_IDs.items():
		tempcodonseq = map(int, np.hstack((np.zeros(60), value, np.zeros(60))).tolist())
		tempX = np.array([tempcodonseq[i:i+121] for i in range(0, len(value))])
		predictY = model.predict(tempX)[0:,0]
		transcript_shape_new[key] = shape_Normalization(predictY)
	return [transcript_shape_new, model]

def user_Klogmse(a, b):
	a_log = K.log(K.clip(a, K.epsilon(), None))
	b_log = K.log(K.clip(b, K.epsilon(), None))
	return K.mean(K.square(a_log - b_log), axis=-1)

# Generate training and validation data: a subfunction of update_Ribosome_Shape_Model().
def generate_train_val_data(transcript_codon_IDs, reads_distributions):
	selected_transcript_list = filter_transcripts(reads_distributions)
	[train_list, val_list] = split_train_val_set(selected_transcript_list)
	[X_val, Y_val] = get_X_and_Y(transcript_codon_IDs, reads_distributions, val_list)
	[X, Y] = get_X_and_Y(transcript_codon_IDs, reads_distributions, train_list)
	return [X, Y, X_val, Y_val]

def get_X_and_Y(transcript_codon_IDs, reads_distributions, transcript_list):
	X = []
	Y = []
	for transID in transcript_list:
		tempcodonseq = map(int, np.hstack((np.zeros(60), transcript_codon_IDs[transID], np.zeros(60))).tolist())
		X.extend([tempcodonseq[i:i+121] for i in range(0, len(transcript_codon_IDs[transID]))])
		Y.extend([shape_Normalization(reads_distributions[transID])[i] for i in range(0, len(reads_distributions[transID]))])
	return [np.array(X), np.array(Y)]

def removeZeros(ori_x, ori_y):
	data_order = np.array([i for i in range(0,len(ori_y)) if ori_y[i] > 0])
	new_y = ori_y[data_order]
	new_x = ori_x[data_order]
	return [new_x, new_y]

def removeOutliers(ori_x, ori_y):
	ori_y_mean = np.mean(np.log10(ori_y))
	ori_y_std = np.std(np.log10(ori_y))
	data_order = np.array([i for i in range(0,len(ori_y)) if (ori_y[i] > np.power(10, ori_y_mean-3*ori_y_std) and ori_y[i] < np.power(10, ori_y_mean+3*ori_y_std))])
	new_y = ori_y[data_order]
	new_x = ori_x[data_order]
	return [new_x, new_y, ori_y_mean, ori_y_std]

def resampleData(ori_y_withoutOutliers, ori_x_withoutOutliers, num_bin, num_data_per_bin, ori_y_mean, ori_y_std):
	binwidth = 6*ori_y_std/num_bin
	data_order = []
	for k in range(0, num_bin):
		y_tempbin_order = np.array([i for i in range(0,len(ori_y_withoutOutliers)) if (ori_y_withoutOutliers[i] > np.power(10, ori_y_mean-3*ori_y_std+k*binwidth) and ori_y_withoutOutliers[i] <= np.power(10, ori_y_mean-3*ori_y_std+(k+1)*binwidth))])
		if len(y_tempbin_order) >= num_data_per_bin:
			y_tempbin_resample_order = np.array(pd.Series(y_tempbin_order).sample(num_data_per_bin,replace=False))
		else:
			y_tempbin_resample_order = np.array(pd.Series(y_tempbin_order).sample(num_data_per_bin,replace=True))
		data_order.extend(list(y_tempbin_resample_order))
	random.shuffle(data_order)
	resampled_y = ori_y_withoutOutliers[data_order]
	resampled_x = ori_x_withoutOutliers[data_order]
	return [resampled_x, resampled_y]

# Filter transcript for model training according to these rules:
# 1. Remove genes with more than 10 continuous zeros (abnormal gene).
# 2. Remove genes with less than 5 reads per codon.
# This is a subfunction of generate_train_val_data().
def filter_transcripts(reads_distributions):
	selected_list = []
	for key, value in reads_distributions.items():
		if is_abnormal_gene(value, 10) == 0 and np.mean(value) >= 5:
			selected_list.append(key)
	return selected_list

# Determine whether a gene is abnormal (more than zero_Len continuous zeros).
def is_abnormal_gene(line_ribocounts, zero_Len):
	longest_zero_len = 0
	this_zero_len = 0
	num_codon = len(line_ribocounts)
	for i in range(0,num_codon):
		if line_ribocounts[i] == 0:
			this_zero_len = this_zero_len + 1
		else:
			this_zero_len = 0
		if this_zero_len > longest_zero_len:
			longest_zero_len = this_zero_len
	if longest_zero_len >= zero_Len:
		return 1
	else:
		return 0

# Split train and validation datasets (90% train, 10% validation)
# This is a subfunction of generate_train_val_data().
def split_train_val_set(selected_transcript_list):
	shuffled_list = selected_transcript_list
	random.shuffle(shuffled_list)
	num_train = int(len(shuffled_list)*0.9)
	return [shuffled_list[0:num_train], shuffled_list[num_train:(len(shuffled_list))]]

# Shape normalization
def shape_Normalization(distribution_before_normalized):
	return distribution_before_normalized/sum(distribution_before_normalized)*len(distribution_before_normalized)

# write result files
def write_runlog(reads_distributions, transcript_TPMs, transcript_shape, model, path_out_dir, loopnum):
	# reads_distributions
	path_reads_distributions = path_out_dir + "/reads_distributions_" + str(loopnum) + ".txt"
	f_out = open(path_reads_distributions, 'w')
	for key, value in reads_distributions.items():
	 templine = key + '\t' + '\t'.join(map(str,value.tolist())) + '\n'
	 f_out.write(templine)
	f_out.close()
	os.system('gzip ' + path_reads_distributions)
	# transcript_TPMs
	f_TPMs = open(path_out_dir + "/TPM_" + str(loopnum) + ".txt", 'w')
	for key, value in transcript_TPMs.items():
		f_TPMs.write(key + '\t' + str(value) + '\n')
	# transcript_shape
	path_transcript_shape = path_out_dir + "/transcript_shape_" + str(loopnum) + ".txt"
	f_out = open(path_transcript_shape, 'w')
	for key, value in transcript_shape.items():
	 templine = key + '\t' + '\t'.join(map(str,value.tolist())) + '\n'
	 f_out.write(templine)
	f_out.close()
	os.system('gzip ' + path_transcript_shape)
	# model weights
	path_weights = path_out_dir + "/weights_"+str(loopnum)+".hdf5"
	model.save_weights(path_weights)
