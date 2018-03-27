import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.layers import Conv1D, MaxPooling1D
from keras.layers import Embedding
from keras.layers.core import Flatten

# initialization of embedding_matrix
# 64 codons, 64 words. One-hot embedding, so vector length is 64.
def setmodel(x_len):
	nb_words = 64
	embedding_vector_length = 64
	embedding_matrix = np.zeros((nb_words+1, embedding_vector_length))
	for i in range(1,nb_words+1):
		embedding_matrix[i,i-1] = 1

	model = Sequential()
	model.add(Embedding(nb_words+1,
						embedding_vector_length,
						weights=[embedding_matrix],
						input_length=x_len,
						trainable=False))
	model.add(Conv1D(40, 1, activation='relu'))
	model.add(Dropout(0.2))
	model.add(Conv1D(40, 1, activation='relu'))
	model.add(Dropout(0.2))
	model.add(Conv1D(40, 21, activation='relu'))
	model.add(Dropout(0.2))
	model.add(Flatten())
	model.add(Dense(20, activation='relu'))
	model.add(Dropout(0.2))
	model.add(Dense(1, activation='relu'))
	return model
