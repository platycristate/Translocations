import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence, ImageEnhance # generally for image processing
import cv2


#_______________________________________useful functions_____________________________________________

def mat2gray(matrix, amin, amax):
	"""MATLAB equivalent of mat2gray"""
	matrix = np.asarray(matrix)
	# converts to float
	matrix = matrix.astype('float32') 
	# set thresholds for minimum intensity value and maximum
	matrix[matrix <= amin] = 0; matrix[matrix >= amax] = 1 
	# rescale the values that are not 0 and 1 to 0/1 scale
	matrix[matrix != 1] /= amax
	return matrix

def imadjust(matrix, perc):
	'''MATLAB equivalent of imadjust. Enhances the contrast of image by saturating
	the bottom percent of pixels and top percent pf pixels by rescaling.'''
	matrix = matrix.astype('float32')
	 # rescale an image to [0,1] range if it is in [0,255]
	if matrix.max() > 1: # if the image is not in 0/1 scale
		matrix /= matrix.max()
		print('here')
	amin = matrix.min()
	amax = matrix.max()
	# find the value amax when perc% of pixels will be greater than amax
	while (np.count_nonzero(matrix >= amax)) / (np.shape(matrix)[0]*np.shape(matrix)[1]) <= (perc/100): 
		amax = amax - 0.001
	# find the value amin when perc% of pixels will be less than amin
	while (np.count_nonzero(matrix <= amin)) / (np.shape(matrix)[0]*np.shape(matrix)[1]) <= (perc/100):
		amin = amin + 0.001
	# rescale everuthing between amin and amax to [0,1] range with amin and amax
	matrix[(matrix > amin) & (matrix < amax)] = (matrix[(matrix > amin) & (matrix < amax)] - amin) / (amax - amin)
	#matrix *=255 # this line is not required, it is just my windows 10 image viewer cannot displa 0/1 images
	contrasted_image = Image.fromarray(matrix)
	return contrasted_image 


def pairwise(iterable):
	''' function for iterating by two elements in a list'''
	a = iter(iterable)
	return zip(a, a)

#________________________________________________________________

# Constants and names 
Name_dir2 = 'D:\\Lab\\Translocations_HPCA\\Cell1'
Name_dir_proc = '\\corr'
Name_seq_t = '\\p'
amin = 0
amax = 17000
MaskA = 'Fluorescence 435nm'
Mask_MasterImg = 'Fluorescence 435nm'
Nfiles = []

# r=root, d=directories, f = files
# for loop gets the name of files in the directory
for r, d, f in os.walk(Name_dir2):
    for file in f:
        Nfiles.append(file)

#os.mkdir(Name_dir2 + Name_dir_proc) # creates a folder for processed images

# This loop counts .tif files in the directory and find the index of Master image file
count_tifs = []
for name in Nfiles:
	if name[-4:] == '.tif' and  name[0:len(MaskA)]:
		count_tifs.append(name)


# this loop works with a particular .tif file
ind = 0
for i in count_tifs: 
	ind +=1 
	Name = i
	Name435 = i[:-4] # Name for a CFP 
	Name505 = i[:-4]
	S = list(Name505)
	# change the name to "FRET", because YFP ch. is a donor of energy
	S[0:18] = 'Fluorescence  FRET'
	Name505 = ''.join(S)
	path_to_tif = Name_dir2 + '\\' + Name

	img_tif = Image.open(path_to_tif)
	# loops over .png imgs in a .tif file: 2 PNG imgs per loop
	c = 0
	for frame1, frame2 in pairwise(ImageSequence.Iterator(img_tif)):
		c +=1
		#_________________imitation of mat2gray from MATLAB_______________________#
		frame1 = mat2gray(frame1, amin, amax)  # not reliable!
		frame2 = mat2gray(frame2, amin, amax)  # not reliable!
		if ind == 1 and c == 1:
			cv2.imshow('image', frame1)
			cv2.waitKey(0)
		# prior normalization for better correction
		max1 = np.max(np.max(frame1))
		max2 = np.max(np.max(frame2)) 

		norm1 = frame1 / max1
		norm2 = frame2 / max2
		
		norm1 = norm1 - 0.018; # I don't know why one needs that, but it was in the script
		norm2 = norm2 - 0.018;

		# enhance the contrast of frames






	
	









