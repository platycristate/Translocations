import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence # generally for image processing


#_______________________________________useful functions_____________________________________________

def mat2gray(matrix, amin, amax):
	"""MATLAB equivalent of mat2gray"""
	matrix = np.asarray(matrix)
	# converts to float
	matrix = matrix.astype('float32') 
	# set thresholds for minimum intensity value and maximum
	matrix[matrix <= amin] = 0; matrix[matrix >= amax] = 1 
	# rescale the values that are not 0 and 1 to 0/1 scale
	matrix[matrix != 1] /= matrix.max()
	return matrix

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
for i in count_tifs:
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
	for frame1, frame2 in pairwise(ImageSequence.Iterator(img_tif)):
		#_________________imitation of mat2gray from MATLAB_______________________#
		frame1 = mat2gray(frame1, amin, amax)  # not reliable!
		frame2 = mat2gray(frame2, amin, amax)  # not reliable!
		# prior normalization for better correction
		max1 = np.max(np.max(Img11))
        max2 = np.max(np.max(Img21))

        norm1 = frame1 / max1
        norm2 = frame2 / max2
       
        norm1 = norm1 - 0.018; # I don't know why one needs that, but it was in the script
        norm2 = norm2 - 0.018;
        




	
	








