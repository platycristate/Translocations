import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence # generally for image processing
import cv2 # for goddamn normalizaion to 0,1 scale


def pairwise(iterable):
	''' function for iterating by two elements in a list'''
	a = iter(iterable)
	return zip(a, a)


#Constants and names 
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
		# converts to float
		frame1 = np.asarray(frame1)
		frame2 = np.asarray(frame2)
		cv.normalize(frame1,  0, 1, cv.NORM_MINMAX)
		print(frame1.max(), frame1.min())
		frame1 = frame1.astype('float32') 
		frame2 = frame2.astype('float32')
		# se thresholds for minimum intensity value and maximum
		frame1[frame1 <= amin] = 0; frame1[frame1 >= amax] = 1 
		frame2[frame2 <= amin] = 0; fram1[frame2 >= amax] = 1



	
	









