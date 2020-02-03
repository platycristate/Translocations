import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence # generally for image processing
import cv2
from image_proc_func import mat2gray,imadjust, pairwise, Align, BGR_correction, roicolor


#_____________________________________________Constants and names________________________________
Name_dir2 = 'D:\\Lab\\Translocations_HPCA\\Cell31'
Name_dir_proc = '\\corr'
Name_seq_t = '\\p'
amin = 0
amax = 17000
MaskA = 'Fluorescence 435nm'
MaskB = 'Fluorescence FRET'
Name_MasterImg = 'Fluorescence 435nm'


g = 0.3
h = np.ones([3,3])/9 # kernel for filter
aMask = 0.5
dXYMask = 80
Nfiles = []
#_________________________________________________________________________________________________
# r=root, d=directories, f = files
# for loop gets the name of files in the directory
for r, d, f in os.walk(Name_dir2):
    for file in f:
        Nfiles.append(file)

os.mkdir(Name_dir2 + Name_dir_proc) # creates a folder for processed images

# This loop counts .tif files in the directory and find the index of Master image file
count_tifs = []
for name in Nfiles:
    if name[-4:] == '.tif' and  name[0:len(MaskA)]:
        count_tifs.append(name)

print(count_tifs)
#________________________________________Lopp over .tif files________________________________________
# for stacking .tif file for 435 nm channel
# for stacking .tif file for FRET (505 nm) channel
ind = 0
for i in count_tifs: 
    imlist1 = [] # for stacking .tif file for 435 nm channel
    imlist2 = []
    ind +=1 
    Name = i
    Name435 = i[:-4] # Name for a CFP 
    Name505 = i[:-4]
    S = list(Name505)
    # change the name to "FRET", because YFP ch. is a donor of energy
    S[0:18] = 'Fluorescence  FRET'
    Name505 = ''.join(S)
    Name505 = Name505 + '.tif'
    path_to_tif = Name_dir2 + '\\' + Name

    Name_seq_temp1  = Name_dir2 + Name_dir_proc + '\\' + Name 
    Name_seq_temp2 = Name_dir2 + Name_dir_proc + '\\' + Name505 
    print(path_to_tif)
    img_tif = Image.open(path_to_tif)

    # loops over .png imgs in a .tif file: 2 PNG imgs per loop
    Img_big_float = np.zeros(np.shape(img_tif))
#_____________________________________________________________________________________________
    for frame in ImageSequence.Iterator(img_tif):  # here should be a vector summation over all frames in a .tif
        Img_big_float += np.array(mat2gray(frame, amin, amax))

    frame_index = 0
    for frame1, frame2 in pairwise(ImageSequence.Iterator(img_tif)):
        frame_index += 1
#_________________preparation for alignment of images_______________________

        frame1 = mat2gray(frame1, amin, amax)  # not reliable!
        frame2 = mat2gray(frame2, amin, amax)  # not reliable!

        norm1 = frame1 / np.max(np.max(frame1))
        norm2 = frame2 / np.max(np.max(frame2))

        norm1 -=  0.018 # I don't know why one needs that, but it was in the script
        norm2 -=  0.018

        # enhance the contrast of frames
        norm1 = imadjust(norm1, 1)
        norm2 = imadjust(norm2, 1)
#____________________________Main part of correction___________________________________________

        # alignment of images (correction)
        if frame_index == 1:
            WarpMatrix = Align(norm1, norm2)[1] # a geometric transformation matrix
            # applies geom. transformation
            frame1 = cv2.warpAffine(frame1, WarpMatrix, frame1.shape[::-1], flags=cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP)
        else:  
           frame2 = cv2.warpAffine(frame1, WarpMatrix, frame1.shape[::-1], flags=cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP)
                # subtraction of background
        imlist1.append(Image.fromarray(frame1)); imlist2.append(Image.fromarray(frame2))




#_________________________Saving the .tif files________________________________________________________

    imlist1[0].save(Name_seq_temp1, save_all=True, append_images=imlist1[1:])
    imlist2[0].save(Name_seq_temp2, save_all=True, append_images=imlist2[1:])