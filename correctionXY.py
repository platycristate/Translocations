import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence # generally for image processing
import cv2
from image_proc_func import mat2gray,imadjust, pairwise, Align


#_____________________________________________Constants and names________________________________
Name_dir2 = 'D:\\Lab\\Translocations_HPCA\\Cell3'
Name_dir_proc = '\\corr'
Name_seq_t = '\\p'
amin = 0
amax = 17000
MaskA = 'Fluorescence 435nm'
Mask_MasterImg = 'Fluorescence 435nm'
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


#________________________________________Lopp over .tif files________________________________________
imlist1 = [] # for stacking .tif file for 435 nm channel
imlist2 = [] # for stacking .tif file for FRET (505 nm) channel
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
    Name505 = Name505 + '.tif'
    path_to_tif = Name_dir2 + '\\' + Name

    Name_seq_temp1  = Name_dir2 + Name_dir_proc + '\\' + Name 
    Name_seq_temp2 = Name_dir2 + Name_dir_proc + '\\' + Name505 
    img_tif = Image.open(path_to_tif)

    # loops over .png imgs in a .tif file: 2 PNG imgs per loop
    frame_index = 0
    for frame1, frame2 in pairwise(ImageSequence.Iterator(img_tif)):
        frame_index +=1
#_________________preparation for alignment of images_______________________
        frame1 = mat2gray(frame1, amin, amax)  # not reliable!
        frame2 = mat2gray(frame2, amin, amax)  # not reliable!
        # prior normalization for better correction
        max1 = np.max(np.max(frame1))
        max2 = np.max(np.max(frame2)) 

        norm1 = frame1 / max1
        norm2 = frame2 / max2

        norm1 = norm1 - 0.018 # I don't know why one needs that, but it was in the script
        norm2 = norm2 - 0.018

        # enhance the contrast of frames
        norm1 = imadjust(norm1, 1)
        norm2 = imadjust(norm2, 1)
#____________________________Main part of correction___________________________________________

        # alignment of images (correction)
        if frame_index == 1:
            WarpMatrix = Align(norm1, norm2)[1] # a geometric transformation matrix
            # applies geom. transformation
            Transformed_frame = cv2.warpAffine(frame1, WarpMatrix, frame1.shape[::-1], flags=cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP)
            imlist1.append(Image.fromarray(Transformed_frame))
        else:  
           Transformed_frame = cv2.warpAffine(frame1, WarpMatrix, frame1.shape[::-1], flags=cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP)
           imlist1.append(Image.fromarray(Transformed_frame))
        imlist2.append(Image.fromarray(frame2))
#_________________________Saving the .tif files________________________________________________________

    imlist1[0].save(Name_seq_temp1, save_all=True, append_images=imlist1[1:])
    print(imlist2[0])
    imlist2[0].save(Name_seq_temp2, save_all=True, append_images=imlist2[1:])










	
	









