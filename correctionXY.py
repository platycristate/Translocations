import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence # generally for image processing
import cv2
from image_proc_func import mat2gray,imadjust, pairwise, Align, BGR_correction, roicolor


#_____________________________________________Constants and names________________________________
Name_dir2 = 'D:\\Lab\\Translocations_HPCA\\Cell2'
Name_dir_proc = '\\corr'
Name_seq_t = '\\p'
amin = 0
amax = 17000
MaskA = 'Fluorescence 435nm'
MaskB = 'Fluorescence FRET'
Name_MasterImg = 'Fluorescence 435nm'


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
#________________________________________Loop over .tif files________________________________________
# for stacking .tif file for 435 nm channel
# for stacking .tif file for FRET (505 nm) channel
ind = 0
for i in count_tifs: 
    imlist1 = [] # for stacking .tif file for 435 nm channel
    imlist2 = []
    ind += 1 
    Name = i
    Name435 = i[:-4] # Name for a CFP 
    Name505 = i[:-4]
    S = list(Name505)
    # change the name to "FRET", because YFP ch. is a donor of energy
    S[0:18] = 'Fluorescence  FRET'
    Name505 = ''.join(S)
    path_to_tif = Name_dir2 + '\\' + Name
    Name_seq_temp1  = Name_dir2 + Name_dir_proc + '\\' + Name435 + '.tif'
    Name_seq_temp2 = Name_dir2 + Name_dir_proc + '\\' + Name505 + '.tif'
    print(path_to_tif)
    img_tif = Image.open(path_to_tif)

    # loops over .png imgs in a .tif file: 2 PNG imgs per loop
    Img_big_float = np.zeros(np.shape(img_tif))
#_____________________________________________________________________________________________
    for frame in ImageSequence.Iterator(img_tif):  # here should be a vector summation over all frames in a .tif
        Img_big_float += np.array(mat2gray(frame, amin, amax))    
    frame_index = 0
    for frame  in ImageSequence.Iterator(img_tif):
        frame_index += 1
#_________________preparation for alignment of images_______________________
        frame = mat2gray(frame, amin, amax)
        if frame_index == 1:
            norm1 = frame/ np.max(np.max(frame))
            norm1 -=  0.018 
            norm1 = imadjust(norm1, 1)
        if frame_index == 2:
            norm2 = frame/ np.max(np.max(frame))
            norm2 -=  0.018 
            norm2 = imadjust(norm2, 1)
        if frame_index == 3:
            WarpMatrix = Align(norm1, norm2)[1]
            break
           
    frame_index = 0
    for frame in ImageSequence.Iterator(img_tif):
        frame_index += 1
        frame = mat2gray(frame, amin, amax)
        if frame_index % 2 == 0:
            frame =  cv2.warpAffine(frame, WarpMatrix, frame.shape[::-1], flags=cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP)
            imlist2.append(Image.fromarray(frame))
        else:
            imlist1.append(Image.fromarray(frame))
                            
    imlist1[0].save(Name_seq_temp1, save_all=True, append_images=imlist1[1:])
    imlist2[0].save(Name_seq_temp2, save_all=True, append_images=imlist2[1:])


