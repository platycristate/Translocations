import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence # generally for image processing
import cv2
from image_proc_func import mat2gray,imadjust, pairwise, Align, BGR_correction, roicolor


#_____________________________________________Constants and names________________________________
Name_dir2 = 'D:\\Lab\\Translocations_HPCA\\Cell31\\corr'
Name_dir_proc = '\\bgr'
Name_seq_t = '\\p'
amin = 0
amax = 17000
MaskA = 'Fluorescence 435nm'
MaskB = 'Fluorescence FRET'
Name_MasterImg = 'Fluorescence 435nm'


g = 0.2
aMask = 0.5
dXYMask = 80
dMax = 10
Nfiles = []
#_________________________________________________________________________________________________
# for loop that gets names of files in the directory
for r, d, f in os.walk(Name_dir2):
    for file in f:
        Nfiles.append(file)

#os.mkdir(Name_dir2 + Name_dir_proc) # creates a folder for processed images

# This loop counts .tif files in the directory 
count_tifs = []
for name in Nfiles:
    if name[-4:] == '.tif' and  name[0:len(MaskA)]:
        count_tifs.append(name)

#print(count_tifs)
#________________________________________Loop over .tif files________________________________________

ind = 0 # for keeping the no. of a .tif image
for i in count_tifs:
    ind += 1 
    Name = i 
    Name505 = i[:-4]
    S = list(Name505)
    # change the name to "FRET"
    S[0:18] = 'Fluorescence  FRET'
    Name505 = ''.join(S)
    Name505 = Name505 + '.tif'
    path_to_tif = Name_dir2 + '\\' + Name

    Name_seq_temp1  = Name_dir2 + Name_dir_proc + '\\' + Name 
    Name_seq_temp2 = Name_dir2 + Name_dir_proc + '\\' + Name505 
    print(Name)
    img_tif = Image.open(path_to_tif)

    # loops over .png images in a .tif file: 2 PNG images per loop
    Img_big_float = np.zeros(np.shape(img_tif)) # it reads only the dimensions of a single image
#_____________________________________________________________________________________________
    frame_ind = 0
    for frame in ImageSequence.Iterator(img_tif):
        frame_ind += 1
        #print("max value of frame:", np.array(frame).max())  # here should be a vector summation over all frames in a .tif
        Img_big_float += np.array(frame)
    #______________________________background correction___________________________________________________________

    bgr = BGR_correction(Img_big_float, d, g)
    low = (bgr + 1.5)/amax # the value can be adjusted in order to extrac the cell
    bw = roicolor(Img_big_float, low, 1)

    # Adjusting bgr coordinates for each .tif file
    if ind == 1:
        NumCellPixels = np.sum(np.sum(bw))
    NumCellPixels_new = np.sum(np.sum(bw))


    while NumCellPixels_new > NumCellPixels:
        low = low + 0.1 / amax
        bw = roicolor(Img_big_float, low, 1)
        NumCellPixels_new = np.sum(np.sum(bw)) 

    while NumCellPixels > NumCellPixels_new:
        low = low - 0.1 / amax
        bw = roicolor(Img_big_float, low, 1)
        NumCellPixels_new = np.sum(np.sum(bw))
    #_________________________coordinates of the soma_________________________
    Img_big_float = ((Img_big_float - bgr*np.ones(np.shape(Img_big_float)))) * bw;
    if ind == 1: # first .tif file
        SomaMax = np.max(np.max(Img_big_float))

        bw_soma = roicolor(Img_big_float, aMask*SomaMax/amax, 1)

        [Xval, Yval] = np.shape(Img_big_float)

        # find boundaries of soma on X-axis
        XD = np.sum(bw_soma, 1)
        X0 = 0
        while XD[X0] == 0:
            X0 += 1

        if X0 > dXYMask:
            X0 -= dXYMask
        else:
            X0 == 0

        X1 = Xval - 1
        while XD[X1] == 0:
            X1 -= 1

        if X1 < Xval-1 - dXYMask:
            X1 += dXYMask
        else:
            X1 = Xval-1

        # find boundaries of soma on Y-axis
        YD = np.sum(bw_soma, 0)
        Y0 = 0
        while YD[Y0] == 0:
            Y0 += 1
        if Y0 > dXYMask:
            Y0 -= dXYMask
        else:
            Y0 = 0

        Y1 = Yval - 1
        while YD[Y1] == 0:
            Y1 -= 1
        if Y1 < Y1 - dXYMask:
            Y1 += dXYMask
        else:
            Y1 = Yval-1


    print(NumCellPixels, NumCellPixels_new)   
    frame_index = 0
    imlist = []
    for frame in ImageSequence.Iterator(img_tif):
        frame_index += 1

        #___________subtraction of background_____________
        frame = frame - (np.ones(np.shape(frame)) * bgr) 
        frame *= bw # we have values zero everywhere but a soma with dendrites

        #____________photobleaching compensation_____________________
        if frame_index == 1 and ind == 1:
            sum0 = np.sum(np.sum(frame)) / np.sum(np.sum(bw)) * (2000/amax) / np.max(frame) # initial intensity, 2000 is the middle of the 
            # camera range, {mean intensity of the image * mean of the camera range / maximum intensity of the image}?
        sum_new = np.sum(np.sum(frame)) / np.sum(np.sum(bw))
        frame = frame * sum0/sum_new
        imlist.append(Image.fromarray(frame))
        if frame_index == 1:
            reference_img = frame
            Delta = []
        Delta.append(np.sum(np.sum(abs(frame - reference_img)))/np.sum(np.sum(reference_img))) # changes in the total intensity relatively to the first frame in the .tiff
    #____________________________ SumTranslocations____________________
        if frame_index == 1:
            img_main_a = np.zeros(np.shape(frame))
        img_main_a += frame

    imlist[0].save(Name_seq_temp1, save_all=True, append_images=imlist[1:])     
    img_main_a = img_main_a / frame_index # the mean intensity frame over all .tif stack
    
    Max, MidxA = max(Delta), Delta.index(max(Delta))
    img_tif.seek(MidxA)
    img_max_a = img_tif # now it is a single frame after the function PIL.seek()
    img_delta_a = img_max_a - img_main_a
    bw_a = roicolor(img_delta_a, dMax/amax, 1)
    if ind == 1:
        bw_aa = bw_a
    frame_index -= 1
    SumTranslA = np.zeros([frame_index, 2])
    # now we will
    img_tif_proc = Image.open(Name_seq_temp1)
    for idx in range(frame_index):
        img_tif_proc.seek(idx)
        frame = img_tif_proc
    #_____________________Sum translocations soma______________________
        Soma_bw = bw_aa[X0:X1,Y0:Y1] # coordinates of soma 
        # calculating dF for soma (F - F_mean)
        img_d_a = frame - img_main_a
        # Calculating dF/F_mean in selected region for soma
        SumTranslA[idx,1] = np.sum(np.sum(img_d_a[X0:X1, Y0:Y1] * Soma_bw)) / np.sum(np.sum(img_main_a[X0:X1, Y0:Y1] * Soma_bw))
    #_____________________Sum translocations dedndrites__________________
        soma_bw = bw_aa
        # extracting dendrite by setting soma coordinates to zero
        Dendrite_bw = np.zeros(np.shape(bw_aa[X0:X1, Y0:Y1]))
        # Calculating dF/F_mean in selected region for dendrites
        SumTranslA[idx,2] = np.sum(np.sum(img_d_a * Dendrite_bw)) / np.sum(np.sum(img_main_a * Dendrite_bw))
        






