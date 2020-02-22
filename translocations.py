# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 16:08:46 2020

@author: Arsentii Ivasiuk
"""
import logging

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from skimage.external import tifffile
from skimage import filters
from skimage import measure
from threshold import cellMask
from skimage.viewer import ImageViewer
from skimage import data, io

from scipy.ndimage import measurements as msr

import os # for work with directories 
import numpy as np
from PIL import Image, ImageSequence # generally for image processing
from image_proc_func import mat2gray,imadjust, pairwise, Align, BGR_correction, roicolor
import pandas as pd
import cv2
import matplotlib.pyplot as plt


#_____________________________________________Constants and names________________________________
Name_dir2 = 'D:\\Lab\\Translocations_HPCA\\Cell2\\corr2'
Name_dir_proc = '\\bgr'
Name_seq_t = '\\p'
amin = 0
amax = 16000
MaskA = 'Fluorescence 435nm'
MaskB = 'Fluorescence FRET'
Name_MasterImg = 'Fluorescence 435nm'


g = 0.2
aMask = 0.5
dXYMask = 80
dMax = 20
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
plots = []; ranges = []
ind = 0 # for keeping the no. of a .tif image
for i in count_tifs:
    Name = i 
    Name505 = i[:-4]
    S = list(Name505)
    # change the name to "FRET"
    S[0:18] = 'Fluorescence  FRET'
    Name505 = ''.join(S)
    Name505 = Name505 + '.tif'
    path_to_tif = Name_dir2 + '\\' + Name

    Name_seq_temp1  = Name_dir2 + Name_dir_proc + '\\' + Name 
    print(Name)
    tiff_tensor = tifffile.imread(path_to_tif)
    print(tiff_tensor.shape)
    c = 0
    for frame in tiff_tensor:
        tiff_tensor[c] = cellMask(tiff_tensor[c], thbreshold_method="percent", percent=95)
        c += 1
    Img_big_float = np.sum(tiff_tensor, axis=0) / len(tiff_tensor)
    Img_big_float2 = np.copy(Img_big_float)
#    fig, axes = plt.subplots(nrows=1,ncols=2)
#    axes[0].imshow(Img_big_float) 

 #   for a in axes:
 #       a.axis('off')
 #   plt.tight_layout()
    
    if  ind == 0: # first .tif file 
        SomaMax = np.max(Img_big_float)
        Img_big_float[Img_big_float > SomaMax*aMask] = 1
        bw_soma = Img_big_float
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

    frame_ind = 0
    for frame in tiff_tensor:
        mask = np.copy(frame)
        mask[mask > 0] = 1
        
        #____________photobleaching compensation_____________________
        if frame_ind == 0 and ind == 0:
            sum0 = (np.sum(np.sum(frame)) / np.sum(np.sum(mask)) * 8000) / np.max(frame) # initial intensity, 2000 is the middle of the 
            # camera range, {mean intensity of the image * mean of the camera range / maximum intensity of the image}?
        sum_new = np.sum(np.sum(frame)) / np.sum(np.sum(mask))
        frame = frame * (sum0/sum_new)
        if frame_ind == 0:
            reference_img = frame
            Delta = []
        Delta.append(np.sum(np.sum(abs(frame - reference_img)))/np.sum(np.sum(reference_img)))
        frame_ind += 1
        
    img_base = np.sum(tiff_tensor[1:3], axis=0) / 2
    max_ind = Delta.index(max(Delta))
    img_max = tiff_tensor[max_ind]
    img_delta = img_max - img_base
    img_delta[img_delta < 20] = 0
    img_delta[img_delta >= 20] = 1
    if ind == 0:
        bw_aa = img_delta
        #Dendrite_bw = np.copy(bw_aa)

        
    SumTranslA = np.zeros([len(tiff_tensor), 2])   
    Soma_bw = bw_aa[X0:X1,Y0:Y1]
   # Dendrite_bw[X0:X1, Y0:Y1] = np.zeros(np.shape(Soma_bw))


    frame_ind = 0
    for frame in tiff_tensor:
        img_d_a = frame - img_base
        SumTranslA[frame_ind,0] = np.sum(np.sum(img_d_a[X0:X1, Y0:Y1] * Soma_bw)) / np.sum(np.sum(img_base[X0:X1, Y0:Y1]))
        
        #Dendrite_bw = bw_aa
        #Dendrite_bw[X0:X1, Y0:Y1] = np.zeros(bw_aa[X0:X1, Y0:Y1].shape)
        #SumTranslA[frame_ind,1] = np.sum(np.sum(img_d_a * Dendrite_bw)) / np.sum(np.sum(img_base * Dendrite_bw))
        frame_ind += 1
        
        
        
    plots.append(SumTranslA[:,0])
    ranges.append(len(tiff_tensor))
    ind += 1
plt.plot(range(ranges[0]), plots[0])
plt.plot(range(ranges[6]), plots[6])
plt.show()
        
    
        
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        