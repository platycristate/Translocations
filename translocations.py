import matplotlib
import matplotlib.pyplot as plt
from skimage.external import tifffile
from skimage import filters
from skimage import measure
from threshold import cellMask
from skimage.viewer import ImageViewer
from skimage import data, io
import os # for work with directories 
import numpy as np

#_____________________________________________Constants and names________________________________
Name_dir = 'D:\\Lab\\Translocations_HPCA\\Cell2\\corr2'
Name_dir_proc = '\\bgr'
MaskA = 'Fluorescence 435nm'
Nfiles = []
#_________________________________________________________________________________________________
# for loop that gets names of files in the directory
for r, d, f in os.walk(Name_dir):
    for file in f:
        Nfiles.append(file)

#os.mkdir(Name_dir2 + Name_dir_proc) # creates a folder for processed images

# This loop counts .tif files in the directory 
tifs = []
for name in Nfiles:
    if name[-4:] == '.tif' and  name[0:len(MaskA)]:
        tifs.append(name)

#________________________________________Loop over .tif files________________________________________
tif_ind = 0 # for keeping the index of a .tif image
for Name in tifs:
    path_to_tif = Name_dir + '\\' + Name   
    path_for_saving = Name_dir + Name_dir_proc + '\\' + Name 
    print(Name)
    
    # reading tif image
    tiff_tensor = tifffile.imread(path_to_tif)
    
    # subtract background from every frame
    c = 0
    for c in range(len(tiff_tensor)-1):
        tiff_tensor[c] = cellMask(tiff_tensor[c], thbreshold_method="percent", percent=95)
        c += 1
    
    Img_big_float = np.sum(tiff_tensor, axis=0) / len(tiff_tensor)
    Img_big_float = cellMask(Img_big_float, thbreshold_method='percent', percent=95)
    #Img_big_float[Img_big_float > 0] = 1

    # Creating coordinates mask for soma
    Soma = cellMask(Img_big_float, thbreshold_method='percent', percent=98)
    Soma[Soma > 0] = 1
    
    # Creating coordinates mask for dendrite
    Img_big_float[Img_big_float > 0] = 1
    Dendrite = Img_big_float - Soma
    
    # photobleaching compensation    
    frame_ind = 0
    for frame in tiff_tensor:
        frame_mask = np.copy(frame)
        frame_mask[frame_mask > 0] = 1
        if frame_ind == 0:
            sum0 = frame.sum() / frame_mask.sum() * 8000 / frame.max()
        sum1 = frame.sum() / frame_mask.sum()
        frame = frame * (sum0/sum1)
        tiff_tensor[frame_ind] = frame
        if frame_ind == 0:
            base_f = frame
            Delta = []
        Delta.append(abs(frame - base_f).sum() / base_f.sum())
        frame_ind += 1

    
    # creating a mean image with baseline fluorescence
    img_base = np.sum(tiff_tensor[:3], axis=0) / 3
    
    # creating an image with maximal change in F
    img_max = tiff_tensor[Delta.index(max(Delta))]
    
    # creating an image with spots of maximal change in intensity
    img_delta = img_max - img_base
    img_delta[img_delta < 600] = 0
    img_delta[img_delta >= 600] = 1
    
    # a matrix for storing dF/F for soma and dendrites
    transl = np.zeros([len(tiff_tensor), 2])   
    
    # calculating translocations 
    frame_ind = 0
    for frame in tiff_tensor:
        dF = frame - img_base
        transl[frame_ind,0] =\
        np.sum(np.sum(dF*Dendrite*img_delta))/np.sum(np.sum(img_base*Dendrite*img_delta))
        frame_ind += 1
    plt.plot(range(len(tiff_tensor)-1), transl[:-1,0])
    plt.title(Name)
    plt.show()
        
        
        
        
                

    
    
    
    
    
    
    
    
    
    
    
    
    
    #fig, axes = plt.subplots(nrows=1,ncols=3, figsize=(8,8))
    #axes[0].imshow(Img_big_float) 
    #axes[0].set_title('Img_big_float')
    #axes[1].imshow(Soma)
    #axes[1].set_title('Soma')
    #axes[2].imshow(Dendrite)
    #axes[2].set_title('Dendrite')
    #axes[0].imshow(img_delta*Dendrite)
    #axes[0].set_title('img_delta*Dendrite')
    #axes[1].imshow(Dendrite)
    #axes[1].set_title('Dendrite')
    #axes[2].imshow(img_delta)
    #axes[2].set_title('img_delta')

    #for a in axes:
    #    a.axis('off')
    #plt.tight_layout()
    




        
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        