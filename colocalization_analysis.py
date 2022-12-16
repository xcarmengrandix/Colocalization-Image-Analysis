#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import random

import numpy as np
import pandas as pd
import numba

import skimage.filters
import skimage.io
import skimage.morphology

from skimage.io import imread, imshow
from skimage.measure import label, regionprops, regionprops_table
from skimage import io
from skimage.color import rgb2gray
from holoviews import dim, opts
from bokeh.themes.theme import Theme

import bebi103
import bokeh_catplot
import colorcet
import holoviews as hv
hv.extension('bokeh')

import bokeh
bokeh.io.output_notebook()
from bokeh.io import export_svgs

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from PIL import Image
from PIL import ImageDraw, ImageOps
import fnmatch
import csv
from pandas import DataFrame


# ## Import RGB image

# In[4]:


data_dir = r'C:\Users\path'

#define file name
fname = os.path.join(data_dir, 'Composite (RGB)_ERGIC3_N_dsRNA.tif')

#load the image
im = skimage.io.imread(fname)

#matadata about interpixel distance
ip = 0.216 #value depends on the objective used

#apply Gaussian blur
#sigma = 0.5
#im = skimage.filters.gaussian(
 #   im, sigma=(sigma,sigma), truncate = 3.5, multichannel = True)

bokeh.io.show(
    bebi103.image.imshow(im, cmap = 'rgb', interpixel_distance=ip, length_units='µm')
)


# In[5]:


#slice the channels we are interested in 

im_r = im[:,:,0].astype(float)
im_g = im[:,:,1].astype(float)
im_b = im[:,:,2].astype(float)


# In[6]:


#code the image to be cyan, magenta and yellow, corresponding to R, G and B channels. 
titles = ['Alpha-Tubulin', 'Nucleocapsid', 'Spike', 'merged'] #change this titles based on your proteins
ims = [im_r, im_g, im_b, im]
plots = [
    bebi103.image.imshow(
        image, frame_height=250, interpixel_distance=ip, length_units='µm', title=title
    )
    for image, title in zip(ims, titles)
]

#lock in zooming
for p in plots[1:]:
    p.x_range = plots[0].x_range
    p.y_range = plots[0].y_range
    
bokeh.io.show(bokeh.layouts.gridplot(plots,ncols=2))


# In[7]:


#import the corresponding mask of infected cells
mask =  os.path.join(data_dir,'mask.tif')
mask = skimage.io.imread(mask)
label_im = label(mask)
io.imshow(label_im)


# In[8]:


#label the regions of the input image depending on the connecitivty of the pixels to each other.    

label_im = label(mask)
regions = regionprops(label_im)
#properties = ['area','convex_area','bbox_area'] #can add more properties if needed


max_number_of_cells = np.amax(label_im)  #find the highest number of cells in all the images
    
#select only the cell of interest by setting its pixel values == 1 and all the others == 0:
cells_array = []
index = 0
for x in range(1, max_number_of_cells+1):
    cells_array.append(np.where(label_im != x, 0, label_im))
    cells_array[index] = (np.where(cells_array[index] == x, 1, cells_array[index]))
    index = index+1


# ## R channel vs G channel 

# #### Here R and G channel colocalization is quantified. You can do the same also for other combinations (e.g. B vs R). If so, you need to change the channels selected below here. 

# In[9]:


#Multiply arrays of mask image and original image (Red Channel)
single_cell_rCh = []
for index in range(0, max_number_of_cells):
    for num1, num2 in zip(cells_array[index], im_r):
        single_cell_rCh.append(num1*num2)
index = index+1

#Multiply arrays of mask image and original image (Green Channel)
single_cell_gCh = []
for index in range(0, max_number_of_cells):
    for num1, num2 in zip(cells_array[index], im_g):
        single_cell_gCh.append(num1*num2)
index = index+1

#show one cell to check that it worked
single_cell_roi_rCh = np.array_split(single_cell_rCh, np.amax(label_im))
bokeh.io.show(
    bebi103.image.imshow(
        single_cell_roi_rCh[1], interpixel_distance=ip, length_units='µm', cmap=colorcet.gray
    )
)


# In[10]:


#check it also for the second channel
single_cell_roi_gCh = np.array_split(single_cell_gCh, np.amax(label_im))
bokeh.io.show(
    bebi103.image.imshow(
        single_cell_roi_gCh[1], interpixel_distance=ip, length_units='µm', cmap=colorcet.gray
    )
)


# In[13]:


#set thresholds to remove background and select only the cells

thresh_r = 20 #change this based on your needs
thresh_g = 20 #change this based on your needs


roi_mask = []
index = 0
for index in range(0, max_number_of_cells):
    roi_mask.append((single_cell_roi_rCh[index]>thresh_r) | (single_cell_roi_gCh[index]>thresh_g))
    roi_mask[index] = skimage.morphology.remove_small_objects(roi_mask[index], min_size=3)
index=index+1


# view a cell to check threshold values
bokeh.io.show(
    bebi103.image.imshow(
        roi_mask[0], interpixel_distance=ip, length_units='µm', cmap=colorcet.gray
    )
)


# In[16]:


#make scatterplots for all cells

from bokeh.plotting import figure, output_file, save


# plot 2D histogram
def plot_hist(ch1, ch2, ch1_name='cyt', ch2_name='Nucleocapsid'):
    '''Plot a 2D histogram from two channels'''
    # get histogram data
    df_hist = pd.DataFrame({ch1_name: ch1.ravel(), ch2_name: ch2.ravel()})
    df_hist = df_hist.groupby(df_hist.columns.tolist()).size().reset_index(name='count')
    df_hist['log count'] = np.log10(df_hist['count'])
    
    # make the plot
    return hv.Points(
        data=df_hist, kdims=['cyt', 'Nucleocapsid'], vdims=['log count'],
    ).opts(
        size=7,
        cmap='magma',
        color='log count',
        colorbar=True,
        colorbar_opts={"title": "log₁₀ count"},
        frame_height=350,
        frame_width=350,
        padding=0.05,
        xlabel='α-Tubulin',
        ylabel='Nucleocapsid',
        fontsize=15,   
        #title='Colocalization G3BP1 and N', 
        #change threshold in hv.VLine(x) based on threshold used
    )*hv.VLine(20).opts(apply_ranges=True, line_width=2, color = 'black', line_dash = 'dashed') * hv.HLine(20).opts(apply_ranges=True, line_width=2, color = 'black', line_dash = 'dashed')

roi_mask_nt = []
index = 0
for index in range(0, max_number_of_cells):
    roi_mask_nt.append((single_cell_roi_rCh[index]>0) | (single_cell_roi_gCh[index]>0))
    roi_mask_nt[index] = skimage.morphology.remove_small_objects(roi_mask_nt[index], min_size=3)
index=index+1

#plot one specific cell only to check it
plot = plot_hist(im_r[roi_mask_nt[0]], im_g[roi_mask_nt[0]])

hv.renderer('bokeh').theme = 'caliber'

#LOOP PLOTS OF ALL CELLS IN THE INDEX
plots = []
index = 0
for index in range(0, max_number_of_cells):
    plots.append(plot_hist(im_r[roi_mask_nt[index]], im_g[roi_mask_nt[index]])) #no threshold set in order to show all pixels, otherwise the ones below threhsold would disappear
    plots[index] = hv.render(plots[index])
    plots[index].output_backend = "svg"
    #export_svgs(plots[index], filename = "cell_{}.svg".format(index))

plot


# In[17]:


#calculate Pearson and Manders coefficients
coloc = [] 
results = []
index = 0
for index in range(0, max_number_of_cells):
    coloc.append(bebi103.image.costes_coloc(im_r, im_g, psf_width=3, n_scramble=1000, thresh_r=0.0, roi=roi_mask[index], roi_method='all', do_manders=True))
    # print results
    results.append(print("""
    Pearson r = {0:.2f}
    prob of colocalization = {1:.2f}
    Manders coefficient for R = {2:.2f}
    Manders coefficient for G = {3:.2f}
    """.format(coloc[index].pearson_r, coloc[index].p_coloc, coloc[index].M_1, coloc[index].M_2)))


# # Export Results

# In[795]:


MandersCoeff_Gchannel = []
MandersCoeff_Rchannel = []
Pearson = []
index = 0
for index in range(0, max_number_of_cells):
    MandersCoeff_Gchannel.append(coloc[index].M_2),
    MandersCoeff_Rchannel.append(coloc[index].M_1),
    Pearson.append(coloc[index].pearson_r)

#export as csv
r = zip(MandersCoeff_Gchannel, MandersCoeff_Rchannel, Pearson) 
with open('coloc_results.csv', 'w') as s:
    w = csv.writer(s)
    w.writerow(['Manders G over R', 'Manders R over G', 'Pearson r'])
    for row in r:
        w.writerow(row)





