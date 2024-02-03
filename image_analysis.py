# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 17:05:29 2024

@author: lxiaoyang
"""
import h5py
import numpy as np
import glob
import os
import matplotlib.pylab as pylab
import pandas as pd
from skimage import io
from sklearn.cluster import KMeans
import math
#%%
#plot histogram on whole image

path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\figure\export_tif_5det_20231201'
os.chdir(path)
color = ['darkgrey','orange','limegreen']  #15deg, 15_cdeg, 45deg (ground truth)
ele = ['K','Mn','Fe','Ca']
pre = 'NNLS-'
mid = '-US_IC_'
end = 'deg.tif'
ang = ['15','15_c_','45']

for e in ele:
    fig = pylab.figure(figsize=(9, 9))
    ax = fig.subplots()
    for i,a in enumerate(ang):
        f = f'{pre}{e}{mid}{a}{end}'
        img = io.imread(f)
        img_f = img.flatten()
        hist, bins = np.histogram(img_f, bins=256)
        if a == '15_c_':
            ax.plot(bins[:-1], hist,color=color[i],linewidth=2.5,label=f'{e}_{a}')
            ax.fill_between(x = bins[:-1], y1 = hist,alpha=0.9,color=color[i])
        else:    
            ax.plot(bins[:-1], hist,color=color[i],linewidth=2.5,label=f'{e}_{a}')
            ax.fill_between(x = bins[:-1], y1 = hist,alpha=0.4,color=color[i])

        #pylab.xlim([0,2])
        #pylab.ylim([0,9000])
        ax.legend()
        pylab.title(f'{e}')
        params = {'legend.fontsize':25,
          'legend.loc':'upper right',
          'legend.frameon':False,
          'legend.facecolor':'white',
          #'figure.figsize':(11,12),
          'axes.labelsize':25,
          'axes.titlesize':25,
          'xtick.labelsize':20,
          'ytick.labelsize':20,
          'axes.linewidth': 1.5,
          'ytick.major.width': 1,
          'xtick.major.width': 1.2,
          'figure.dpi': 100,
         'xtick.major.size':9,
         'ytick.major.size':8}
        pylab.rcParams.update(params)  
        pylab.tight_layout()
#%%
#remove hotregion on Fitted K 45 deg image
path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\figure\export_tif_5det_20231201'
os.chdir(path)
f = 'Fitted-K-US_IC_45deg.tif'
fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()
img = io.imread(f)
h,w = img.shape
img_1d = img.reshape(h*w,1)
kmeans = KMeans(n_clusters=2, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(img_1d)
center = kmeans.cluster_centers_
#ax2.plot(bin_edges[0:-1], histogram,linewidth=3)
his_fit = kmeans.fit_predict(img_1d)
his_fit_2d = his_fit.reshape(h,w)
ax.imshow(his_fit_2d)
io.imsave('Fitted-K-US_IC_45deg_hot.tif',his_fit_2d)
y,x = np.where(his_fit_2d==1)
ele = ['K','Mn','Fe','Ca']
pre = 'NNLS-'
mid = '-US_IC_'
end = 'deg.tif'
ang = ['15','15_c_','45']
sp = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\figure\export_tif_5det_20231201\processed_img'
for e in ele:
    for a in ang:
        ff = f'{pre}{e}{mid}{a}{end}'
        img = io.imread(ff)
        for i,yy in enumerate(y):
            img[yy,x[i]] = 0
        io.imsave(f'{pre}{e}{mid}{a}deg_removehot.tif',img)
f = 'Fitted-K-US_IC_45deg_removehot.tif'
img = io.imread(f)
img_1d = img.reshape(h*w,1)
kmeans = KMeans(n_clusters=2, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(img_1d)
center = kmeans.cluster_centers_
#ax2.plot(bin_edges[0:-1], histogram,linewidth=3)
his_fit = kmeans.fit_predict(img_1d)
his_fit_2d = his_fit.reshape(h,w)
fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()
ax.imshow(his_fit_2d)
io.imsave('Fitted-K-US_IC_45deg_removehot_segcell.tif',his_fit_2d)
y2,x2 = np.where(his_fit_2d==1)
yb,xb = np.where(his_fit_2d == 0)
#%%
#plot average and std
path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\figure\export_tif_5det_20231201'
os.chdir(path)
ang = ['15','15_c_','45']
ele = ['K','Fe','Ca']
pre = 'NNLS-'
mid = '-US_IC_'
end = 'deg_removehot.tif'
color = ['darkgrey','orange','limegreen']
df = pd.DataFrame()
for e in ele:
    fig = pylab.figure(figsize=(9, 9))
    ax = fig.subplots()
    avg_list = []
    std_list = []
    for ii,a in enumerate(ang):
        f = f'{pre}{e}{mid}{a}{end}'
        img = io.imread(f)
        for i,ybb in enumerate(yb):
            img[ybb,xb[i]] = 0
        io.imsave(f'{pre}{e}{mid}{a}deg_removehot_maskbkg.tif',img)
        d = img.reshape(h*w,1)
        cell_pixels_value = np.delete(d,np.where(d==0))
        avg = np.average(cell_pixels_value)
        std = np.std(cell_pixels_value)
        print(f'{e}_{a}:{avg},{std}')
        avg_list.append(avg)
        std_list.append(std)
        hist, bins = np.histogram(cell_pixels_value, bins=256)
        if a == '15_c_':
            ax.plot(bins[:-1], hist,color=color[ii],linewidth=2.5,label=f'{e}_{a}')
            ax.fill_between(x = bins[:-1], y1 = hist,alpha=0.9,color=color[ii])
        else:
            ax.plot(bins[:-1], hist,color=color[ii],linewidth=2.5,label=f'{e}_{a}')
            ax.fill_between(x = bins[:-1], y1 = hist,alpha=0.4,color=color[ii])            
    df[f'{e}_avg_cell'] = avg_list
    df[f'{e}_std_cell'] = std_list
    ax.legend()
    pylab.title(f'{e}')
    #pylab.xlim([0,11])
    params = {'legend.fontsize':25,
              'legend.loc':'upper right',
              'legend.frameon':False,
              'legend.facecolor':'white',
              #'figure.figsize':(11,12),
              'axes.labelsize':25,
              'axes.titlesize':25,
              'xtick.labelsize':20,
              'ytick.labelsize':20,
              'axes.linewidth': 1.5,
              'ytick.major.width': 1,
              'xtick.major.width': 1.2,
              'figure.dpi': 100,
             'xtick.major.size':9,
             'ytick.major.size':8}
    pylab.rcParams.update(params)  
    pylab.tight_layout()
        #io.imsave(f'{pre}{e}{mid}{a}deg_removehot_maskbkg.tif',img)
#%%
#remove hot spot region use 45 deg K as reference, then apply to all other images
from sklearn.cluster import KMeans
fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()
path = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//figure//export_tiff_5det_20231129//'
f = 'NNLS-K-US_IC_45deg.tif'
img = io.imread(f)
h,w = img.shape
img_1d = img.reshape(h*w,1)
kmeans = KMeans(n_clusters=2, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(img_1d)
center = kmeans.cluster_centers_
#ax2.plot(bin_edges[0:-1], histogram,linewidth=3)
his_fit = kmeans.fit_predict(img_1d)
his_fit_2d = his_fit.reshape(h,w)
ax.imshow(his_fit_2d)

y,x = np.where(his_fit_2d==1)
'''
#apply to one image
ff = f'NNLS-Fe-US_IC_15deg_corrected.tif'
img = io.imread(ff)
for i,yy in enumerate(y):
    img[yy,x[i]] = 0
io.imsave(f'NNLS-Fe-US_IC_15deg_corrected_removehot.tif',img)
'''
#apply to all images
ele = ['K','Mn','Fe','Ca']
pre = 'NNLS-'
mid = '-US_IC_'
end = 'deg.tif'
ang = ['15','15_c','45']
sp = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//figure//export_tiff_5det_20231129//processed_img//remove_hotspot//'
for e in ele:
    for a in ang:
        ff = f'{pre}{e}{mid}{a}{end}'
        img = io.imread(ff)
        for i,yy in enumerate(y):
            img[yy,x[i]] = 0
        io.imsave(f'{sp}{pre}{e}{mid}{a}deg_removehot.tif',img)



#%%
#segmentation on remove hot spot K-45deg image, get coordinates apply to other images
from sklearn.cluster import KMeans
import cv2 as cv
from scipy import ndimage
from skimage import io
path = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//figure//export_tiff_5det_20231129//processed_img//remove_hotspot//'
os.chdir(path)
ang = ['15','15_c','45']
ele = ['K','Mn','Fe','Ca']
pre = 'NNLS-'
mid = '-US_IC_'
end = 'deg_removehot.tif'
'''
num_clu = 2
fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()
c = ['green','orange']
fr = 'NNLS-K-US_IC_45deg_removehot.tif'
img_r = io.imread(fr)
h,w = img_r.shape
img_1d = img.reshape(h*w,1)
kmeans = KMeans(n_clusters=num_clu, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(img_1d)
center = kmeans.cluster_centers_
his_fit = kmeans.fit_predict(img_1d)
his_fit_2d = his_fit.reshape(h,w)
ax.imshow(his_fit_2d)

'''
for i,a in enumerate(axo_ang):
    #fig = pylab.figure(figsize=(9, 9))
    #ax = fig.subplots()
    img_f = f'{pre}{ele}{mid}{a}{end}'
    img = io.imread(img_f)
    h,w = img.shape

    img_1d = img.reshape(h*w,1)
    kmeans = KMeans(n_clusters=2, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(img_1d)
    center = kmeans.cluster_centers_
    #ax2.plot(bin_edges[0:-1], histogram,linewidth=3)
    his_fit = kmeans.fit_predict(img_1d)
    his_fit_2d = his_fit.reshape(h,w)
    #his_fit_2d = cv.blur(his_fit_2d,(2,2))
    #his_fit_2d = ndimage.median_filter(his_fit_2d, 2)
    #ax.imshow(his_fit_2d)
   
    y,x = np.where(his_fit_2d==1)
    #cal bkg avg
    yb,xb = np.where(his_fit_2d==0)
    new_bkg = []
    for iib, yyb in enumerate(yb):
        new_bkg.append(img[yyb,xb[iib]]) #bkg region
        img[yyb,xb[iib]] = 0
        #io.imsave(f'{pre}{ele}{mid}{a}deg_removehot_bkgmasked.tif',img)
    bkg_avg = np.average(new_bkg)
    #--------------------------------
    #for cell region
    new = []
    for ii,yy in enumerate(y):
        new.append(img[yy,x[ii]]) #cell region
        #s = np.subtract(img[yy,x[ii]],bkg_avg)
        #new.append(s) #cell region-bkg_avg
    #new = np.delete(new, np.where(new==0))  #only hist cell region
    hist, bins = np.histogram(new, bins=256)
    ax.plot(bins[:-1], hist,color=c[i],linewidth=2,label=f'{a}')
    ax.fill_between(x = bins[:-1], y1 = hist,alpha=0.4,color=c[i])      
    #---------------------------
    #statistic of cell-bkg_avg
    cb_avg = np.average(new)
    cb_std = np.std(new)
    print(f'{a}deg:cell_region-bkg_avg average value{cb_avg}; std value{cb_std}')
ax.legend()
params = {'legend.fontsize':25,
          'legend.loc':'upper right',
          'legend.frameon':False,
          'legend.facecolor':'white',
          #'figure.figsize':(11,12),
          'axes.labelsize':25,
          'axes.titlesize':25,
          'xtick.labelsize':20,
          'ytick.labelsize':20,
          'axes.linewidth': 1.5,
          'ytick.major.width': 1,
          'xtick.major.width': 1.2,
          'figure.dpi': 100,
         'xtick.major.size':9,
         'ytick.major.size':8}
pylab.rcParams.update(params)  
pylab.tight_layout()

#%%
#apply K bkg mask to Fe
c = ['green','orange','royalblue']
ref_K = 'NNLS-K-US_IC_45deg_removehot_bkgmasked.tif'
Fe_15 = 'NNLS-Fe-US_IC_15deg_removehot.tif'
Fe_15_c = 'NNLS-Fe-US_IC_15deg_corrected.tif'
Fe_45 = 'NNLS-Fe-US_IC_45deg_removehot.tif'

img_K = io.imread(ref_K)
img_Fe_15 = io.imread(Fe_15)
img_Fe_15_c = io.imread(Fe_15_c)
img_Fe_45 = io.imread(Fe_45)

h,w = np.shape(img_K)
y,x = np.where(img_K==0) #K bkg coordinates
for i,yi in enumerate(y):
    img_Fe_15[yi,x[i]] = 0
    img_Fe_15_c[yi,x[i]] = 0
    img_Fe_45[yi,x[i]] = 0
    #io.imsave('NNLS-Fe-US_IC_15deg_removehot_Kbkgmask.tif',img_Fe_15)
    #io.imsave('NNLS-Fe-US_IC_15deg_corrected_Kbkgmask.tif',img_Fe_15_c)
    #io.imsave('NNLS-Fe-US_IC_45deg_removehot_K15bkgmask.tif',img_Fe_45)
img_Fe_1d = img_Fe_15.reshape(h*w,1)
img_Fe_1d = np.delete(img_Fe_1d,np.where(img_Fe_1d==0))
img_Fe_15_avg = np.average(img_Fe_1d)
std_15 = np.std(img_Fe_1d)

img_Fe_c_1d = img_Fe_15_c.reshape(h*w,1)
img_Fe_c_1d = np.delete(img_Fe_c_1d,np.where(img_Fe_c_1d==0))
img_Fe_15_c_avg = np.average(img_Fe_c_1d)
std_15_c = np.std(img_Fe_c_1d)

img_Fe_45_1d = img_Fe_45.reshape(h*w,1)
img_Fe_45_1d = np.delete(img_Fe_45_1d,np.where(img_Fe_45_1d==0))
img_Fe_45_avg = np.average(img_Fe_45_1d)
std_45 = np.std(img_Fe_45_1d)

print(f'avg 15, 15c, 45:{img_Fe_15_avg},{img_Fe_15_c_avg},{img_Fe_45_avg}; std: {std_15},{std_15_c},{std_45}')
fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()

hist15, bins15 = np.histogram(img_Fe_1d, bins=256)
hist15_c, bins15_c = np.histogram(img_Fe_c_1d, bins=256)
hist45, bins45 = np.histogram(img_Fe_45_1d, bins=256)

ax.plot(bins15[:-1], hist15,color='green',linewidth=2.5,label=f'{Fe_15}')
ax.fill_between(x = bins15[:-1], y1 = hist15,alpha=0.4,color='green')     

ax.plot(bins15_c[:-1], hist15_c,color='royalblue',linewidth=2.5,label=f'{Fe_15_c}')
ax.fill_between(x = bins15_c[:-1], y1 = hist15_c,alpha=0.4,color='royalblue')   

ax.plot(bins45[:-1], hist45,color='orange',linewidth=2.5,label=f'{Fe_45}')
ax.fill_between(x = bins45[:-1], y1 = hist45,alpha=0.4,color='royalblue')   
ax.legend()
pylab.xlim([0,0.35])
params = {'legend.fontsize':25,
          'legend.loc':'upper right',
          'legend.frameon':False,
          'legend.facecolor':'white',
          #'figure.figsize':(11,12),
          'axes.labelsize':25,
          'axes.titlesize':25,
          'xtick.labelsize':20,
          'ytick.labelsize':20,
          'axes.linewidth': 1.5,
          'ytick.major.width': 1,
          'xtick.major.width': 1.2,
          'figure.dpi': 100,
         'xtick.major.size':9,
         'ytick.major.size':8}
pylab.rcParams.update(params)  
pylab.tight_layout()
#%%
#background subtraction of projections
path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\Tomo_merk\20231228'
os.chdir(path)
fn_img = 'Ca_original_aligned_withfly42_crop.tif'
fn_ref_img = 'Ca_original_aligned_withfly42_crop_fji_autoseg.tif'
img = io.imread(fn_img)
ref_img = io.imread(fn_ref_img)
z,y,x = np.where(ref_img==0)
for i,zi in enumerate(z):
    img[zi,y[i],x[i]] = 0
io.imsave('Ca_original_aligned_withfly42_crop_bkg0.tif',np.float32(img))
#%%
path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\Tomo_merk\20231228'
os.chdir(path)
fn = 'Ca_partcle_values.csv'
#fn = 'Ca_original_recon_partcle_values.csv'
df = pd.read_csv(fn)

y1 = sorted(df['particle 3 (slice 59-74)'])
y2 = sorted(df.loc[:,'particle 2 (slice 34-43)'])
y3 = sorted(df.loc[:,'Particle 1 (slice 18-26)'])
'''
y1 = sorted(df['particle 3 (slice 63-72)'])
y2 = sorted(df.loc[:,'particle 2 (slice 36-44)'])
y3 = sorted(df.loc[:,'particle 1 (slice 19-27)'])
'''
#counts, bins = np.histogram(y1)
filtered_y1 = [item for item in y1 if not math.isnan(item)]
filtered_y2 = [item for item in y2 if not math.isnan(item)]
filtered_y3 = [item for item in y3 if not math.isnan(item)]
y1_tot = sum(filtered_y1)
y2_tot = sum(filtered_y2)
y3_tot = sum(filtered_y3)
x1 = [1]*len(filtered_y1)
x2 = [2]*len(filtered_y2)
x3 = [3]*len(filtered_y3)
cumulative_probabilities_y1 = [sum(filtered_y1[:i+1]) / y1_tot for i in range(len(filtered_y1))]
cumulative_probabilities_y2 = [sum(filtered_y2[:i+1]) / y2_tot for i in range(len(filtered_y2))]
cumulative_probabilities_y3 = [sum(filtered_y3[:i+1]) / y3_tot for i in range(len(filtered_y3))]

fig = pylab.figure(figsize=(5, 8))
ax = fig.subplots()
#pylab.hist(bins[:-1], bins=256)
#pylab.plot(x1,filtered_y1,'.',markersize=10)
#pylab.plot(x2,filtered_y2,'.',markersize=10)
#pylab.plot(x3,filtered_y3,'.',markersize=10)

ax.plot(filtered_y1,np.multiply(cumulative_probabilities_y1,100),'.',label='A',color='tomato',markersize=20)
ax.plot(filtered_y2,np.multiply(cumulative_probabilities_y2,100),'.',label='B',color='royalblue',markersize=20)
ax.plot(filtered_y3,np.multiply(cumulative_probabilities_y3,100),'.',label='C',color='lime',markersize=20) 
    

#ax.plot(filtered_y1,np.multiply(cumulative_probabilities_y1,100),'v',label='A',color='darksalmon',markersize=11)
#ax.plot(filtered_y2,np.multiply(cumulative_probabilities_y2,100),'v',label='B',color='royalblue',markersize=8)
#ax.plot(filtered_y3,np.multiply(cumulative_probabilities_y3,100),'v',label='C',color='lime',markersize=8)   

#pylab.xlim([0.045,0.3])
pylab.xlim([0.03,0.5])
pylab.ylim([0,102])
ax.set_xlabel(u'Concentration (\u03bcg/cm$^3$)')
ax.set_ylabel('Cumulative concentration probability (%)')  
ax.legend()  
params = {'legend.fontsize':22,
          'legend.loc':'lower right',
          'legend.frameon':False,
          #'figure.figsize':(11,12),
          'axes.labelsize':20,
          'axes.titlesize':29,
          'xtick.labelsize':25,
          'ytick.labelsize':25,
          'axes.linewidth': 1.8,
          'ytick.major.width': 1.2,
          'xtick.major.width': 1.2,
          'figure.dpi': 100,
         'xtick.major.size':6,
         'ytick.major.size':6}
pylab.rcParams.update(params) 
pylab.tight_layout()
#%%
#compare each particles for corrected and original
path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\Tomo_merk\20231228'
os.chdir(path)
fig = pylab.figure(figsize=(6, 8))
ax = fig.subplots()
pylab.axhline(xmin=0,xmax=0.3,y=5,linestyle='--',color = 'gray',linewidth=2)
pylab.axhline(xmin=0,xmax=0.3,y=5,linestyle='--',color = 'gray',linewidth=2)
pylab.axhline(xmin=0,xmax=0.3,y=5,linestyle='--',color = 'gray',linewidth=2)

fn_corrected = 'Ca_partcle_values.csv'
df_c = pd.read_csv(fn_corrected)
fn_original = 'Ca_original_recon_partcle_values.csv'
df_o = pd.read_csv(fn_original)

yc = sorted(df_c['particle 3 (slice 59-74)'])
yo = sorted(df_o['particle 3 (slice 63-72)'])

#yc = sorted(df_c.loc[:,'particle 2 (slice 34-43)'])
#yo = sorted(df_o.loc[:,'particle 2 (slice 36-44)'])

#yc = sorted(df_c.loc[:,'Particle 1 (slice 18-26)'])
#yo = sorted(df_o.loc[:,'particle 1 (slice 19-27)'])


color = ['tomato','darksalmon']
#color = ['royalblue','cornflowerblue']
#color = ['forestgreen','lightgreen']
marker = ['.','v']
ms = [20,11]
avg = []
std = []
quantile = []
for i,y in enumerate([yc,yo]):
    filtered_y1 = [item for item in y if not math.isnan(item)]
    y1_tot = sum(filtered_y1)
    cumulative_probabilities_y1 = [sum(filtered_y1[:i+1]) / y1_tot for i in range(len(filtered_y1))]
    ax.plot(filtered_y1,np.multiply(cumulative_probabilities_y1,100),linewidth=0,marker=marker[i],label='A',color=color[i],markersize=ms[i])
    a = sum(filtered_y1)/len(filtered_y1)
    avg.append(a)
    s = np.std(filtered_y1)
    std.append(s)
    q = np.quantile(filtered_y1, 0.8)
    quantile.append(q)
pylab.xlim([0.03,0.5])
pylab.ylim([0,102])
ax.set_xlabel(u'Concentration (\u03bcg/cm$^3$)')
ax.set_ylabel('Cumulative concentration probability (%)')  
ax.legend()  
params = {'legend.fontsize':22,
          'legend.loc':'lower right',
          'legend.frameon':False,
          #'figure.figsize':(11,12),
          'axes.labelsize':22,
          'axes.titlesize':25,
          'xtick.labelsize':21,
          'ytick.labelsize':21,
          'axes.linewidth': 1.8,
          'ytick.major.width': 1.2,
          'xtick.major.width': 1.2,
          'figure.dpi': 100,
         'xtick.major.size':6,
         'ytick.major.size':6}
pylab.rcParams.update(params) 
pylab.tight_layout()