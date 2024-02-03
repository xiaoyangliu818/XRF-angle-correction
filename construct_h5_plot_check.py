# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 09:34:36 2023

@author: lxiaoyang
"""

import h5py
import numpy as np
import glob
import os
import matplotlib.pylab as pylab
import pandas as pd
#%%
inpath = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20230906//test_angular_BIO_100000_10to10_97perc'
os.chdir(inpath)
f = glob.glob('*.csv')
fig = pylab.figure(figsize=(12, 12))
new = np.zeros([3,2048])
axes = fig.subplots()
c = ['black','red','gold','green','blue']

for i,fi in enumerate(f):     
    name = os.path.splitext(os.path.basename(fi))[0]
    df = pd.read_csv(fi,names=['ind','eng','1st','2nd','3rd','4th'])
    x = df['eng']
    y = df['1st']
    '''   
    FeKL2 = np.where(x==6.39)[0][0]
    yFe = y[FeKL2]
    ynorm = y/y[FeKL2]
    '''
    inci = np.where(x==2.01)[0][0]
    yinci = y[inci]
    ynorm = y/y[inci]  
    
   # new[i,:] = ynorm  #order: 45-->60-->90
    axes.plot(x,ynorm,linestyle='-',linewidth=2,
              label=f'{name}',alpha=1)
'''
    else:
        name = os.path.splitext(os.path.basename(fi))[0]
        df = pd.read_csv(fi,names=['ind','eng','1st','2nd','3rd','4th'])
        x = df['eng']
        y = df['1st']
        FeKL2 = np.where(x==6.39)[0][0]
        yFe2 = y[FeKL2]
        scale = yFe2/yFe
        #yall = df['1st']+df['2nd']+df['3rd']+df['4th']
        ynorm = y/scale
        axes.plot(x,ynorm,linestyle='-',linewidth=2,
                  label=f'{name}',alpha=1)        

for i,fi in enumerate(f):
    name = os.path.splitext(os.path.basename(fi))[0]
    if '100_' in name:
        df = pd.read_csv(fi,names=['ind','eng','1st','2nd','3rd','4th'])
        x = df['eng']
        y = df['1st']
        FeKL2 = np.where(x==6.39)[0][0]
        #yall = df['1st']+df['2nd']+df['3rd']+df['4th']
        #ynorm = y/y[FeKL2]
        axes.plot(x,y,linestyle='-',linewidth=2,color=c[i//2],
                  label=f'{name}',alpha=1)
    else:
        df = pd.read_csv(fi,names=['ind','eng','1st','2nd','3rd','4th'])
        x = df['eng']
        y = df['1st']
        FeKL2 = np.where(x==6.39)[0][0]
        #yall = df['1st']+df['2nd']+df['3rd']+df['4th']
        #ynorm = y/y[FeKL2]
        axes.plot(x,y,linestyle='--',linewidth=1,color=c[i//2],
                  label=f'{name}',alpha=1)
'''
axes.set_yscale('log')
axes.set_xlim([0,10.5])
#axes.set_ylim(bottom=0.0001,top=100)
axes.set_ylim(bottom=0.0001,top=10)
axes.set_xlabel('Energy (keV)')
axes.set_ylabel('Counts/channel')
axes.legend()
params = {'legend.fontsize':20,
          'legend.loc':'lower right',
          'legend.frameon':False,
          'legend.facecolor':'white',
          #'figure.figsize':(11,12),
          'axes.labelsize':20,
          'axes.titlesize':20,
          'xtick.labelsize':16,
          'ytick.labelsize':16,
          'axes.linewidth': 1.5,
          'ytick.major.width': 1,
          'xtick.major.width': 1.2,
          'figure.dpi': 100,
         'xtick.major.size':9,
         'ytick.major.size':8}
pylab.rcParams.update(params)  

pylab.tight_layout()


#%%
#calculate logscale norm diff
diff45_90 = new[0] - new[2]
diff60_90 = new[1] - new[2]
diff45_60 = new[0] - new[1]
fig = pylab.figure(figsize=(12, 12))
axes = fig.subplots()
axes.plot(x,abs(diff45_90),linestyle='-',linewidth=2,label='45-90deg')
axes.plot(x,abs(diff60_90),linestyle='-',linewidth=2,label='60-90deg')
axes.plot(x,abs(diff45_60),linestyle='-',linewidth=2,label='45-60deg')
axes.set_xlim([0,10.5])
axes.set_ylim(bottom=-0.1,top=3)
axes.set_xlabel('Energy (keV)')
axes.set_ylabel('log(count/channel) angular difference')
axes.legend()
params = {'legend.fontsize':20,
          'legend.loc':'upper right',
          'legend.frameon':False,
          'legend.facecolor':'white',
          #'figure.figsize':(11,12),
          'axes.labelsize':20,
          'axes.titlesize':20,
          'xtick.labelsize':16,
          'ytick.labelsize':16,
          'axes.linewidth': 1.5,
          'ytick.major.width': 1,
          'xtick.major.width': 1.2,
          'figure.dpi': 100,
         'xtick.major.size':9,
         'ytick.major.size':8}
pylab.rcParams.update(params)  
pylab.tight_layout()

#%%
#construct h5
import h5py
import pandas as pd
import numpy as np
import os
import glob
inpath = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20230906//h5_construct//20230927_BIO//'
os.chdir(inpath)
#fn = 'bnp_fly0003.mda.h5'
pre = 'BIO_'
end = 'deg.mda.h50'
fns = [15,30,45,60,90]
time = np.ones([2,2])*33.7717825
Nv = np.ones([54,2,2])
Nv[31,:,:] = 140157.96784607437 #change US_IC to be the same as standard
Nv[42,:,:] = 440658.0596590909 #change DS_IC to be the same as standard

simudata_f = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20230906//test_angular_BIO_100000_10to10_97perc//'
sim_f = glob.glob(simudata_f+'*.csv')

for i,f in enumerate(sim_f):
    sd = pd.read_csv(f,names=['bin','eng','1','2','3','4'])
    y = sd['1']  #spectrum

    new = np.zeros([2048,2,2])
    for xi in range(2):
        for yi in range(2):       
            new[:,xi,yi] = y
#fbnp = h5py.File(fn,'r')
#mca_arr_bnp = fbnp['MAPS/mca_arr'][:]
    ysum = new[:,0,0]+new[:,0,1]+new[:,1,0]+new[:,1,1]
    fsimu = h5py.File(pre+str(fns[i])+end,'r+') #Access and modify the desired dataset
    del fsimu['MAPS/Spectra/mca_arr']
    del fsimu['MAPS/Spectra/Integrated_Spectra/Spectra']
    del fsimu['/MAPS/Spectra/Elapsed_Livetime']
    del fsimu['/MAPS/Spectra/Input_Counts']
    del fsimu['/MAPS/Spectra/Output_Counts']
    del fsimu['/MAPS/Scalers/Values']
    del fsimu['/MAPS/Spectra/Elapsed_Realtime']
    
    new_dataset = fsimu.create_dataset('MAPS/Spectra/mca_arr', data=new, shape=new.shape)  # Adjust dtype if necessary
    new_dataset_2 = fsimu.create_dataset('MAPS/Spectra/Integrated_Spectra/Spectra', data=ysum,
                                     shape=y.shape) #integrated spectra 
    new_data_time = fsimu.create_dataset('/MAPS/Spectra/Elapsed_Livetime',data=time,shape=time.shape)
    new_data_in = fsimu.create_dataset('/MAPS/Spectra/Input_Counts',data=time,shape=time.shape)
    new_data_out = fsimu.create_dataset('/MAPS/Spectra/Output_Counts',data=time,shape=time.shape)
    new_data_US = fsimu.create_dataset('/MAPS/Scalers/Values',data=Nv,shape=Nv.shape)
    new_data_time_real = fsimu.create_dataset('/MAPS/Spectra/Elapsed_Realtime',data=time,shape=time.shape)
   
   
    fsimu.flush()
    fsimu.close()

#%%
#for one h5 scan construct
import h5py
import pandas as pd
import numpy as np
import os
import glob
inpath = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//step_by_step//constructh5//'
os.chdir(inpath)
fn = 'realAXO_angle15.mda.h50'

time = np.ones([2,2])*33.7717825
Nv = np.ones([54,2,2])
Nv[31,:,:] = 140157.96784607437 #change US_IC to be the same as standard
Nv[42,:,:] = 440658.0596590909 #change DS_IC to be the same as standard

simudata_f = 'C://Research\OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//step_by_step//'
sim_f = simudata_f+'realAXO_angle15_csv.csv'


sd = pd.read_csv(sim_f,names=['bin','eng','1','2','3','4'])
y = sd['1']  #spectrum

new = np.zeros([2048,2,2]) 
for xi in range(2): 
    for yi in range(2):       
        new[:,xi,yi] = y

ysum = new[:,0,0]+new[:,0,1]+new[:,1,0]+new[:,1,1]
fsimu = h5py.File(fn,'r+') #Access and modify the desired dataset
del fsimu['MAPS/Spectra/mca_arr']
del fsimu['MAPS/Spectra/Integrated_Spectra/Spectra']
del fsimu['/MAPS/Spectra/Elapsed_Livetime']
del fsimu['/MAPS/Spectra/Input_Counts']
del fsimu['/MAPS/Spectra/Output_Counts']
del fsimu['/MAPS/Scalers/Values']
del fsimu['/MAPS/Spectra/Elapsed_Realtime']
    
new_dataset = fsimu.create_dataset('MAPS/Spectra/mca_arr', data=new, shape=new.shape)  # Adjust dtype if necessary
new_dataset_2 = fsimu.create_dataset('MAPS/Spectra/Integrated_Spectra/Spectra', data=ysum,
                                 shape=y.shape) #integrated spectra 
new_data_time = fsimu.create_dataset('/MAPS/Spectra/Elapsed_Livetime',data=time,shape=time.shape)
new_data_in = fsimu.create_dataset('/MAPS/Spectra/Input_Counts',data=time,shape=time.shape)
new_data_out = fsimu.create_dataset('/MAPS/Spectra/Output_Counts',data=time,shape=time.shape)
new_data_US = fsimu.create_dataset('/MAPS/Scalers/Values',data=Nv,shape=Nv.shape)
new_data_time_real = fsimu.create_dataset('/MAPS/Spectra/Elapsed_Realtime',data=time,shape=time.shape)
    
fsimu.flush()
fsimu.close()

#%%
#check
import matplotlib.pylab as pylab

fns = 'AXO_3nm_15deg_times10.mda.h50'
fsimu = h5py.File(fns,'r')
a = fsimu['MAPS/Spectra/mca_arr'][:]
b = fsimu['MAPS/Spectra/Integrated_Spectra/Spectra'][:]
fsimu.close()

t = pd.read_csv(sim_f[0],names=['bin','eng','1','2','3','4'])
y = t['1']

pylab.plot(a[:,0,0],label='[0,0]')
pylab.plot(a[:,0,1],label='[0,1]')
pylab.plot(a[:,1,0],label='[1,0]')
pylab.plot(a[:,1,1],label='[1,1]')

pylab.plot(b,label='integrated')

#pylab.plot(y)
#compare with original integrated spectrum
#fig = pylab.figure(figsize=(12, 12))
#ax = fig.subplots()
#p = pd.read_csv('//micdata//data1//bnp//2023-1//simu_test//output//AXO_3nm_15deg.mda.h5_NNLS.csv')

#x = p['Energy']
#ys = p['Spectrum']
#pylab.plot(x,ys,'.',label='AXO_3nm_15deg.mda.h5_NNLS')
#pylab.plot(x,b,label='integrated')

#pylab.yscale('log')
#pylab.xlim([0,11]) #overview
#pylab.ylim([0.1,10000])
#pylab.legend()
#pylab.yscale('log')


pylab.xlabel('bin')
pylab.ylabel('counts/channel')
pylab.xlim([0,2048]) #overview
pylab.ylim([0.01,10000])
pylab.legend()

#%%
#to check integrated spec
inpath = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20230906//h5_construct//20230927_BIO//'
os.chdir(inpath)
#fn = 'bnp_fly0003.mda.h5'
pre = 'BIO_'
end = 'deg.mda.h50'
fns = [15,30,45,60,90]
for i,f in enumerate(sim_f):
    fsimu = h5py.File(pre+str(fns[i])+end,'r') #Access and modify the desired dataset
    b = fsimu['MAPS/Spectra/Integrated_Spectra/Spectra'][:]
    pylab.plot(b,label=f)
    fsimu.close()
pylab.yscale('log')
pylab.xlabel('bin')
pylab.ylabel('counts/channel')
pylab.xlim([0,1100]) #overview
pylab.ylim([0.1,11000])
#%%
import os
f = 'C://Research//OneDrive - Argonne National Laboratory//anl//Yu-Chung//wegreen//Supporting_Materials_final//'
os.chdir(f)

for i in np.arange(1,20):
    os.mkdir(f'Exhibit {i}')

#%%
import glob
import os
import pandas as pd
f = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//quantification//20230927_bio_noice//'
os.chdir(f)
files = glob.glob('*.csv')
#fig = pylab.figure(figsize=(12, 12))
#ax = fig.subplots()
#fns = 'AXO_3nm_15deg.mda.h5'
#fsimu = h5py.File(fns,'r')
#b = fsimu['MAPS/Spectra/Integrated_Spectra/Spectra'][:]
#ax.plot(x,b,'-',label='15deg_integrated_spec',linewidth=2)
#c = ['b','r','green','gray','pink']
for i,fi in enumerate(files):
    a = pd.read_csv(fi)
    n = os.path.splitext(os.path.basename(fi))[0]
    fig = pylab.figure(figsize=(10, 10))
    ax = fig.subplots()
    x = a['Energy']
    ys = a['Spectrum']
    yf = a['Fitted']
    ax.plot(x,ys,'.',label=f'{n}')
    ax.plot(x,yf,'-',label=f'{n}_fit',color='r')
    pylab.yscale('log')
    pylab.xlim([0,11]) #overview
    pylab.ylim([0.1,10000])
    pylab.legend()
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('log(count/channel)')
params = {'legend.fontsize':22,
          'legend.loc':'upper right',
          'legend.frameon':False,
          'legend.facecolor':'white',
          #'figure.figsize':(11,12),
          'axes.labelsize':22,
          'axes.titlesize':22,
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
# to check
f = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//quantification//after fit//'
os.chdir(f)
fns = 'AXO_3nm_15deg.mda.h5'
fsimu = h5py.File(fns,'r')
c = fsimu['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec'][:]
a = fsimu['MAPS/Scalers/us_amp_num'][:]






xlist = []
ylist = []
speclist = []
for x in range(a.shape[2]+1):
    for y in range(a.shape[1]+1):
        spec = list(a[:,y,x])
        xlist.append(x)
        ylist.append(y)
        speclist.append(spec)
        print(f'done{x,y}')
df = pd.DataFrame(zip(xlist,ylist,speclist),columns=['x','y','xrf spectrum'])




#%%
#compare simulated and experiment axo
import os
import matplotlib.pylab as pylab
import h5py
import numpy as np

path = '//micdata//data1//bnp//2023-1//AXO_varyangle//img.dat//'
os.chdir(path)
exp = 'bnp_fly0076.mda.h5' #15deg
simu = '//micdata//data1//bnp//2023-1//simu_realAXO//img.dat//realAXO_angle15.mda.h50'

exp_h5 = h5py.File(exp,'r')
simu_h5 = h5py.File(simu,'r')

exp_eng = exp_h5['MAPS/Spectra/Energy'][:]
exp_mca = exp_h5['MAPS/Spectra/mca_arr'][:] #(2048,11,11)
exp_int = exp_h5['MAPS/Spectra/Integrated_Spectra/Spectra'][:]
y= exp_mca.shape[1]
x= exp_mca.shape[2]

exp_int_avg = np.divide(exp_int,x*y)

simu_int_avg = simu_h5['MAPS/Spectra/Integrated_Spectra/Spectra'][:]/4/5
simu_eng = simu_h5['MAPS/Spectra/Energy'][:]

fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()

ax.plot(simu_eng,simu_int_avg,'-',color='red',linewidth=2.5,label='simu_15')
ax.plot(exp_eng,exp_int_avg,'-',color='forestgreen',linewidth=2.5,label='exp_15')
pylab.yscale('log')
pylab.xlim([1,11]) #overview
pylab.ylim([0.01,200])
pylab.legend()
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (counts/s)')
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
#plot exp AXO fitted Fe and Cu counts
ang = [0,15,45]
ele = ['Fe','Cu']
Fe_avg_7det = [47.71,69.25,96.45,2.85,3.69,3.6] #[0,15,45,std_0deg,std_15deg,std_45deg]
Cu_avg_7det = [49.76,74.44,98.93,2.95,4.05,4.23] #[0,15,45,std_0deg,std_15deg,std_45deg]
Fe_avg_5det = [61.75,67.08,93.86,3.77,4.06,4.18]
Cu_avg_5det = [62.93,70.6,94.77,3.86,4.09,4.69]
Fe_det0 = [78.18,81.92,113.15,7.94,8.58,8.75]
Cu_det0 = [79.59,83.33,112.3,9.21,8.57,9.16]
Fe_simu = [67.02,72.38,101.21]
Cu_simu = [49.15,51.91,72.37]

y1 = Fe_det0#Fe_avg_5det
y2 = Fe_simu
y3 = Cu_det0#Cu_avg_5det
y4 = Cu_simu

fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()
ax.errorbar(ang, y1[0:3], yerr=y1[-3:], fmt='o--',linewidth=2.5, markersize=15,color='royalblue',label='Fe_5det',capsize=5)
#ax.errorbar(ang, y2[0:3], fmt='o-', linewidth=2.5, markersize=15,color='royalblue',label='Fe_simu')
ax.errorbar(ang, y3[0:3], yerr=y3[-3:], fmt='o--',linewidth=2.5, markersize=15,color='orange',label='Cu_5det',capsize=5)
#ax.errorbar(ang, y4[0:3], fmt='o-', linewidth=2.5, markersize=15,color='orange',label='Cu_simu')

pylab.legend()
pylab.xlim([-5,50]) #overview
#pylab.ylim([40,120])
pylab.ylim([40,130])
pylab.legend()
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Intensity (counts/s)')
params = {'legend.fontsize':25,
          'legend.loc':'upper left',
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
#check fit results
path = r'\\micdata\data1\bnp\2022-3\Merk_xyl_fly86\output'
os.chdir(path)
fig = pylab.figure(figsize=(9, 9))
ax = fig.subplots()
f = 'bnp_fly0086.mda.h51_Fitted.csv'
df = pd.read_csv(f)
x = df['Energy']
yr = df['Spectrum']
yf = df['Fitted']
ax.plot(x,yr,'-',label='Raw')
ax.plot(x,yf,'-',label='Fitted')
pylab.legend()
pylab.xlim([1,11])
pylab.yscale('log')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (counts/s)')
params = {'legend.fontsize':25,
          'legend.loc':'lower left',
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
#plot histogram on whole image
from skimage import io

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
from sklearn.cluster import KMeans
#correction Fe
sv = 1.366 #the value respect to 15 deg
f = 'NNLS-Fe-US_IC_15deg.tif'
img = io.imread(f)
img_c = np.divide(img,sv)  #correct 45 deg sample image quantified by 15 deg AXO
#io.imsave('NNLS-Fe-US_IC_15deg_corrected.tif',img_c)
h,w = img_c.shape
img_1d = img_c.reshape(h*w,1)
kmeans = KMeans(n_clusters=2, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(img_1d)
center = kmeans.cluster_centers_
#ax2.plot(bin_edges[0:-1], histogram,linewidth=3)
his_fit = kmeans.fit_predict(img_1d)
his_fit_2d = his_fit.reshape(h,w)

ax.imshow(his_fit_2d)

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
#percentage error (PErr): according to Mingyuan's self-absorption paper
#PErr = |Img_correction - Img_groundtruth|/Img_groundtruth
#generate PErr image
#Need to do:
    #if we consider 45 deg AXO quantified image as ground truth, 
    #1. PErr between 15 deg and 45 deg for K, Ca and Fe
    #2. After correction, PErr between 15 deg and 45 deg for K, Ca and Fe
import matplotlib as mpl    
path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\figure\export_tif_5det_20231201'
os.chdir(path)
ele = ['K']
ang = ['15','15_c_']
ra = '45'
pre = 'NNLS-'
mid = '-US_IC_'
end = 'deg_removehot_maskbkg.tif'

for e in ele:
    rf = f'{pre}{e}{mid}{ra}{end}'
    imgf = io.imread(rf)
    y,x = imgf.shape
    for a in ang:
        f = f'{pre}{e}{mid}{a}{end}'#ref img: 15deg
        img = io.imread(f)
        Pu = np.abs(np.subtract(img,imgf)) 
        PErr = np.divide(Pu,imgf)
        PErr[np.isnan(PErr)] = 0
        PErr = np.multiply(PErr,100)
        #io.imshow(PErr)
        PErr_1d = PErr.reshape(y*x)
        PErr_1d = np.delete(PErr_1d,np.where(PErr_1d==0))
        avg = np.average(PErr_1d)
        print(avg)
        fig = pylab.figure(figsize=(9, 9))
        ax = fig.subplots()
        ax.imshow(PErr,vmin=0,vmax=25,cmap='viridis')
        #ax.set_title(f'{e}-{a}-45')
        norm = mpl.colors.Normalize(vmin=0.0, vmax=3.0)
        cmap2 = mpl.cm.ScalarMappable(norm=norm, cmap='viridis')
        cmap2.set_array([])
        cbar = fig.colorbar(cmap2,orientation='vertical')
   #io.imsave(f'PErr_{e}_45-15deg.tif',PErr)


#%%
#change ug/cm2 to weight percent between K, Ca, Fe
path = 'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//figure//5det_fly90_compare_tif_img//'
os.chdir(path)
axo_ang = [15,45]
ele = ['K','Ca','Fe']
pre = 'NNLS-'
mid = '-US_IC_'
end = 'deg.tif'


for a in axo_ang:
    f1 = f'{pre}{ele[0]}{mid}{a}{end}'
    imgK = io.imread(f1)
    f2 = f'{pre}{ele[1]}{mid}{a}{end}'
    imgCa = io.imread(f2)
    f3 = f'{pre}{ele[2]}{mid}{a}{end}'
    imgFe = io.imread(f3)
    img_s = np.add(imgK,imgCa,imgFe)
    w_K = np.divide(imgK,img_s)
    w_Ca = np.divide(imgCa,img_s)
    w_Fe = np.divide(imgFe,img_s)
    io.imsave(f'{pre}{ele[0]}{mid}{a}deg_wtp.tif',w_K)
    io.imsave(f'{pre}{ele[1]}{mid}{a}deg_wtp.tif',w_Ca)
    io.imsave(f'{pre}{ele[2]}{mid}{a}deg_wtp.tif',w_Fe)
#%%
#modify AXO flyscan Fe and Cu counts/s
#change flyscan MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec Fe and Cu value for each pixel (Fe:10, Cu:12)
import h5py

path = '//micdata//data1//bnp//2023-1//AXO_varyangle//img.dat//'
os.chdir(path)
scanid = [76]
det = [0,1,2,3,4]
pre = 'bnp_fly'
end = '.mda.h5'
Fe_45 = 1.366
Cu_45 = 1.366

for d in det:
    f = f'{pre}{str(scanid[0]).zfill(4)}{end}{d}'
    fi = h5py.File(f,'r+')
    cps = fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec']
    Fe = fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec'][10,:,:]#original Fe intensity in each pixel
    Fe_c = np.multiply(Fe,Fe_45)
    Cu = fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec'][12,:,:]#original Cu intensity in each pixel
    Cu_c = np.multiply(Cu,Cu_45)
    cps[10,:,:] = Fe_c
    cps[12,:,:] = Cu_c
    del fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec']
    new_cps = fi.create_dataset('MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec',data=cps,shape=cps.shape)
   
   
    fi.flush()
    fi.close()
'''
#for h5 scan
f = f'{pre}{str(scanid[0]).zfill(4)}{end}'
fi = h5py.File(f,'r+')
cps = fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec']
Fe = fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec'][10,:,:]#original Fe intensity in each pixel
Fe_c = np.multiply(Fe,Fe_45)
Cu = fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec'][12,:,:]#original Cu intensity in each pixel
Cu_c = np.multiply(Cu,Cu_45)
cps[10,:,:] = Fe_c
cps[12,:,:] = Cu_c
del fi['MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec']
new_cps = fi.create_dataset('MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec',data=cps,shape=cps.shape)
   
   
fi.flush()
fi.close()
'''
#%%
#Check Merk 2022-3 images
import os
import matplotlib.pyplot as plt 
path = r'\\micdata\data1\bnp\2022-3\Compare_Merk\img'
os.chdir(path)
files = sorted(glob.glob('*.tif'))
fig = pylab.figure(figsize=(25, 12))
ax = fig.subplots(2,5)
for i in np.arange(0,len(files),2):
    img_c = io.imread(files[i]) #corrected image
    img_o = io.imread(files[i+1])
    s1 = ax[0,int(i/2)].imshow(img_c,vmin=0.0, vmax=1.2)
    ax[0,int(i/2)].title.set_text(f'{files[i][-15:]}')
    s2 = ax[1,int(i/2)].imshow(img_o,vmin=0.0, vmax=1.2)
    ax[1,int(i/2)].title.set_text(f'{files[i+1][-15:]}')
#fig.colorbar('all')
    
#%%
#create folder
#s = np.hstack((np.arange(52,53,1),np.arange(60,77,2),np.arange(80,99,2),np.arange(104,107,2)))
s = np.hstack((np.arange(114,115,2),np.arange(117,128,2),np.arange(131,140,4),np.arange(143,146,1),np.arange(154,157,2),np.arange(163,165,1),np.arange(180,187))) #negative angle

path = r'\\micdata\data1\bnp\xyl_tomo'
for si in s:
    directory = f'fly{si}'
    path2 = os.path.join(path, directory) 
    os.mkdir(path2)
               
img = 'img.dat'
fx = 'flyXRF'
mda = 'mda'
for si in s:
    d = f'fly{si}'
    p = os.path.join(path,d,fx)
    os.mkdir(p)

#%%
#copy files to folder
import glob
import shutil
path = r'\\micdata\data1\bnp\2022-3\Merk\mda'
os.chdir(path)
s = np.hstack((np.arange(52,53,1),np.arange(60,77,2),np.arange(80,99,2),np.arange(104,107,2))) #positive angle
#s = np.hstack((np.arange(114,115,2),np.arange(117,128,2),np.arange(131,140,4),np.arange(143,146,1),np.arange(154,157,2),np.arange(163,165,1),np.arange(180,187)))#negative angle

'''
for si in s:
    files = glob.glob(f'bnp_fly_{si}_???.nc')
    sf = '//micdata\\data1\\bnp\\xyl_tomo\\'
    sp = f'{sf}fly{si}\\flyXRF\\'
    for f in files:      
        shutil.copy2(f, sp+f)
    files_ref = glob.glob('bnp_fly_3_???.nc')
    for fr in files_ref:
        shutil.copy2(fr,sp+fr)
    print(f'Done{si}')


#for mda file
for si in s:
    sif = str(si).zfill(4)
    f = f'bnp_fly{sif}.mda'
    sf = '//micdata\\data1\\bnp\\xyl_tomo\\'
    sp = f'{sf}fly{si}\\mda\\'    
    shutil.copy2(f, sp+f)
    files_ref = f'bnp_fly0003.mda'
    shutil.copy2(files_ref,sp+files_ref)
    print(f'Done{si}')
'''
#for override file
name = 'maps_fit_parameters_override'
end = ['.txt','.txt0','.txt1','.txt2','.txt3','.txt4']
o = []
stdinfo = 'maps_standardinfo.txt'
for e in end:
    of = f'{name}{e}'
    o.append(of)
for si in s:
    #for oo in o:
        #sf = '//micdata\\data1\\bnp\\xyl_tomo\\'
        #sp = f'{sf}fly{si}\\' 
        #shutil.copy2(oo,sp+oo) 
    sf = '//micdata\\data1\\bnp\\xyl_tomo\\'
    sp = f'{sf}fly{si}\\' 
    shutil.copy2(stdinfo,sp+stdinfo)
    print(f'Done{si}')

    
#%%
#check angle for scans
import h5py
import math
path = r'\\micdata\data1\bnp\2022-3\Merk\img.dat.5det.ele'
os.chdir(path)
s = np.hstack((np.arange(52,53,1),np.arange(60,77,2),np.arange(80,99,2),np.arange(104,107,2))) #positive angle
sn = np.hstack((np.arange(114,115,2),np.arange(117,128,2),np.arange(131,140,4),np.arange(143,146,1),np.arange(154,157,2),np.arange(163,165,1),np.arange(180,187)))
#files = glob.glob('*.h5')   
ang_list = []
for si in sn:
    sif = str(si).zfill(4)
    f = f'bnp_fly{sif}.mda.h5'
    ff = h5py.File(f,'r')
    ang = ff['MAPS/Scan/Extra_PVs/Values'][:][3].decode("utf-8")
    ang = round(float(ang))
    ang_list.append(ang)

#%%
Cu_ratiolist = [0.8721044456207909,0.9267666901161332,0.9556148703918238,0.9730661410968833,0.9889243605438298,1.0249504255858402,
 1.0469414212803092,1.0723638328358656,1.101753603435355,1.1356964775859066,1.2201088181383712,1.2724140157793917,
 1.3330614992919578,1.4036716863987435,1.4863425563239157,1.58383834074124,1.6998742585425035,1.8395603056194345,
 2.0101188641094248,2.222096144648276,1.174874206076842,2.4915155132367035]
Fe_ratiolist = [0.8262840067342844,0.9053120883555604,0.9462745902751929,0.9688846842660424,
 0.9876148019360047,1.0268062879149917,1.0498027701130097,1.07604226899489,1.1061319352609407,
 1.1407031526270923,1.2262763701298578,1.2791575663099903,1.3404023558790366,1.4116487590363118,
 1.495013822286133,1.5932840407574311,1.7102023005956704,1.8509148257109997,2.0226936863290055,
 2.236156358962175,1.1804693387474237,2.5074333663234314]
#ang_pos = [2,4,7,10,13,19,22,25,28,31,37,40,43,46,49,52,55,58,61,64,34,67]
sn = np.hstack((np.arange(114,115,2),np.arange(117,128,2),np.arange(131,140,4),np.arange(143,146,1),np.arange(154,157,2),np.arange(163,165,1),np.arange(180,187)))

ang_pos_neg = [-31,-37,-40,-43,-46,-49,-52,-58,-64,-70,-34,-55,-67,0,-3,-6,-9,-61,-76,-79,-82,-11,-14,-73]
s = np.hstack((np.arange(52,53,1),np.arange(60,77,2),np.arange(80,99,2),np.arange(104,107,2))) #positive angle
Fe_ratio_neg_list = [1.1701619693788374,1.253328772368374,1.3055593917255683,1.3664622925703247,1.4376561229999556,
 1.5212549675123896,1.6200583560846444,1.8797994278213748,2.269068193739679,2.902552849462617,1.2085195466055567,
 1.737840028239216,2.5434612225691144,0.6961067208807252,1.1277447523501307,1.0614851981674336,
 1.0463642225988348,2.0532928025257133,4.089855182076953,5.1702719706060725,7.052813945347865,1.0447700063695473,
 1.049197490500609,3.390691365512972]

Cu_ratio_neg_list=[1.1546017908995447,1.2374699578869506,1.2893576914552896,1.3497858025349476,1.420362280818111,
 1.503183207424446,1.6010211812067108,1.8580974839735986,2.2432178808838232,2.8698051829793454,1.19287562880158,
 1.717611222925764,2.5146370098118567,0.7772437330726222,1.0677589388001176,1.0286556519229388,1.0217020939390056,
 2.0297563690792115,4.044017998410129,5.112456264242147,6.974064344939213,1.022927728169784,1.0299455732151424,
 3.352577216442493]
Fe = 0.51
Cu = 0.23

Cu_c = np.divide(Cu,Cu_ratiolist)
Fe_c = np.divide(Fe,Fe_ratiolist)
Fe_c_neg = np.divide(Fe,Fe_ratio_neg_list)
Cu_c_neg = np.divide(Cu,Cu_ratio_neg_list)
sf = '//micdata\\data1\\bnp\\xyl_tomo\\'
for i,si in enumerate(sn):
    sp = f'{sf}fly{si}\\' 
    std = sp+'maps_standardinfo.txt'
    rf = open(std,'r')
    data = rf.read()
    dc = data.replace(f'{Cu}',f'{Cu_c_neg[i]}')
    rf.close()
    rf = open(std,'w')
    rf.write(dc)
    rf.close()
    print(f'Done{si}')

#%% show image projections of original (quantified by 15 deg AXO) and corrected tomo

path = r'C:\Research\OneDrive - Argonne National Laboratory\anl\xrf\sim_test\20231020\batch\Tomo_merk\20231228'
os.chdir(path)
f = ['Ca_aligned_crop.tif']
for fi in f:
   # fig= pylab.figure(figsize=(12, 12))
    #ax = fig.subplots()
    #ax = fig.subplots(6,7)

    #ax = axes.ravel()
    img = io.imread(fi)
    z,y,x = img.shape
    for zi in range(z):
        fig= pylab.figure(figsize=(6, 9))
        ax = fig.subplots()
        ax.imshow(img[zi],vmin=0,vmax=3)
        #row_idx = zi // 7
        #col_idx = zi % 7
        #ax[row_idx,col_idx].imshow(img[zi],vmin=0,vmax=1.5)
        #ax[row_idx,col_idx].set_title(f'{zi}')
        ax.set_title(f'{zi}')
        pylab.tight_layout()
        pylab.savefig(f'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//Tomo_merk//20231228//Ca_projections_aligned_saved//{zi}.tif')
        pylab.close()
        #ax.imsave(f'C://Research//OneDrive - Argonne National Laboratory//anl//xrf//sim_test//20231020//batch//Tomo_merk//20231228//Ca_projections_aligned_saved//{zi}.tif',img[zi])
    #pylab.show()   

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
#plot average and std of Fe,Ca and K for 15 and 45 deg
fig= pylab.figure(figsize=(7, 7))
K_avg = [1.229437,0.9067044,0.8982655]#15,45,15_c
K_std = [0.57112813,0.4212045,0.41728422]
Ca_avg = [0.27803403,0.20504889,0.20314042]
Ca_std = [0.23449878,0.17294185,0.17133221]
Fe_avg = [0.2612451,0.19266722,0.1908739]
Fe_std = [3.3866146,2.4976149,2.4743676]
c = ['blue','orange','green']
x = [1,2,3]
for i,ele in enumerate([K_avg,Ca_avg,Fe_avg]):
    pylab.plot(x,ele,'.--',color=c[i],markersize=25,linewidth=2)
    
    
    