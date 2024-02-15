# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 13:40:15 2024

@author: lxiaoyang
"""

#test iterative method for ice thickness estimation
#%%
from skimage import io
import os
import matplotlib.pylab as pylab
import numpy as np
import math
import glob
import h5py
import xraydb
from scipy.optimize import minimize,fsolve
#%%
#get ug/cm2 values from h50 file
path = r'\\micdata\data1\bnp\2023-1\simu_test_bio_varyice\img.dat'
os.chdir(path)

fn = 'BIO5um_ice30um.mda.h50'
bio = h5py.File(fn,'r')
a = bio['MAPS/Spectra/mca_arr'][:]
cali = bio['MAPS/Quantification/Calibration/NNLS/Calibration_Curve_US_IC'][:][0] #calibration curve in K edge
cali_nm = bio['MAPS/Quantification/Calibration/NNLS/Calibration_Curve_Labels'][:][0]
US_IC = bio['MAPS/Quantification/Standard0/Scalers/US_IC'][:][0]

bio.close()

#%%
#get xray mass attenuation coefficient (total atten coeff) for each element at certain energy from Elam tables
def mu_tot_e(e,eng): #for list of elements
    ei_mu_list = []
    for ei in e:
        ei_mu = xraydb.mu_elam(ei,eng,kind='total') #cm2/gr, mass attenuation coefficient
        ei_mu_list.append(ei_mu)
    return ei_mu_list
#%%
#cal total attenation coefficient for the sample at certain energy
def cal_mu_tot(ele,wr,eng):
    ei_mu = mu_tot_e(ele,eng)
    v = np.sum(np.multiply(ei_mu,wr))
    return v
#%%
#get xrf energy
def xrf_eng(e,edge,line,inc_eng):
    fy = xraydb.fluor_yield(e,edge,line,inc_eng) #fluorescence yield, weighted average fluorescence energy, net_probability
    return fy[1]
#%%
#get fluorescence yield of a set of elements
def fluo_y(e,edge,line,inc_eng):
    fy = xraydb.fluor_yield(e,edge,line,inc_eng)  #fluorescence yield, weighted average fluorescence energy, net_probability
    return fy[0]   
#%%
#calculate xrf intensity, forward geometry, no ice layer
def cal_xrf_int(I0, e, xrf_y, in_ang, out_ang, solid_ang, mu_tot_in, mu_tot_xrf, t):
    print(f'Calculate {e}')
    up = I0*xrf_y*solid_ang
    down = mu_tot_in + np.divide(math.sin(math.radians(in_ang)),math.sin(math.radians(out_ang)))*mu_tot_xrf
    a = np.divide(up,down)
    xrf_int_noice = a*(math.exp(np.divide(mu_tot_xrf,math.sin(math.radians(out_ang)))*t)-math.exp(-np.divide(mu_tot_in,math.sin(math.radians(in_ang)))*t))
    return xrf_int_noice
#%%
#optimized sample thickness from xrf intensity calculation, no ice layer
def cal_sp_t(e,real_xrf_int,I0, xrf_y, in_ang, out_ang, solid_ang, mu_tot_in, mu_tot_xrf,initial_t):
    #objective function to min the difference bewteen the result and target
    objective_function = lambda t: abs(cal_xrf_int(I0, e, xrf_y, in_ang, out_ang, solid_ang, mu_tot_in, mu_tot_xrf, t) - real_xrf_int)
    #optimization to find the value of t
    result = fsolve(objective_function, initial_t)
    #extract optimized value
    #optimized_t = result.x[0]
    print(f'optimized sample thickness: {result[0]} based on {e}')
    return result[0]
#%%
def cal_sp_I0(e,real_xrf_int,t, xrf_y, in_ang, out_ang, solid_ang, mu_tot_in, mu_tot_xrf,initial_I0,method):
    #objective function to min the difference bewteen the result and target
    objective_function = lambda I0: abs(cal_xrf_int(I0, e, xrf_y, in_ang, out_ang, solid_ang, mu_tot_in, mu_tot_xrf, t) - real_xrf_int)
    #optimization to find the value of t
    result = minimize(objective_function, initial_I0,method=method)
    #extract optimized value
    optimized_I0 = result.x[0]
    print(f'optimized I0: {optimized_I0} based on {e}')
    return optimized_I0
#%%
#calculate weight ratio for all elements from intensity
def cal_wr(e,xrf_int,US_IC,cali):  #ug/cm2 = xrf_int/cali/US_IC from Arthur's paper
    ug_cm2_list = []
    for i,ei in enumerate(e):
        ug_cm2 = np.divide(np.multiply(xrf_int[i],US_IC),cali[i])
        ug_cm2_list.append(ug_cm2)
    wt_percent = [ug / sum(ug_cm2_list) for ug in ug_cm2_list]
    return wt_percent
#%%
#validation sample thickness t
def val_t(e,I0,xrf_y,solid_ang,mu_tot_in,mu_tot_xrf,in_ang,out_ang,real_int):
    up = I0*xrf_y*solid_ang
    down = mu_tot_in + np.divide(math.sin(math.radians(in_ang)),math.sin(math.radians(out_ang)))*mu_tot_xrf
    a = np.divide(up,down)
    aa = np.divide(real_int,a)
    return aa
#%%
#calculate ice attenuation (assume ice at the backside of the sample)
def atten_ice(e,out_ang,t_ice,mu_tot_e):
    A_out = math.exp(-np.divide(mu_tot_e,math.sin(math.radians(out_ang)))*t_ice)
    return A_out
    

        
#%%
#true weight fraction for each element in simulation
vali = {'H':0.09027,
        'C':0.5015,
        'N':0.12036,
        'O':0.22066,
        'Na':0.01003,
        'Mg':0.00502,
        'P':0.02006,
        'S':0.01003,
        'Cl':0.00502,
        'K':0.01003,
        'Ca':0.00502,
        'Fe':0.00201}
#weight ratio without H, C, N, O, Na, Mg
vali2 = {'P':0.3845,
         'S':0.1923,
         'Cl':0.0962,
         'K':0.1923,
         'Ca':0.0962,
         'Fe':0.0385}
#%%
#input parameters #use forward geometry
#I0 = 10000000000 #photons/s
I0 = 302251131.0389207 #calculated from true weight ratio, thickness, no consider of Be, Si, need to figure out the incident beam or add the effect of Be and Si or others
inc_eng = 10000 #eV
out_ang = 15 #deg
in_ang = 90-15 #deg
d_det = 1.5 #cm, detector distance
a_det = 0.5 #cm2, detector active area
solid_ang = np.divide(0.5,4*math.pi*1.5*1.5) #0.01768
real_sp_t = 0.0005 #cm
ele = ['P','S','Cl','K','Ca','Fe'] #elements can be measured by xray
int_xrf = [78.2926, 106.972, 110.5, 672.188, 514.982, 990.782] #xrf intensity, NNLS
ug_per_cm2 = [17.4023,13.9059,9.48069,25.1849,14.0776,6.85693] #read from UprobeX, NNLS-US_IC
cali = [3.209938305726523e-05,5.488493631665106e-05,8.315814574643702e-05,
        0.0001904286529741502,0.00026100285815410657,0.001030932709108794]
wt_percent = [ug / sum(ug_per_cm2) for ug in ug_per_cm2]  #initialize weight ratio
#initial wt_percent:
    #[0.20023744562085655,
     #0.16000654482792903,
     #0.10908840488459562,
     #0.2897869847213707,
     #0.1619821899675428,
     #0.07889842997770524]
US_IC = 140157.96784607437
Fe_xrf_eng = xrf_eng('Fe','K','Ka',inc_eng)
P_xrf_eng = xrf_eng('P','K','Ka',inc_eng)
S_xrf_eng = xrf_eng('S','K','Ka',inc_eng)
Cl_xrf_eng = xrf_eng('Cl','K','Ka',inc_eng)
K_xrf_eng = xrf_eng('K','K','Ka',inc_eng)
Ca_xrf_eng = xrf_eng('Ca','K','Ka',inc_eng)

mu_tot_inc = cal_mu_tot(ele=ele,wr=wt_percent,eng=inc_eng)
mu_tot_Fek = cal_mu_tot(ele=ele,wr=wt_percent,eng=Fe_xrf_eng)
Fe_xrf_y = fluo_y('Fe','K','Ka',inc_eng)
P_xrf_y = fluo_y('P','K','Ka',inc_eng)
S_xrf_y = fluo_y('S','K','Ka',inc_eng)
Cl_xrf_y = fluo_y('Cl','K','Ka',inc_eng)
Ca_xrf_y = fluo_y('K','K','Ka',inc_eng)
Ca_xrf_y = fluo_y('Ca','K','Ka',inc_eng)

mu_tot_Pk = cal_mu_tot(ele=ele,wr=wt_percent,eng=P_xrf_eng)
ele_full = ['H','C','N','O','Na','Mg','P', 'S', 'Cl', 'K', 'Ca','Fe']
mu_tot_inc = cal_mu_tot(ele=ele,wr=list(vali2.values()),eng=inc_eng)
mu_tot_Fek = cal_mu_tot(ele=ele,wr=list(vali2.values()),eng=Fe_xrf_eng)

#%%
#get t_initial
t_ini = cal_sp_t(e='Fe',real_xrf_int=990.7821,I0=I0, xrf_y=Fe_xrf_y, in_ang=in_ang, out_ang=out_ang, solid_ang=solid_ang, mu_tot_in=mu_tot_inc, mu_tot_xrf=mu_tot_Fek,initial_t=0.001)
#%%
#calculate xrf intensity for other elements
xrf_int_list = []
wr_list = []
t = 0.0005
#wr_percent = [0.2770372930947969, 0.2573910275260115, 0.22462358284414186, 0.10774806636894253, 0.10364224148693205, 0.029557788679175124]
for i,ei in enumerate(ele):
    if ei != 'Fe':
        ei_xrf_eng = xrf_eng(ei,'K','Ka',inc_eng)
        mu_tot_inc = cal_mu_tot(ele=ele,wr=wt_percent,eng=inc_eng)
        mu_tot_eik = cal_mu_tot(ele=ele,wr=wt_percent,eng=ei_xrf_eng)
        ei_xrf_y = fluo_y(ei,'K','Ka',inc_eng)
        ei_int = cal_xrf_int(I0=I0, e=ei, xrf_y=ei_xrf_y, in_ang=in_ang, out_ang=out_ang, solid_ang=solid_ang, mu_tot_in=mu_tot_inc, mu_tot_xrf=mu_tot_eik, t=t)
        xrf_int_list.append(ei_int)
    else:
        xrf_int_list.append(int_xrf[i])
wr_update = cal_wr(e=ele,xrf_int=xrf_int_list,US_IC=US_IC,cali=cali)
print(wr_update)

mu_tot_inc = cal_mu_tot(ele=ele,wr=wr_update,eng=inc_eng)
mu_tot_Fek = cal_mu_tot(ele=ele,wr=wr_update,eng=Fe_xrf_eng)
#t_update = cal_sp_t(e='Fe',real_xrf_int=990.7821,I0=I0, xrf_y=Fe_xrf_y, in_ang=in_ang, out_ang=out_ang, solid_ang=solid_ang, mu_tot_in=mu_tot_inc, mu_tot_xrf=mu_tot_Fek,initial_t=t)
#%%
#if we know the thickness of sample, step size of t_ice: 0.5cm
step_t = 0.001 #cm

Fe_xrf_eng = xrf_eng('Fe','K','Ka',inc_eng)  #not use
P_xrf_eng = xrf_eng('P','K','Ka',inc_eng)
S_xrf_eng = xrf_eng('S','K','Ka',inc_eng)
Cl_xrf_eng = xrf_eng('Cl','K','Ka',inc_eng)
K_xrf_eng = xrf_eng('K','K','Ka',inc_eng)
Ca_xrf_eng = xrf_eng('Ca','K','Ka',inc_eng)

e_w = ['H','O']
e_w_wr = [0.11189,0.88811]

mu_ice_P = cal_mu_tot(ele=e_w,wr=e_w_wr,eng=P_xrf_eng)
mu_ice_S = cal_mu_tot(ele=e_w,wr=e_w_wr,eng=S_xrf_eng)
mu_ice_Cl = cal_mu_tot(ele=e_w,wr=e_w_wr,eng=Cl_xrf_eng)
mu_ice_K = cal_mu_tot(ele=e_w,wr=e_w_wr,eng=K_xrf_eng)
mu_ice_Ca = cal_mu_tot(ele=e_w,wr=e_w_wr,eng=Ca_xrf_eng)
mu_ice = [mu_ice_P,mu_ice_S,mu_ice_Cl,mu_ice_K,mu_ice_Ca]

ele = ['P','S','Cl','K','Ca','Fe'] #elements can be measured by xray

cali = [3.209938305726523e-05,5.488493631665106e-05,8.315814574643702e-05,
        0.0001904286529741502,0.00026100285815410657,0.001030932709108794]

int_xrf = [78.2926, 106.972, 110.5, 672.188, 514.982, 990.782] #xrf intensity, NNLS
US_IC = 140157.96784607437
I0 = 302251131.0389207 #calculated from true weight ratio, thickness, no consider of Be, Si, need to figure out the incident beam or add the effect of Be and Si or others
inc_eng = 10000 #eV
out_ang = 15 #deg
in_ang = 90-15 #deg
d_det = 1.5 #cm, detector distance
a_det = 0.5 #cm2, detector active area
solid_ang = np.divide(0.5,4*math.pi*1.5*1.5) #0.01768
real_sp_t = 0.0005 #cm
Fe_xrf_y = fluo_y('Fe','K','Ka',inc_eng)
fig = pylab.figure(figsize=(7,7))
ax = fig.subplots()
Fe_xrf_update = 1159.9468
for i in range(31):
    Fe_xrf_old = Fe_xrf_update
    t_ice_update = i*step_t
    att_e_int_update = []
    for ii,e in enumerate(ele):
        if e != 'Fe':
            e_int_update = np.divide(int_xrf[ii], atten_ice(e=e,out_ang=out_ang,t_ice=t_ice_update,mu_tot_e=mu_ice[ii]))
            att_e_int_update.append(e_int_update)
        else:
            att_e_int_update.append(1159.9468)
    print(att_e_int_update)
    wr_update = cal_wr(e=ele,xrf_int=att_e_int_update,US_IC=US_IC,cali=cali)
    print(wr_update)
    mu_tot_inc = cal_mu_tot(ele=ele,wr=wr_update,eng=inc_eng)
    mu_tot_Fe = cal_mu_tot(ele=ele,wr=wr_update,eng=Fe_xrf_eng)
    Fe_xrf_update = cal_xrf_int(I0=I0, e='Fe', xrf_y=Fe_xrf_y, in_ang=in_ang, out_ang=out_ang, solid_ang=solid_ang, 
                mu_tot_in=mu_tot_inc, mu_tot_xrf=mu_tot_Fe, t=real_sp_t) 
    #diff_1 = np.divide(np.abs(Fe_xrf_update-Fe_xrf_old),Fe_xrf_old)
    diff = np.divide(np.abs(np.subtract(Fe_xrf_update,1159.94682)),1159.9468)
    print(f'loop: {i}; Updated Fe XRF: {Fe_xrf_update}; ice thick: {t_ice_update}')
    pylab.plot(Fe_xrf_update)
    if diff < 0.01:# and Fe_xrf_update < 990.782:
        print(f'loop{i}: ice_thick is {t_ice_update}, cal Fe XRF: {Fe_xrf_update}')
        break
    else:
        continue
                                 
    
