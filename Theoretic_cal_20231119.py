# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 21:22:23 2023

@author: lxiaoyang
"""

import os
import xraydb
import numpy as np
import math
#%%
inci_eng = 10000 #ev
ele_AXO = ['Pb','La','Pd','Mo','Cu','Fe']
wr = np.divide([21.82320442,50,4.972375691,2.762430939,6.35359116,14.08839779],100)
ele_siN = ['Si','N']
FeK = xraydb.fluor_yield('Fe','K','Ka',inci_eng)[1] #eV
CuK = xraydb.fluor_yield('Cu','K','Ka',inci_eng)[1]
thick = 3e-7 #cm
ang_pos = [0,10,15,20,30,40,45,50,60,70,80,89]
#ang_pos = [7.45,2,4,7,10,13,19,22,25,28,31,37,40,43,46,49,52,55,58,61,64,34,67]
ang_neg = [0,-10,-15,-20,-30,-40,-45,-50,-60,-70,-80,-89]
#ang_neg = [-31,-37,-40,-43,-46,-49,-52,-58,-64,-70,-34,-55,-67,1,-3,-6,-9,-61,-76,-79,-82,-11,-14,-73]  #change 0 to 1

#%%
#calculate xray mass attenuation coefficient (total atten coeff) from Elam tables
def cal_mu_tot(e,eng):
    ei_mu_list = []
    for ei in e:
        ei_mu = xraydb.mu_elam(ei,eng,kind='total') #cm2/gr, mass attenuation coefficient
        ei_mu_list.append(ei_mu)
    return ei_mu_list
#%%
#cal total attenation coefficient for axo at energy
def cal_muAXO(ele,wr,eng):
    ei_mu = cal_mu_tot(ele,eng)
    v = np.sum(np.multiply(ei_mu,wr))
    return v
#%%
#cal based on Optimizing... paper
#postive ang respect to 15 deg
ang_ref = 15
fluo_eng = FeK
u = 1-math.exp(-thick*((cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng)/math.sin(math.radians(90-ang_ref)))+(cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng)/math.sin(math.radians(ang_ref)))))
d = cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng)+(np.divide(math.sin(math.radians(90-ang_ref)),math.sin(math.radians(ang_ref))))*cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng)
int_ref = np.divide(u,d)  
ratio_list = []
for ang in ang_pos:
    u = 1-math.exp(-thick*((cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng)/math.sin(math.radians(90-ang)))+(cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng)/math.sin(math.radians(ang)))))
    d = cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng)+(np.divide(math.sin(math.radians(90-ang)),math.sin(math.radians(ang))))*cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng)
    int_ang = np.divide(u,d)  
    r = int_ang / int_ref
    ratio_list.append(r)
    print(f'{ang} done,ratio {r}')
#%%
#cal based on Optimizing... paper
#negative ang respect to 15 deg
#rotatio angle is phi', phi is 90-rotation angle
ratio_neg_list = []
ang_ref = 15
fluo_eng = CuK
u = 1-math.exp(-thick*((cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng)/math.sin(math.radians(90-ang_ref)))+(cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng)/math.sin(math.radians(ang_ref)))))
d = cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng)+(np.divide(math.sin(math.radians(90-ang_ref)),math.sin(math.radians(ang_ref))))*cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng)
int_ref = np.divide(u,d)  
for ang in ang_neg:
    u = math.exp(np.divide(cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng),math.sin(math.radians(-ang)))*thick)-math.exp(-np.divide(cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng),math.sin(math.radians(90+ang)))*thick)     
    d = cal_muAXO(ele=ele_AXO,wr=wr,eng=inci_eng)+(np.divide(math.sin(math.radians(90+ang)),math.sin(math.radians(-ang))))*cal_muAXO(ele=ele_AXO,wr=wr,eng=fluo_eng)
    int_ang_neg = np.divide(u,d)
    r = np.divide(int_ang_neg,int_ref)
    ratio_neg_list.append(r)
    print(f'{ang}done,ratio{r}')