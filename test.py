# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:32:09 2015

@author: xfz
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
stats=pickle.load(open('save500.p','rb'))
energies=stats['energy']
T_array=stats['temperature']
EVector=np.mean(np.mean(energies,axis=2),axis=1)
EStd=np.std(np.mean(energies,axis=1),axis=1)
# plot it!
fig, ax = plt.subplots(1)
ax.plot(1./T_array,EVector, 'bo-',lw=2, label='mean population 1', color='blue')
ax.plot(1./T_array,-20*(np.tanh(1./T_array)+1./np.tanh(1./T_array)),'-ro',lw=2,label='exact theory')

ax.fill_between(1./T_array, EVector+EStd, EVector-EStd, facecolor='blue',alpha=0.5)
ax.set_xscale('log')
