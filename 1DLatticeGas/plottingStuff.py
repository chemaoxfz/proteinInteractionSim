# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:56:05 2015

@author: xfz
"""
import numpy as np
import matplotlib.pyplot as plt
import random
import pdb
import matplotlib.animation as animation
import pickle


from scipy.special import comb

def plotFromFile(fileName,figName):  
    stats=pickle.load(open(fileName,'rb'))
    energies=stats['energy']
    T_array=stats['temperature']
    EVector=np.mean(np.mean(energies,axis=2),axis=1)
    EStd=np.std(np.mean(energies,axis=1),axis=1)
    fig, ax = plt.subplots(1)
    ax.plot(1./T_array, EVector,'-bo', lw=2, markersize = 10, label='simulation')
    ax.plot(1./T_array,theoryE(1./T_array),'-ko',lw=2,label='Exact 1D Gas Theory')
    
#    ax.plot(1./T_array,-20*(np.tanh(1./T_array)),'-ro',lw=2,label='Large N 1D Ising Theory')
    ax.fill_between(1./T_array, EVector+EStd, EVector-EStd, facecolor='blue',alpha=0.5)
    ax.set_title('1D Lattice Gas Model')
    ax.legend(loc='upper left')
    ax.set_xlabel('1/T')
    ax.set_ylabel('E')
    ax.set_xscale('log')
    plt.savefig(figName)
    plt.show()
    

def theoryE(beta,epsilon=1,m=2,N=4):
    if N==m: return -N*np.ones(len(beta))
    else:
        kmax = min(N-m+1,m)
        E=np.zeros(len(beta))
        temp1=np.zeros([kmax, len(beta)])
        temp2=np.zeros([kmax, len(beta)])
        for t in xrange(kmax):
            k=t+np.float64(1.)
            temp=1./k * comb(m-1,k-1)*comb(N-m-1,k-1)*np.exp(-k*epsilon*beta)
            temp1[t]=temp*(m-k)*epsilon
            temp2[t]=temp
        tempSum1=np.sum(temp1,axis=0)
        tempSum2=np.sum(temp2,axis=0)
    #    tempSum2+=np.ones(len(beta))*1e-10 #for numerical stability
    #    tempSum1+=np.ones(len(beta))*9e-10
    #    pdb.set_trace()
        E=- tempSum1/tempSum2
    return E


if __name__ == "__main__":
    fileName='1DGasTestN4m2_1e3'
    plotFromFile(fileName+'.p',fileName+'.pdf')