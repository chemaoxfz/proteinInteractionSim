# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 07:31:26 2015

@author: xfz
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 23:59:47 2015

@author: xfz

Simulation of 1-D protein interaction.

The rule is currently as follows:
0 Provided
0.1 n, number of slots in this universe
0.2 m_11, m_12,..., m_kk, total k(k+1)/2 parameters,
this specifies number of proteins of each type.
0.3 H_ij, i,j range from 0 to k-1. k is the number of domain types.
0.4 temperature T.

1 Setup
1.1 Create the universe (n slots)
1.2 Create the protein portfolio (m_11, ..., m_kk number of proteins
of each type of proteins)
1.3 fill the universe with m= SUM m_i proteins, randomly.

2 Sample the state space for a given set of proteins,
which has parameters of
x_1,...,x_m position of each protein
type_1,...,type_m, types of each protein


2.1 randomly choose a protein, move to next spot or rotates
with probability proportioinal to
exp(-(E'-E)/T).

where T is some temperature parameter.
E' is sum of interaction energy according to H for after the change
E is sum of interaction energy according to H before the change.


Potential generalizations:
1. explore m_array space (i.e. changing the 
bunch of proteins given in the simulation 
process)

2. explore rotation space

3. use some fitness function other than E.
"""
import numpy as np
import matplotlib.pyplot as plt
import random
import pdb
import matplotlib.animation as animation
import pickle

class InteractionSpace:
    def __init__(self,
                 initParams={'n':50,
                 'T':1,
                 'k':2,
                 'm':45,
                 'timestep':10,
                 'isCircular':False,
                 'isAlt':False,
                 'altProb':0.5,
                 'transRotRatio':0.5,
                 'm_array':[[20,10],[10,5]],
                 'H_array':[[2,1,0],[1,3,0],[0,0,0]]}):
        initParams['m_array']=np.array(initParams['m_array'])
        initParams['H_array']=np.array(initParams['H_array'])
        position=np.array(random.sample(xrange(n),initParams['m']))
        k=initParams['k']
        proteinType=np.repeat(np.reshape(np.ravel(np.indices((k,k)),1),[k**2,2]),initParams['m_array'].flatten(),axis=0)
        self.init_state=(position,proteinType,initParams['m_array'])
        self.params=initParams
        self.iteration=0
        self.state=self.init_state
        self.currentE=self.energy(self.state)
    
    def energy(self,state):
        #here state could be the current state or some hypothetical state
        idx=np.argsort(state[0])
        state_sorted=(state[0][idx],state[1][idx])
        E=0
        for p in np.arange(len(state_sorted[0])-1):
            #check if there is protein at the next position
            if state_sorted[0][p+1]-state_sorted[0][p]==1:
                E+=self.params['H_array'][state_sorted[1][p][1]][state_sorted[1][p+1][0]]
            else:
                # now since there is no protein at the next position,
                # add the two end energy 
                # end energy of the right-end domain of protein p
                E+=self.params['H_array'][state_sorted[1][p][1]][self.params['k']]
                # end energy of the left-end domain of protein p+1
                E+=self.params['H_array'][state_sorted[1][p+1][0]][self.params['k']]
        if self.params['isCircular'] and state_sorted[0][-1]==self.params['n']-1 and state_sorted[0][0]==0:
                # ends only need to be treated as linked
                # if both end positions has proteins there
                # and the scenario wanted is Circular.
                E+=self.params['H_array'][state_sorted[1][-1][1]][state_sorted[1][0][0]]
        
        else:
            #Now takes into account of the two ends
            E+=self.params['H_array'][state_sorted[1][0][0]][self.params['k']]
            E+=self.params['H_array'][state_sorted[1][-1][1]][self.params['k']]
        
        
        return E
    
    def change(self,altProb=0.5):
        '''
        suggests a potential change to state
        return a changed state
        '''
        ratio=self.params['transRotRatio']
        position_idx=random.choice(xrange(self.params['m']))
        if np.random.rand()<=ratio:
            #translation.
            changedPosition=self.state[0].copy()
            if self.params['isAlt']:
                #Alternative dynamics, where proteins only move if nearby position is open
                u=np.random.uniform()
                if u<= self.params['altProb'] and self.state[0][position_idx]-1 in self.state[0]:
                    changedPosition[position_idx]=changedPosition[position_idx]-1
                elif u> self.params['altProb'] and self.state[0][position_idx]+1 in self.state[0]:
                    changedPosition[position_idx]=changedPosition[position_idx]+1
            else:
                #where jumping through proteins is allowed.
                emptySlots = np.array([x for x in xrange(n) if x not in self.state[0]])
                changedPosition[position_idx]=random.choice(emptySlots)
            
            changedState=(changedPosition,self.state[1],self.state[2])
        else:
            changedType=self.state[1].copy()
            changedType[position_idx]=changedType[position_idx][::-1]
            changedM=self.state[2].copy()
            changedM[self.state[1][position_idx][1],self.state[1][position_idx][0]]-=1
            changedM[changedType[position_idx][1],changedType[position_idx][0]]+=1
            #note the above changedM, idx 1 is for row, while idx 0 is for column in the Type array. This is due to np.ravel of np.indices
            changedState=(self.state[0],changedType,changedM)
        return changedState

    
    def step(self):
        for t in xrange(self.params['timestep']):
            temp=np.random.uniform()
            changedState=self.change()
            changedE=self.energy(changedState)
            if temp <=np.exp((self.currentE-changedE)/float(self.params['T'])):
                self.state=changedState
                self.currentE=changedE
        self.iteration+=1

    
def update_plot(i):
    global intSp
    intSp.step()
    x = np.tile(intSp.state[0],2)+np.concatenate((np.zeros(m),np.ones(m)*deviation))
    x = np.ravel(np.reshape(x,[2,m]),1)
    y = np.zeros(2*intSp.params['m'])
    z = np.vstack((x,y)).transpose()
    points.set_offsets(z)
    colors = np.lib.pad(np.repeat(np.reshape(intSp.state[1].flatten()*1/float(k-1),[2*m,1]),2,axis=1),(0,1),'constant',constant_values=(0))
    points.set_facecolor(colors)
#    points.set_lw(2)
    time_text.set_text('time=%.1f' % float(intSp.iteration*intSp.params['timestep']))
    energy_text.set_text('energy=%.1f' % intSp.energy(intSp.state))
    return points,time_text,energy_text

def TE_simulation(fileName,initParams,T_array,simPerPt=1,obsStart=50,obsDur=50):
    energies=np.zeros([len(T_array),simPerPt,obsDur])
    for idx_T in range(len(T_array)):
        T=T_array[idx_T]
        initParams['T']=T
        for repeat in xrange(simPerPt):
            intSp=InteractionSpace(initParams)
            for t in xrange(obsStart): 
                intSp.step()
            for t in xrange(obsDur):
                energies[idx_T][repeat][t]=intSp.currentE
                intSp.step()
    stats={'energy':energies,'temperature':T_array}    
    pickle.dump(stats,open(fileName,'wb'))

 
if __name__ == "__main__":
    
    n=2
    T=1e1
    k=2
    m=1
    timestep=100
    isCircular=True
    transRotRatio=0.5
    m_array=[[0,0],[1,0]]
    H_array=[[-1,-1,0],
             [-1,-1,0],
             [0,0,0]]
    initParams={'n':n,'m':np.sum(m_array),'k':k,'T':T,'m':m,
                'timestep':timestep,'isCircular':isCircular,
                'isAlt':False,'altProb':0.5,
                'transRotRatio':transRotRatio,
                'm_array':m_array,
                'H_array':H_array}
    intSp=InteractionSpace(initParams)
    
    fig=plt.figure()
    ax=fig.add_subplot(111,autoscale_on=False,xlim=(-1,n),ylim=(-n/5.,n/5.))
    deviation=0.3 #deviation is the distance between the center of right domain from left domain of the same protein
    x = np.tile(intSp.state[0],2)+np.concatenate((np.zeros(m),np.ones(m)*deviation))
    x = np.ravel(np.reshape(x,[2,m]),1)
    y = np.zeros(2*m)
    z = np.vstack((x,y)).transpose()
    colors = np.lib.pad(np.repeat(np.reshape(intSp.state[1].flatten()*1/float(k-1),[2*m,1]),2,axis=1),(0,1),'constant',constant_values=(0))
    points=ax.scatter(x,y,s=np.pi*30,c=colors)
    time_text=ax.text(0.02,0.95,'',transform=ax.transAxes)
    energy_text=ax.text(0.02,0.90,'',transform=ax.transAxes)
    ani=animation.FuncAnimation(fig,update_plot,
                                interval=100,blit=False)
    plt.show()
#    T_array=np.logspace(-5,1,num=20)
#    TE_simulation('save.p',initParams,T_array,simPerPt=1,obsStart=50,obsDur=50)