# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 21:17:18 2015

@author: xfz
"""
import numpy as np
import random
import pdb
import pickle
#import cProfile
import sys
from itertools import product
from multiprocessing import Pool


class InteractionSpace:
    def __init__(self,initParams,H):
#        initParams['m_array']=np.array(initParams['m_array'])
#        position=np.array([0,1,2,3,4,5,6,7,8,9,11,13,15,17,19,28,26,24,22,20])
        position=np.array(range(initParams['m']))
#        position=np.random.choice(range(initParams['N']),initParams['m'],replace=False)
        k=initParams['k']
        #m[0,1] entry gives number of proteins that starts with [1,0] domains.
        proteinOrientation=np.repeat(np.reshape(np.ravel(np.indices((k,k)),1),[k**2,2]),initParams['m_array'].flatten(),axis=0)
#        proteinForm=np.floor(np.random.rand(initParams['m'])*initParams['k'])
        proteinForm=np.zeros(initParams['m'])+np.lib.pad(np.ones(initParams['m']/2),(initParams['m']-initParams['m']/2,0),'constant',constant_values=0)
        self.init_state=(position,proteinOrientation,proteinForm,initParams['m_array'])
        self.params=initParams
        self.state=self.checkATP(self.init_state)
        self.H=H
        self.currentE=self.energy(self.state)
    
    def energy(self,state):
        #here state could be the current state or some hypothetical state
        idx=np.argsort(state[0])
        state_sorted=(state[0][idx],state[1][idx],state[2][idx])
        E=0
        d=int(self.params['dist_cutoff'])
        for i in xrange(self.params['m']):
            for j in xrange(d):
                rIdx=(i+j+1)%self.params['m']
                a=[state_sorted[0][i],state_sorted[2][i],state_sorted[1][i]]
                b=[state_sorted[0][rIdx],state_sorted[2][rIdx],state_sorted[1][rIdx]]
                E+=self.H(a,b)
        return E
        
    def change(self):
        '''
        suggests a potential change to state
        return a changed state
        '''
        ratio=self.params['transRotRatio']
        position_idx=int(np.random.rand()*(self.params['m']))
        noChangeFlag=1
        if np.random.rand()<=ratio:
            noChangeFlag=0
            #translation.
            changedPosition=self.state[0].copy()
            if self.params['isAlt']:
                #Alternative dynamics, where proteins only move if nearby position is open
                u=np.random.rand()
                #Perform change of position
                N=self.params['N']
                if u<= self.params['altProb'] and (self.state[0][position_idx]-1)%N not in self.state[0]:
                    changedPosition[position_idx]=(changedPosition[position_idx]-1)%N
                    
                elif u> self.params['altProb'] and (self.state[0][position_idx]+1)%N not in self.state[0]:
                    changedPosition[position_idx]=changedPosition[position_idx]+1%N
                    
            else:
                # where jumping through proteins is allowed.
                # to be more physical, we force proteins to get into empty slots
                # not connected to any other proteins before they can get back
                # into a slot next to other proteins.
                
                # flag of whether chosen protein is currently connected to 
                # other proteins
                if (self.state[0][position_idx]+1)%self.params['N'] in self.state[0] or (self.state[0][position_idx]-1)%self.params['N'] in self.state[0]:
                    # protein is currently connected to other proteins, 
                    # so we force it to go into empty slots not connected to proteins
                    connectedEmptySlots=np.union1d(np.union1d(self.state[0],(self.state[0]+1)%self.params['N']),(self.state[0]-1)%self.params['N'])
                    slotsToChoose=np.array([x for x in xrange(self.params['N']) if x not in connectedEmptySlots])
                else:
                    # if protein is currently not conencted to proteins, any slot not occupied can be taken
                    slotsToChoose = np.array([x for x in xrange(self.params['N']) if x not in self.state[0]])
                if len(slotsToChoose)!=0:
                    changedPosition[position_idx]=random.choice(slotsToChoose)
            
            changedState=(changedPosition,self.state[1],self.state[2],self.state[3])
        else:
            #Rotation
            changedType=self.state[1].copy()
            changedType[position_idx]=changedType[position_idx][::-1]
            changedM=self.state[3].copy()
            changedM[self.state[1][position_idx][1],self.state[1][position_idx][0]]-=1
            changedM[changedType[position_idx][1],changedType[position_idx][0]]+=1
            #note the above changedM, idx 1 is for row, while idx 0 is for column in the Type array. This is due to np.ravel of np.indices
            changedState=(self.state[0],changedType,self.state[2],changedM)
#        changedState=self.checkATP(changedState)
            noChangeFlag=0
        return changedState,noChangeFlag

    
    def sweep(self):
        for i in xrange(self.params['m']):
            self.step()
    
    def step(self):
        temp=np.random.rand()
        self.TtoD()
        #need to recalculate current energy because of hydrolysis
        self.currentE=self.energy(self.state)
        changedState,noChangeFlag=self.change()
        if not noChangeFlag:
            changedE=self.energy(changedState)
            if temp <=np.exp((self.currentE-changedE)/self.params['T']):
                self.state=changedState
                self.currentE=changedE
            self.state=self.checkATP(self.state)

    def TtoD(self):
        #p_h, hydroplysis probability for each ATP-form protein
        p_h=self.params['probFormChange']
        posIdx=np.array(range(self.params['m']))[self.state[2]==0]
        for idx in posIdx:
            if np.random.rand()<=p_h:
                self.state[2][idx]=1

    
    def checkATP(self,state):
        #turn all ADP form monomers into ATP
        # 0 is ATP state. 1 is ADP state
        hydro=state[2].copy()
        for t in xrange(self.params['m']):
            if state[2][t]!=0:
                if (state[0][t]+1)%self.params['N'] not in state[0] and (state[0][t]-1)%self.params['N'] not in state[0]:
                    hydro[t]=0
        return (state[0],state[1],hydro,state[3])


def TE_simulation(fileName,initParams,T_array,xi_array,eps_array,simPerPt=1,obsStart=50,nPts=50):
    filler = np.frompyfunc(lambda x: list(), 1, 1)
    numInstance=len(xi_array)*len(T_array)*len(eps_array)
    states=np.empty([numInstance,simPerPt],dtype=np.object)
    filler(states,states)

    energies=np.zeros([numInstance,simPerPt,nPts])
    stepSizeVec=np.zeros([numInstance])
    HParamStack=np.stack(np.meshgrid(T_array,eps_array,xi_array),3).reshape(-1,3)
    stats={'initParams':initParams,'simPerPt':simPerPt,'obsStart':obsStart,'nPts':nPts,'xi_array':xi_array,
           'energy':energies,'temperature':T_array, 'states':states,'stepSize':stepSizeVec,'eps_array':eps_array,
           'HParamStack':HParamStack}    
#    for idx in xrange(numInstance):
#        T,eps,xi=HParamStack[idx]
#        initParams['T']=T
#        H=lambda a,b: H_actin(a,b,eps,xi,T,initParams['N'],initParams['dist_cutoff'])
#        charTime=max(min(np.exp(-xi/T),1e3),1)
#        stepSize=np.floor(charTime)
#        stepSizeVec[idx]=stepSize
#        for repeat in xrange(simPerPt):
#            intSp=InteractionSpace(initParams,H)
#            for t in xrange(obsStart): 
#                intSp.sweep()
#            for t in xrange(nPts):
#                energies[idx][repeat][t]=intSp.currentE
#                states[idx][repeat].append(intSp.state)
#                stats['states']=states
#                for step_idx in xrange(int(stepSize)):
#                    intSp.step()
#        pickle.dump(stats,open(fileName+str(idx)+'.p','wb'))
    pool=Pool()
    chunksize=1
    iterList=[[x,initParams,simPerPt,obsStart,nPts] for x in HParamStack]
#    for ind, res in enumerate(pool.imap(parFun, iterList),chunksize):
#        energies[ind-1],states[ind-1],stepSizeVec[ind-1]=res
    parFun(iterList[0])
    stats['energy']=energies
    stats['states']=states
    stats['stepSizeVec']=stepSizeVec
    pickle.dump(stats,open(fileName,'wb'))

def parFun(pList):
    param,initParams,simPerPt,obsStart,nPts=pList
    energy=np.zeros([simPerPt,nPts])
    filler = np.frompyfunc(lambda x: list(), 1, 1)
    state=np.empty([simPerPt],dtype=np.object)
    filler(state,state)
    
    T,eps,xi=param
    params=initParams
    params['T']=T
    H=lambda a,b: H_actin(a,b,eps,xi,T,initParams['N'],initParams['dist_cutoff'])
    charTime=max(min(np.exp(-xi/T),1e3),1)
    stepSize=np.floor(charTime)
    for repeat in xrange(simPerPt):
        intSp=InteractionSpace(params,H)
        for t in xrange(obsStart): 
            intSp.sweep()
        for t in xrange(nPts):
            energy[repeat][t]=intSp.currentE
            state[repeat].append(intSp.state)
            for step_idx in xrange(int(stepSize)):
                intSp.step()
    return energy,state,stepSize

def H_actin(a,b,eps,xi,T,N,dist_cutoff=1):
    # a=(pos,form, domain type, protein type)
    # form: 1 is ADP, 0 is ATP.
    # ORDER MATTERS. (a,b) is an ORDERED pair, meaning protein a is on the left
    #   side of interaction. This is important when distance is calculated.

    h=k=2
    H_array=np.zeros([h,h,k,k])
    H_array[0][0]=np.array([[0,eps],
                           [eps,0]])
    H_array[0][1]=np.array([[0,eps],
                           [xi,0]])
    H_array[1][0]=H_array[0][1].T
    H_array[1][1]=np.array([[0,xi],
                           [xi,0]])


    #Make H_array from H_f1, H_f2, H_f1f2, H_f2f1.
    # 1st dim: left domain form, t=0, d=1
    # 2nd dim: right domain dim
    # 3rd dim: left domain type, p=0, b=1
    # 4th dim: right domain type

    
    #ising and lattice gas set-ups.
#    eps=xi
#    H_ising=np.zeros([h,h,k,k])
#    H_ising[0][0]=np.array([[-eps,eps],
#                            [eps,-eps]])
#    H_ising[0][1]=H_ising[0][0]
#    H_ising[1][0]=H_ising[0][0]
#    H_ising[1][1]=H_ising[0][0]
#    
#    H_gas=np.zeros([h,h,k,k])
#    H_gas[0][0]=np.array([[eps,eps],
#                            [eps,eps]])
#    H_gas[0][1]=H_gas[0][0]
#    H_gas[1][0]=H_gas[0][0]
#    H_gas[1][1]=H_gas[0][0]    
#    
#    H_array=H_gas

    dist=(b[0]-a[0])%N
    if np.abs(dist)<=dist_cutoff:
        return H_array[a[1]][b[1]][a[2][1]][b[2][0]]/np.abs(float(dist))
    else: return 0.
        

    


def main(fName,start,nPts,hydroFlag=True):
    
    N=30
    k=2
    h=2
    m_array=np.array([[0,20],[0,0]])
    
    
    if hydroFlag==False:
        probFormChange=0.
    else: probFormChange=1e-3

    # T, xi, eps, to be filled in.
    initParams={'N':N,'m':np.sum(m_array),'k':k,'h':h,
                'T':None,'xi':None, 'eps':None,
                'isCircular':True,
                'isAlt':False,'altProb':0.5,
                'transRotRatio':1.,
                'probFormChange':probFormChange,
                'm_array':m_array,
                'dist_cutoff':2.}

    #write to animation
#    TE_animate('1dActin',initParams,50,1000)

                
    T_array=np.array([1.])
    xi_array=-1*np.linspace(1,8,num=4)
    eps_array=-1*np.linspace(2,100,num=10)
#    xi_array=np.array([-5.])
#    eps_array=np.array([-100.])
    TE_simulation(fName,initParams,T_array,xi_array,eps_array,simPerPt=1,obsStart=start,nPts=nPts)



if __name__ == "__main__":
#    cProfile.run('main()','profile.tmp')
#    python OneD_actin.py asdf.p 150 100
    fName = sys.argv[1]
    start = int(sys.argv[2])
    nPts = int(sys.argv[3])
    if len(sys.argv)>4:
        hydroFlag=int(sys.argv[4])
        main(fName,start,nPts,hydroFlag)
        
    else: main(fName,start,nPts)