# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 22:12:16 2015

@author: xfz
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 18:10:22 2015

@author: xfz
"""

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
import random
import pdb
import pickle
#import cProfile
import sys

class InteractionSpace:
    def __init__(self,initParams):
        initParams['m_array']=np.array(initParams['m_array'])
        initParams['H_array']=np.array(initParams['H_array'])
#        position=np.array(random.sample(xrange(initParams['N']),initParams['m']))
        position=np.array(range(initParams['m']))
        k=initParams['k']
        #m[0,1] entry gives number of proteins that starts with [1,0] type.
        proteinType=np.repeat(np.reshape(np.ravel(np.indices((k,k)),1),[k**2,2]),initParams['m_array'].flatten(),axis=0)
#        proteinForm=np.floor(np.random.rand(initParams['m'])*initParams['k'])
        proteinForm=np.zeros(initParams['m'])+np.lib.pad(np.ones(initParams['m']/2),(initParams['m']-initParams['m']/2,0),'constant',constant_values=0)
        self.init_state=(position,proteinType,proteinForm,initParams['m_array'])
        self.params=initParams
        self.iteration=0
        self.state=self.checkATP(self.init_state)
        self.currentE=self.energy(self.state)
    
    def energy(self,state):
        #here state could be the current state or some hypothetical state
        idx=np.argsort(state[0])
        state_sorted=(state[0][idx],state[1][idx],state[2][idx])
        E=0
        for pos in xrange(self.params['m']-1):
            #check if there is protein at the next position
            if state_sorted[0][pos+1]-state_sorted[0][pos]==1:
                E+=self.params['H_array'][state_sorted[2][pos]][state_sorted[2][pos+1]][state_sorted[1][pos][1]][state_sorted[1][pos+1][0]]
            else:
                # now since there is no protein at the next position,
                # add the two end energy 
                # end energy of the right-end domain of protein p
                E+=self.params['H_end'][state_sorted[2][pos]][state_sorted[1][pos][1]]
                # end energy of the left-end domain of protein p+1
                E+=self.params['H_end'][state_sorted[2][pos]][state_sorted[1][pos+1][0]]
        if self.params['isCircular'] and state_sorted[0][-1]==self.params['N']-1 and state_sorted[0][0]==0:
                # ends only need to be treated as linked
                # if both end positions has proteins there
                # and the scenario wanted is Circular.
                E+=self.params['H_array'][state_sorted[2][-1]][state_sorted[2][0]][state_sorted[1][-1][1]][state_sorted[1][0][0]]
        
        else:
            #Now takes into account of the two ends
            E+=self.params['H_end'][state_sorted[2][0]][state_sorted[1][0][0]]
            E+=self.params['H_end'][state_sorted[2][-1]][state_sorted[1][-1][1]]
        
        return E
    
    def change(self):
        '''
        suggests a potential change to state
        return a changed state
        '''
        ratio=self.params['transRotRatio']
        position_idx=random.choice(xrange(self.params['m']))
        if self.params['TtoDIsPhysical'] and np.random.rand()<=self.params['probFormChange'] and self.state[2][position_idx]==0:
                print('WENT THROUGH THE WRONG PLACE!!!!')                
                changedForm=self.state[2].copy()
                changedForm[position_idx]=1
                changedState=(self.state[0],self.state[1],changedForm,self.state[3])
        elif np.random.rand()<=ratio:
            #translation.
            changedPosition=self.state[0].copy()
            if self.params['isAlt']:
                #Alternative dynamics, where proteins only move if nearby position is open
                u=np.random.uniform()
                #Perform change of position
                if u<= self.params['altProb'] and (self.state[0][position_idx]-1)%self.params['N'] in self.state[0]:
                    changedPosition[position_idx]=changedPosition[position_idx]-1
                elif u> self.params['altProb'] and (self.state[0][position_idx]+1)%self.params['N'] in self.state[0]:
                    changedPosition[position_idx]=changedPosition[position_idx]+1
            else:
                #where jumping through proteins is allowed.
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
        return changedState

    
    def step(self):
        for t in xrange(self.params['timestep']):
            temp=np.random.uniform()
            if not self.params['TtoDIsPhysical']:
                self.unphysicalTtoD()
            
            self.currentE=self.energy(self.state)
            changedState=self.change()
            changedE=self.energy(changedState)
            if temp <=np.exp((self.currentE-changedE)/float(self.params['T'])):
                self.state=changedState
                self.currentE=changedE
            self.state=self.checkATP(self.state)
        self.iteration+=1

    def unphysicalTtoD(self):
        #p_h, hydroplysis probability for each ATP-form protein
        p_h=self.params['probFormChange']/3300.
#        p_h=0
        posIdx=np.array(range(self.params['m']))[self.state[2]==0]
        for idx in posIdx:
            if np.random.rand()<=p_h:
                self.state[2][idx]=1

    
    def checkATP(self,state):
        #turn all ADP form monomers into ATP
        # 0 is ATP state. 1 is ADP state
    
        for t in xrange(self.params['m']):
            if state[2][t]!=0:
                if (state[0][t]+1)%self.params['N'] not in state[0] and (state[0][t]-1)%self.params['N'] not in state[0]:
                    state[2][t]=0
        return state
    
    
    
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
    

    xcenters=intSp.state[0]+deviation/2
    ycenters=np.zeros(m)
    zcenters = np.vstack((xcenters,ycenters)).transpose()
    centers.set_offsets(zcenters)
    colorscenters = np.zeros([3,m])
    colorscenters[2]=intSp.state[2]*1/float(h-1)
    colorscenters[0]=1-colorscenters[2]
    colorscenters=colorscenters.T
    centers.set_facecolor(colorscenters)        
    
    
    time_text.set_text('time=%.1f' % float(intSp.iteration*intSp.params['timestep']))
    energy_text.set_text('energy=%.1f' % intSp.energy(intSp.state))
    return points,time_text,energy_text

def TE_simulation(fileName,initParams,T_array,xi_array,simPerPt=1,obsStart=50,obsDur=50,xi=0.05):
    filler = np.frompyfunc(lambda x: list(), 1, 1)
    energies=np.empty([len(xi_array),len(T_array),simPerPt],dtype=np.object)
    states=np.empty([len(xi_array),len(T_array),simPerPt],dtype=np.object)
    filler(energies,energies)
    filler(states,states)
    stepSizeVec=np.zeros([len(xi_array),len(T_array)])
    xiTStack=np.dstack(np.meshgrid(xi_array,T_array)).reshape(-1,2)
    stats={'initParams':initParams,'simPerPt':simPerPt,'obsStart':obsStart,'obsDur':obsDur,'xi_array':xi_array,
           'energy':energies,'temperature':T_array, 'states':states,'stepSize':stepSizeVec,'xiTStack':xiTStack}    
    for idx in xrange(len(xiTStack)):
        xi,T=xiTStack[idx]
        initParams['T']=T
        initParams['H_array'],initParams['H_end']=H_array_gen(initParams['h'],initParams['k'],xi,-1.,T)
        initParams['H_end']
        charTime=max(min(np.exp(-xi/T),1e5),10)
        stepSize=np.floor(charTime)
        stepSizeVec[idx]=stepSize
        numDataPts=np.int(np.ceil(obsDur/stepSize))
        for repeat in xrange(simPerPt):
            intSp=InteractionSpace(initParams)
            for t in xrange(obsStart): 
                intSp.step()
            for t in xrange(numDataPts):
                energies[idx][repeat].append(intSp.currentE)
                states[idx][repeat].append(intSp.state)
                stats['states']=states
                pickle.dump(stats,open(fileName,'wb'))
                for idx in np.arange(stepSize):
                    intSp.step()
    pickle.dump(stats,open(fileName,'wb'))

def H_array_gen(h=2,k=2,xi=-0.05,eps=-1.,T=1e-2):
    H_array=np.zeros([h,h,k,k])
    H_array[0][0]=np.array([[0,eps],
                           [eps,0]])
    H_array[0][1]=np.array([[0,eps],
                           [xi,0]])
    H_array[1][0]=H_array[0][1].T
    H_array[1][1]=np.array([[0,xi],
                           [xi,0]])
    H_end=np.zeros([h,k])
    return H_array, H_end


def main(fName,start,dur,hydroFlag=True):
    
    N=100
#    N=4
    
    T=1e-2
    k=2
    h=2
    
    m=20
#    m=2
    
    timestep=1
    isCircular=True
    transRotRatio=1.0
    if hydroFlag==False:
        probFormChange=0.
    else: probFormChange=0.33
    m_array=[[0,20],[0,0]]
#    m_array=[[0,2],[0,0]]
    
    #Make H_array from H_f1, H_f2, H_f1f2, H_f2f1.
    # 1st dim: left domain form, t=0, d=1
    # 2nd dim: right domain dim
    # 3rd dim: left domain type, p=0, b=1
    # 4th dim: right domain type

    
    #ising and lattice gas set-ups.
#    xi=eps
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
#    H_array=H_ising
    
    H_array=None
    H_end=None
    xi=0.5
    
    initParams={'N':N,'m':np.sum(m_array),'k':k,'h':h,'T':T,'xi':xi,'m':m,
                'timestep':timestep,'isCircular':isCircular,
                'TtoDIsPhysical':False,
                'isAlt':False,'altProb':0.5,
                'transRotRatio':transRotRatio,
                'probFormChange':probFormChange,
                'm_array':m_array,
                'H_array':H_array,
                'H_end'  :H_end}

    #write to animation
#    TE_animate('1dActinUnphysical_moreD',initParams,50,1000)

                
    #Save result to file
    # T_array=np.logspace(-3,1,num=10)
    T_array=np.array([1e-2])
    xi_array=np.logspace(1e-2,1,num=5)
    TE_simulation(fName,initParams,T_array,xi_array,simPerPt=1,obsStart=start,obsDur=dur,xi=xi)


if __name__ == "__main__":
#    cProfile.run('main()','profile.tmp')
    fName = sys.argv[1]
    start = int(sys.argv[2])
    dur = int(sys.argv[3])
    if len(sys.argv)>3:
        hydroFlag=int(sys.argv[4])   
        main(fName,start,dur,hydroFlag)
        
    else: main(fName,start,dur)