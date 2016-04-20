# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 01:55:18 2015

@author: chemaoxfz
"""

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


class Unit:
    def __init__(self,dic):
        self.pos=dic['pos']
        self.domains=dic['domains']
        self.state=dic['state']
        
class InteractionSpace:
    def __init__(self,initParams,H):
#        initParams['m_array']=np.array(initParams['m_array'])
        position=np.array([0,1,2,3,4,5,6,7,8,9,11,13,15,19,20,22,25,27,28,29])
#        position=np.array(range(initParams['m']))
#        position=np.random.choice(range(initParams['N']),initParams['m'],replace=False)
        k=initParams['k']
        #m[0,1] entry gives number of proteins that starts with [1,0] domains.
        proteinOrientation=np.repeat(np.reshape(np.ravel(np.indices((k,k)),1),[k**2,2]),initParams['m_array'].flatten(),axis=0)
#        proteinForm=np.floor(np.random.rand(initParams['m'])*initParams['k'])
        proteinForm=np.zeros(initParams['m'])+np.lib.pad(np.ones(initParams['m']/2),(initParams['m']-initParams['m']/2,0),'constant',constant_values=0)
    

        self.init_state=[position,proteinOrientation,proteinForm]
        self.params=initParams
        self.state=self.checkATP(self.init_state)
        self.newState=[x.copy() for x in self.state]
        self.newStateIdx=np.arange(self.params['m']) #the idx s.t. protein identity satisfy newState[idx]=state
        self.H=H
#        self.currentE=self.energy(self.state)
        
        self.nonCyclicPos=position
    
        self.al,self.ald=self.alConstruct()    
        self.event=None
        self.time=0.
    
    def getIdx(self,pos,statePos=None):
        if statePos==None:
            statePos=self.state[0]
        idxFun=lambda pos: np.searchsorted(statePos,pos)
        return map(idxFun,pos)
    
    #al is list of a, it's actually a dictionary.
    #ald is description of physically what each a is.
    def alConstruct(self):
        
        alDiff,aldDiff=self.aDiffConstruct(state=self.state,params=self.params)
        self.multimerIdx=np.where(np.array(map(len,aldDiff['pos']))>1)[0]
        alHydro,aldHydro=self.aHydroConstruct(aldDiff,state=self.state,params=self.params)
        alEnd,aldEnd=self.aEndConstruct(aldDiff,state=self.state,params=self.params)
        al={'diff':alDiff,'hydro':alHydro,'end':alEnd}
        ald={'diff':aldDiff,'hydro':aldHydro,'end':aldEnd}
        return al,ald
    
    
    def drawEvent(self):
        r1,r2=np.random.rand(2)
        alFlat=np.concatenate(self.al.values())
        alCumsum=np.cumsum(alFlat)
        alNormed=alCumsum/alCumsum[-1]
#        event=np.where((alNormed[:-1]<r2) * (alNormed[1:]>r2)==True)[0] <-- not correct when r2 is small than all elements
        event=np.searchsorted(alNormed,r2)
        time=1/alCumsum[-1]*np.log(1/r1)
        
        #-1 at the end as np.searchsorted is <=.
        lens=np.cumsum(map(len,self.al.values()))-1
        
        #create the event dictionary. {type of event: subject index}
        temp=np.searchsorted(lens,event)
        event={self.al.keys()[temp]:event-([0]+list(lens[:-1]))[temp]}
        return event,time
    
    
    def execEvent(self,event,time):
        self.time+=time
        failFlag=False #failFlag indicate whether the event happens or not (could fail due to other constraints)
        if event.keys()[0]=='diff':
            # diffusion event, for oligos or monomers
            leftRight=((event['diff']%2)*2-1)
            pos=self.ald['diff']['pos']
            posIdx=self.getIdx(pos[event['diff']/2])
            newOligoPos=(self.state[0][posIdx]+leftRight)%self.params['N']
            
            self.newState[0][self.newStateIdx[posIdx]]=newOligoPos
            
            
            # update the state, reset newStateIdx
            self.state=[x.copy() for x in self.newState]
            self.newStateIdx=np.arange(self.params['m'])            
            
            
            ##############
            #update al and ald
            
            endCared=-(leftRight+1)/2 #-1 if right, 0 if left.
            tempFunc = lambda x: (newOligoPos[endCared]+leftRight)%self.params['N'] in x
        
#                tempFunc = lambda x: newOligoPos[-1]+1 in x
            exist=np.where(map(tempFunc,self.ald['diff']['pos']))[0]
            if exist.size!=0: # there is new linking, exist should have only 1 element
                eventIdx=event['diff']/2

                #update ald['end'] first, before ald['diff'] is changed
                # 4 cases: multimer + multimer,  multimer + monomer, monomer + multimer, monomer+monomer
                if len(newOligoPos)>1:
                    numMultimerLeft=sum(self.multimerIdx<event['diff']/2)
                    if len(self.ald['diff']['pos'][exist])>1:
                        # multimer + multimer, delete between ends
                        #merged multimer is immediately to the right
                        #endCared is -1 if right, 0 if left. we want 1 if right, 0 if left.
                        idx_delete=(2*numMultimerLeft-endCared,2*numMultimerLeft+2*leftRight+endCared)
                        self.ald['end']=np.delete(self.ald['end'],idx_delete)
                        self.al['end']=np.delete(self.al['end'],idx_delete)
                        #now move the other unmerged end as well, due to diffusion
                        self.ald['end'][idx_delete-leftRight]+=leftRight
                    else:
                        # multimer + monomer
                        # substitute interacting end with new end energy
                        idx_sub=2*numMultimerLeft-endCared
                        #substitute position, using knowledge they must be next to each other
                        #times 2 because we need to diffuse first, then merge.
                        self.ald['end'][idx_sub]=(self.ald['end'][idx_sub]+leftRight*2)%self.params['N']
                        self.ald['end'][idx_sub-leftRight]=(self.ald['end'][idx_sub-leftRight]+leftRight)%self.params['N']
                        # substitute rate                        
                        newAlpha, _ =self.aEndCalc([newOligoPos],self.state,self.params)
                        self.al['end'][idx_sub]=newAlpha[-endCared]
                else:
                    # note this numMultimerLeft is left to Multimer that does not move.
                    numMultimerLeft=sum(self.multimerIdx<exist)
                    if len(self.ald['diff']['pos'][exist])>1:
                        
                        # monomer + multimer
                        # substitute interacting end with new end energy

                        endCared=-endCared #as here the end is wrt non-moving multimer, left gives 1, corresponding to right-end of non-moving multimer, right gives 0.
                        idx_sub=2*numMultimerLeft-endCared
                        #substitute position, using knowledge they must be next to each other
                        self.ald['end'][idx_sub]=(self.ald['end'][idx_sub]+leftRight)%self.params['N']
                        # substitute rate                        
                        newAlpha, _ =self.aEndCalc([newOligoPos],self.state,self.params)
                        self.al['end'][idx_sub]=newAlpha[endCared]
                    
                    else:
                        # monomer + monomer
                        # create a new pair of ends emerging at 2*numMultimerLeft
                        if leftRight==-1:
                            pairPos=np.concatenate((newOligoPos-1,newOligoPos))
                        else:
                            pairPos=np.concatenate((newOligoPos,newOligoPos+1))
                        self.ald['end']=np.insert(self.ald['end'],2*numMultimerLeft,pairPos)
                        
                        newAlpha, _ = self.aEndCalc([pairPos],self.state,self.params)
                        self.al['end']=np.insert(self.al['end'],2*numMultimerLeft,newAlpha)
                    
                    
                #update ald['diff']
                if endCared: #endCared = -1, right.
                    self.ald['diff']['pos'][eventIdx]=np.concatenate((newOligoPos,self.ald['diff']['pos'][exist]))
                else: #left
                    self.ald['diff']['pos'][eventIdx]=np.concatenate((self.ald['diff']['pos'][exist],newOligoPos))
                    
                del(self.ald['diff']['pos'][exist])
                #unnecessary full version:
#                    self.ald['diff']['leftRight']=np.delete(self.ald['diff']['leftRight'],(exist*2,exist*2+1))
                #faster short version:
                self.ald['diff']['leftRight']=self.ald['diff']['leftRight'][:-2]
                
                #update al['diff']
                self.al['diff']=np.delete(self.al['diff'],(exist*2,exist*2+1))
                newDiff=self.aDiffCalc([self.ald['diff']['pos'][eventIdx]])
                self.al['diff'][eventIdx]=newDiff
                self.al['diff'][eventIdx+1]=newDiff
                 
                
                #alternatively, update ald['end'] later, reconstruct                
#                self.al['end'],self.ald['end']=self.aEndConstruct(self.ald['diff'],self.state,self.params)
                
        
        elif event.keys()[0]=='hydro':
            # hydrolysis event, change the selected subunit/unit into ADP
            monomer=np.array(self.state[0][self.ald['hydro'][event['hydro']]])
            tempFunc=lambda x:monomer in self.ald['diff']['pos'][x]
            
            if any(map(tempFunc,np.where(np.array(map(len,self.ald['diff']['pos']))==1)[0])):
                failFlag=True
            else: 
                self.state[2][event['hydro']]=1
                self.newState[2][event['hydro']]=1
        
            #update al and ald for hydro
            np.delete(self.ald['hydro'],event['hydro'])
            np.delete(self.al['hydro'],event['hydro'])
        
        
        elif event.keys()[0]=='end':
            
            # end-breaking event
            endIdx=event['end']
            leftRight=((endIdx%2)*2-1)
            endCared=(leftRight+1)/2 # 1 if right, 0 if left
            pos=self.ald['end'][endIdx]
            idx=self.getIdx([pos])
            self.newState[0][idx]=(self.newState[0][idx]+leftRight)%self.params['N']
            
            #update al and ald, this involves end, diff, and hydro
            
            #update end, six cases; 2-mer or multimer --> monomer, monomer + monomer, monomer+multimer
            # for diffusion. six cases
            oligoNonMoving_idx=np.where([pos in x for x in self.ald['diff']['pos']])[0]
            oligoNonMoving=self.ald['diff']['pos'][oligoNonMoving_idx]
            
            # mergingEnd_idx is used to differentiate whether something is delete or not
            if len(oligoNonMoving)>2:
                # oligoNonMoving remains multimer.
                # substitute breaking end

                # update newly exposed end position
                self.ald['end'][endIdx]=(self.ald['end'][endIdx]-leftRight)%self.params['N'] 
                if endCared:
                    newOligoNonMoving=oligoNonMoving[:-1]
                else: newOligoNonMoving=oligoNonMoving[1:]
                newAlpha, _ = self.aEndCalc([newOligoNonMoving],self.state,self.params)
                self.al['end'][endIdx]=newAlpha[endCared] # update newly exposed end energy
                movingEnd_idx=endIdx

                    
            else:
                # oligoNonMoving becomes a monomer 
                # delete nonmoving end
                nonMovingEndIdx=(endIdx-leftRight)%self.params['m']
                self.ald['end']=np.delete(self.ald['end'],nonMovingEndIdx) #delete non-moving new monomer
                self.al['end']=np.delete(self.al['end'],nonMovingEndIdx)
                
                # if going right, left one is deleted, so the next right one to interact is +0 for mergingEnd_idx, -1 for movingEnd_idx;
                # if going left, right one is deleted does not affect movingEnd_idx, mergingEnd_idx is the more left one, -1.
                movingEnd_idx=(endIdx-endCared)%self.params['m'] # right, -1; left, 0
            
            
            mergingPos=(pos+leftRight*2)%self.params['N']
            oligoMerging_idx=np.where([mergingPos in x for x in self.ald['diff']['pos']])[0]
            
            if oligoMerging_idx.size==0:
                # going into empty space
            
                # multimer case, do nothing.
                if len(oligoNonMoving)<3:
                # monomer case, delete movingEnd
                    pdb.set_trace()
                    self.ald['end']=np.delete(self.ald['end'],movingEnd_idx)
                    self.al['end']=np.delete(self.al['end'],movingEnd_idx)
                
                
            else:
                oligoMerging=self.ald['diff']['pos'][oligoMerging_idx]
                if len(oligoMerging)==1:
                    # oligoMerging is monomer
                    # substitute movingEnd position and new energy
                    # insert the mergedEnd.
                
                    # in this case, mergingEnd_idx = movingEnd_idx +0 if left, +1 if right.
                    mergingEnd_idx=(movingEnd_idx+endCared)%self.params['N']
                    # update movingEnd's ald['end'] position
                    self.ald['end'][movingEnd_idx]=(mergingPos-leftRight)%self.params['N']
                    # insert new end position emerged from merged two monomers forming 2-mer.
                    self.ald['end']=np.insert(self.ald['end'],mergingEnd_idx,mergingPos)
                    
                    # find newly formed 2-mer's position
                    if endCared: #right
                        newOligoMerging=np.array([mergingPos-1,mergingPos])%self.params['N']
                    else: newOligoMerging=np.array([mergingPos,mergingPos+1])%self.params['N']
                    
                    # calculate new alpha
                    newAlpha,_=self.aEndCalc([newOligoMerging],self.state,self.params)
                    
                    #they should be the same. same pair.
                    self.al['end'][movingEnd_idx]=newAlpha[1-endCared]
                    np.insert(self.al['end'],mergingEnd_idx,newAlpha[endCared])
                else:
                    # oligoMerging is multimer
                    # delete movingEnd, substitute mergingEnd
                    self.ald['end']=np.delete(self.ald['end'],movingEnd_idx)
                    self.al['end']=np.delete(self.al['end'],movingEnd_idx)
                    #delete first means mergingEnd_idx should be movingEnd_idx-1 if left, +0 if right.
                    mergingEnd_idx = (movingEnd_idx+(endCared-1))%self.params['N']
                    self.ald['end'][mergingEnd_idx]=(self.ald['end'][mergingEnd_idx]-leftRight)%self.params['N']
                    
                    if endCared: # right
                        newOligoMerging=np.concatenate(([pos+1],oligoMerging))
                    else: newOligoMerging=np.concatenate((oligoMerging,[pos-1]))
                    newAlpha,_=self.aEndCalc([newOligoMerging],self.state,self.params)
                    self.al['end'][mergingEnd_idx]=newAlpha[1-endCared]
                    
            pdb.set_trace()
#            if self.state[0]            
            map(tempFunc,np.where(np.array(map(len,self.ald['diff']['pos']))==1)[0])
            aa=1
            
        
        
        # update state and newState. in this case we just need one operation
        self.state[0][idx]=(self.state[0][idx]+leftRight)%self.params['N']
        return event,failFlag
    
    
    def alUpdate(self,event=None,failFlag=None):
        if event==None:
            event=self.event
        
        if event.keys()[0]=='diff':
            alDiff,aldDiff=self.aDiffConstruct(state=self.state,params=self.params)
            alEnd,aldEnd=self.aEndConstruct(aldDiff,state=self.state,params=self.params)
            self.al['diff']=alDiff
            self.ald['diff']=aldDiff
            self.al['end']=alEnd
            self.ald['end']=aldEnd
        
        elif event.keys()[0]=='hydro' and not failFlag:
            np.delete(self.al['hydro'],event['hydro'])
        
        elif event.keys()[0]=='end':
            pass
        return al
    
    def step(self):
        event,time=self.drawEvent()
        event,failFlag=self.execEvent(event,time)
        self.alUpdate(event,failFlag)
    
    
    def stateSort(self):
        sortIdx=np.argsort(self.state[0])
        tempFunc=lambda x: x[sortIdx]
        self.state=map(tempFunc,self.state)
        return sortIdx
    

    def aDiffCalc(self,segments):
        #use sqrt scaling law
        alDiff=1/np.sqrt(map(len,segments))*self.params['D_mono']*self.params['T']
        return alDiff    


    def consecutive(self,state,N):
        pos=state[0]
        endI=np.where(np.diff(pos)!=1)[0]
        segments=np.split(pos, endI+1)
        
        if (segments[0][0]-segments[-1][-1])%N==1:
            #i.e. first and last segment is one segment
        
            #before join, first modify state[1] and state[2] to reflect change in position, state[0]
            start=self.getIdx([segments[-1][0]])[0]
            
            permFunc=lambda x:np.concatenate((x[start:],x[:start]))
            newState=map(permFunc,state)
            self.newState=newState
            reverseStart=self.params['m']-start
            self.newStateIdx=np.concatenate((self.newStateIdx[reverseStart:],self.newStateIdx[:reverseStart]))
            #did not modify state[0] for future getIdx usage, such as in constructing aldEnds
            
            joined_ends=np.concatenate((segments[-1],segments[0])) 
            #[0,1,2] and [N-2,N-1] becomes [B-2,N-1,0,1,2]
            segments[0]=joined_ends
            del(segments[-1])
        else: 
            endIdx=np.append(endI,len(pos)-1)
        return segments
    
    def aDiffConstruct(self,state=None,params=None):
        #construct diffusion alpha list.
        
        self.stateSort()
        segments=self.consecutive(self.state,self.params['N'])
        aldDiff={'pos':segments}
        aldDiff['leftRight']=np.tile([0,1],len(segments))
        alDiff=np.tile(self.aDiffCalc(segments),(2,1)).T.flatten()
        return alDiff,aldDiff


    def aHydroConstruct(self,aldDiff,state=None,params=None):
        # construct hydrolysis alpha list. 
        # return list of indices that are ATP.
        if params==None:
            params=self.params
        if state==None:
            state=self.state
        
        k_h=params['rateHydro']
        pdb.set_trace()
        posInMultimer=np.concatenate(np.array(aldDiff['pos'])[self.multimerIdx])
        idxATP=np.where(state[2]==0)[0] #0 is ATP
        
        alHydro=np.ones(len(aldHydro))*k_h
        return alHydro, aldHydro
        
        
        
    
    def aEndConstruct(self,aldDiff,state=None,params=None):    
        # construct end-breaking (or, disassembly) a list.
                
        if state==None:
            state=self.State
        if params==None:
            params=self.params
    
        # only consider oligos of length more than 1
        aldDiff_unique=aldDiff['pos']
        oligos=[aldDiff_unique[x] for x in self.multimerIdx]
        
        alEnd,aldEnd=self.aEndCalc(oligos,state,params)
        return alEnd,aldEnd
        
    def aEndCalc(self,oligos,state,params):
    #take the ends and the monomer immediately adjacent to it to calculate change in Energy
        endFun=lambda seg: [seg[:2],seg[-2:]]
        ends=np.reshape(map(endFun,oligos),(len(oligos)*2,2))
        endLastFun=lambda seg:[seg[0],seg[-1]]
        aldEnd=np.array(map(endLastFun,oligos)).flatten()
        getDomains=lambda idx: [state[1][idx[0]][1],state[1][idx[1]][0]]
        getForm=lambda idx:state[2][idx]
        endIdx=np.reshape(self.getIdx(ends.flatten()),ends.shape)
        
        endDomains=np.reshape(map(getDomains,endIdx),(len(oligos)*2,2))
        endForms=np.reshape(map(getForm,endIdx.flatten()),(len(oligos)*2,2))
        
        delE=self.H(endDomains,endForms)
        
        #delE is the adjacent units' current interaction energy
        # taking them apart would result in -delE energy change
        # since probability calc is using -(change in energy), 
        # it is -(-delE) = delE
        alEnd=np.exp(delE/params['T'])
        
        return alEnd,aldEnd
        
    
#    def oligoStep(self):
#        
#        
#            
#        def filament(state,params,minLen):
#            segments,sortIdx,endIdx=consecutive(state,params)
#            mask=np.array(map(len,segments))>=minLen
#            fila=np.array(segments)[mask]
#            return fila,sortIdx,endIdx,mask  
#        N=self.params['N']
#        m=self.params['m']
#        oligo,sortIdx,endIdx,mask=filament(self.state,self.params,2)
#        if list(oligo):
##            print(self.state[0])
#            moveBool=np.random.rand(len(oligo))<1./np.array(map(len,oligo))
#            for j in xrange(len(moveBool)):
#                if moveBool[j]:
##                    if len(oligo)>1: pdb.set_trace()
#                    if np.random.rand()<0.5 and (oligo[j][0]-1)%N not in self.state[0]:
#                        # move to the left.
#                    
##                        if len(oligo)>1:pdb.set_trace()
#                        changeIdx=np.roll(sortIdx,m-endIdx[j]-1)[-len(oligo[j]):]
#                        self.state[0][changeIdx]=(self.state[0][changeIdx]-1)%N
#                        if self.checkOverlap(): pdb.set_trace()
#                    elif (oligo[j][-1]+1)%N not in self.state[0]:
##                        if len(oligo)>1:pdb.set_trace()
#                        changeIdx=np.roll(sortIdx,m-endIdx[j]-1)[-len(oligo[j]):]
#                        self.state[0][changeIdx]=(self.state[0][changeIdx]-1)%N
#                        if self.checkOverlap(): pdb.set_trace()
#
#
#
#
# 
#    def energy(self,state):
#        #here state could be the current state or some hypothetical state
#        idx=np.argsort(state[0])
#        state_sorted=(state[0][idx],state[1][idx],state[2][idx])
#        E=0
#        d=int(self.params['dist_cutoff'])
#        for i in xrange(self.params['m']):
#            for j in xrange(d):
#                rIdx=(i+j+1)%self.params['m']
#                a=[state_sorted[0][i],state_sorted[2][i],state_sorted[1][i]]
#                b=[state_sorted[0][rIdx],state_sorted[2][rIdx],state_sorted[1][rIdx]]
#                if (a[0]-b[0])%self.params['N']==0:pdb.set_trace()
#                E+=self.H(a,b)
#        return E
#        
#   
#    def checkOverlap(self):
#        elt,repeats=np.unique(self.state[0],return_counts=True)
#        return (repeats>1).any()
#    
#
#
#    def TtoD(self):
#        #p_h, hydroplysis probability for each ATP-form protein
#        p_h=self.params['rateHydro']
#        posIdx=np.array(range(self.params['m']))[self.state[2]==0]
#        for idx in posIdx:
#            if np.random.rand()<=p_h:
#                self.state[2][idx]=1
    
    def checkATP(self,state):
        #turn all ADP form monomers into ATP
        # 0 is ATP state. 1 is ADP state
        hydro=state[2].copy()
        for t in xrange(self.params['m']):
            if state[2][t]!=0:
                if (state[0][t]+1)%self.params['N'] not in state[0] and (state[0][t]-1)%self.params['N'] not in state[0]:
                    hydro[t]=0
        return (state[0],state[1],hydro)


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
    for ind, res in enumerate(pool.imap(parFun, iterList),chunksize):
        energies[ind-1],states[ind-1],stepSizeVec[ind-1]=res
#    parFun(iterList[0])
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
    H=lambda endDomainVec,endFormVec: H_actin(endDomainVec,endFormVec,eps,xi)
    charTime=max(min(np.exp(-xi/T),1e3),1)
    stepSize=np.floor(charTime)
    for repeat in xrange(simPerPt):
        intSp=InteractionSpace(params,H)
        for t in xrange(obsStart): 
            intSp.overallStep()
        for t in xrange(nPts):
            energy[repeat][t]=intSp.currentE
            state[repeat].append(intSp.state)
            for step_idx in xrange(int(stepSize)):
                intSp.overallStep()
    return energy,state,stepSize

def H_actin(endDomainVec,endFormVec,eps,xi):
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

    getE=lambda idx:H_array[endFormVec[idx][0]][endFormVec[idx][1]][endDomainVec[idx][0]][endDomainVec[idx][1]]
    return np.array(map(getE,xrange(len(endDomainVec))))
        

    


def main(fName,start,nPts,hydroFlag=True):
    
    N=30
    k=2
    h=2
    m_array=np.array([[0,20],[0,0]])
    
    
    if hydroFlag==False:
        rateHydro=0.
    else: rateHydro=1e-3

    # T, xi, eps, to be filled in.
    initParams={'N':N,'m':np.sum(m_array),'k':k,'h':h,
                'T':None,'xi':None,'eps':None,
                'isCircular':True,
                'transRotRatio':1.,
                'rateHydro':rateHydro,
                'm_array':m_array,
                'dist_cutoff':2.,
                'D_mono':1.}
                # D_mono/T should be on the scale of 1e-8, but we change the unit to make it 1 here.
                # T is similar. Instead some number between 200 and 300, we set it to 1.

    #write to animation
#    TE_animate('1dActin',initParams,50,1000)

                
    T_array=np.array([1.])
#    xi_array=-1*np.linspace(1,8,num=4)
#    eps_array=-1*np.linspace(2,10,num=4)
    xi_array=np.array([-4.])
    eps_array=np.array([-8.])
    TE_simulation(fName,initParams,T_array,xi_array,eps_array,simPerPt=1,obsStart=start,nPts=nPts)



if __name__ == "__main__":
#    cProfile.run('main()','profile.tmp')
#    python OneD_actin.py asdf.p 150 100

#    fName = sys.argv[1]
#    start = int(sys.argv[2])
#    nPts = int(sys.argv[3])
#    if len(sys.argv)>4:
#        hydroFlag=int(sys.argv[4])
#        main(fName,start,nPts,hydroFlag)
#        
#    else: main(fName,start,nPts)
    N=30
    k=2
    h=2
    m_array=np.array([[0,20],[0,0]])
    
    hydroFlag=True
    if hydroFlag==False:
        rateHydro=0.
    else: rateHydro=1e-3
    initParams={'N':N,'m':np.sum(m_array),'k':k,'h':h,
                'T':1.,'xi':None,'eps':None,
                'isCircular':True,
                'transRotRatio':1.,
                'rateHydro':rateHydro,
                'm_array':m_array,
                'dist_cutoff':2.,
                'D_mono':1.}
    aa=lambda x,y:H_actin(x,y,-1.,-10.)
    aa=InteractionSpace(initParams,aa)
#    pdb.set_trace()
#    aa.step()
    aa.drawEvent()
    
#    aa.execEvent({'end':1},2.)
    aa.execEvent({'end':2},2.) 
