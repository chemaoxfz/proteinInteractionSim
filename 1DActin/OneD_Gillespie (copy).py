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
import pandas as pd

class InteractionSpace:
    def __init__(self,initParams,H):
#        initParams['m_array']=np.array(initParams['m_array'])
#        position=np.array([0,1,2,3,4,5,6,7,8,9,12,13,15,16,20,22,25,27,28,29])
#{'diff': 0}
#[ 1  2  9 10 12 13 14 15 17 19 20 21 22 23 24 25 26 27 28 29]      
#        position=np.array([ 1 , 2,  9,  10 , 12 ,13 ,14 ,15, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
#    {'end': 0}
#        position=np.array([ 0 , 1 , 9 ,11 ,12 ,13 ,14, 15, 17 ,19, 20, 21, 22 ,23 ,24 ,25 ,26 ,27 ,28 ,29])
        
        #################### END BREAKING TEST CASES ##########################3
        ### movingPos=0
        ## LEFT, {'end':0}
        # nonmoving: 2-mer; merging: empty
#        position=np.array([ 1,  2,4, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26, 27])
        # nonmoving: 2-mer; merging: monomer
#        position=np.array([ 1,  2,4, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26, 29])
        # nonmoving: 2-mer; merging: multimer
#        position=np.array([ 1,  2,4, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,28, 29])
        # nonmoving: multimer; merging: empty
#        position=np.array([ 1,  2,3, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26, 27])
        # nonmoving: multimer; merging: monomer
#        position=np.array([ 1,  2,3, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26, 29])
        # nonmoving: multimer; merging: multimer        
#        position=np.array([ 1,  2,3, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,28, 29])
        
        ## Right, {'end':5}
        # nonmoving: 2-mer; merging: empty
#        position=np.array([ 2,  3,10, 11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,28, 29])
        # nonmoving: 2-mer; merging: monomer
#        position=np.array([ 1,  3,4, 11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,28, 29])
        # nonmoving: 2-mer; merging: multimer
#        position=np.array([ 1,  2,3, 11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,28, 29])
        # nonmoving: multimer; merging: empty
#        position=np.array([ 2,  3,4, 11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 27,28, 29])
        # nonmoving: multimer; merging: monomer
#        position=np.array([ 1,  3,4, 11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 27,28, 29])
        # nonmoving: multimer; merging: multimer        
#        position=np.array([ 1,  2,3, 11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 27,28, 29])        
        
        
        
        ### mergingPos=0
        ## LEFT, {'end':0}
        # nonmoving: 2-mer; merging: empty
#        position=np.array([ 2,  3,9, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26, 27])
        # nonmoving: 2-mer; merging: monomer
#        position=np.array([ 0,  2,  3,9, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26])
        # nonmoving: 2-mer; merging: multimer
#        position=np.array([0, 2,  3,9, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24,28, 29])
        # nonmoving: multimer; merging: empty
#        position=np.array([ 2,  3,4, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26, 27])
        # nonmoving: multimer; merging: monomer
#        position=np.array([0, 2,  3,4, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,26])
        # nonmoving: multimer; merging: multimer        
#        position=np.array([0, 2,  3,4, 10, 11,12, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24,28, 29])
        ## Right, {'end':5}
        # nonmoving: 2-mer; merging: empty
#        position=np.array([ 2,  3,10, 11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,27, 28])
        # nonmoving: 2-mer; merging: monomer
#        position=np.array([ 0,  3,4,   11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,   27, 28])
        # nonmoving: 2-mer; merging: multimer
#        position=np.array([ 0,1,2,   11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24, 25,   27, 28])
        # nonmoving: multimer; merging: empty
#        position=np.array([ 2,  3,4,    11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24,    26,27, 28])
        # nonmoving: multimer; merging: monomer
#        position=np.array([ 0,  3,4,    11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24,    26,27, 28])
        # nonmoving: multimer; merging: multimer        
#        position=np.array([ 0,1,2,    11, 12,13, 14,15, 16, 17 ,18, 19, 20, 21, 22, 23, 24,      26,27, 28]) 
        
        
        ################ DIFFUSION TEST CASES ####################
        
        position=np.array([ 0,  8, 10, 11, 13, 15, 16 ,17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
#        position_1hot=np.zeros(initParams['N'])
#        position_1hot[position]=1

#        position=np.array(range(initParams['m']))
#        position=np.random.choice(range(initParams['N']),initParams['m'],replace=False)
        k=initParams['k']
        #m[0,1] entry gives number of proteins that starts with [1,0] domains.
        proteinOrientation=np.repeat(np.reshape(np.ravel(np.indices((k,k)),1),[k**2,2]),initParams['m_array'].flatten(),axis=0)
#        proteinOrientation_1hot=-np.ones(initParams['N'])
        ################################## 1 HOT stuff not finished. #############################
        
#        proteinForm=np.floor(np.random.rand(initParams['m'])*initParams['k'])
        proteinForm=np.floor(np.random.rand(initParams['N'])*initParams['k'])
        mask=np.ones(initParams['N'])
        mask[position]=0
        proteinForm[np.where(mask)]=-1 #-1 indicate no protein, 0 is ATP, 1 is ADP
        self.proteinForm_1hot=proteinForm
        proteinForm_pos=proteinForm[position]

        self.init_state=[position,proteinOrientation,proteinForm_pos]
        self.params=initParams
        self.state,self.proteinForm_1hot=self.checkATP(self.init_state,self.proteinForm_1hot)
        self.newState=[x.copy() for x in self.state]
        self.stateSort()
        self.H=H
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
        event=np.searchsorted(alNormed,r2)
        time=1/alCumsum[-1]*np.log(1/r1)
        
        #-1 at the end as np.searchsorted is <=.
        lens=np.cumsum(map(len,self.al.values()))-1e-7
        #create the event dictionary. {type of event: subject index}
        temp=np.searchsorted(lens,event)
        eventDic={self.al.keys()[temp]:event-int(round(([0]+list(lens[:-1]))[temp]))}
#        pdb.set_trace()
        return eventDic,time
    
    
    def execEvent(self,event,time):
        self.time+=time
        if event.keys()[0]=='diff':
            # diffusion event, for oligos or monomers
            leftRight=((event['diff']%2)*2-1)
            pos=self.ald['diff']['pos']
            oldOligoPos=pos[event['diff']/2]
            posIdx=self.getIdx(oldOligoPos)
            newOligoPos=(self.state[0][posIdx]+leftRight)%self.params['N']
            self.newState[0][posIdx]=newOligoPos
            
            ##############
            #update al and ald
            
            endCared=-(leftRight+1)/2 #-1 if right, 0 if left.
            oppEnd=-1-endCared # the other end, 0 if right, -1 if left
            tempFunc = lambda x: (newOligoPos[endCared]+leftRight)%self.params['N'] in x
        
            exist=np.where(map(tempFunc,self.ald['diff']['pos']))[0]
            
            eventIdx=event['diff']/2
            if exist.size!=0: # there is new linking, exist should have only 1 element
                mergingOligoPos=self.ald['diff']['pos'][exist]
                
                if endCared:
                    # endCared=-1, i.e. right.
                    endCalcOligo=np.concatenate((newOligoPos,mergingOligoPos))
                    endCalcOligoOldPos=np.concatenate((oldOligoPos,mergingOligoPos))
                else:
                    endCalcOligo=np.concatenate((mergingOligoPos,newOligoPos))
                    endCalcOligoOldPos=np.concatenate((mergingOligoPos,oldOligoPos))
                #update ald['end'] first, before ald['diff'] is changed
                # 4 cases: multimer + multimer,  multimer + monomer, monomer + multimer, monomer+monomer
                if len(newOligoPos)>1:
                    numMultimerLeft=sum(self.multimerIdx<event['diff']/2)
                    if len(mergingOligoPos)>1:
                        # multimer + multimer
                        ###### END #####
                        # delete between ends
                        # merged multimer is right next to moving multimer
                        
                        #first move the unmerged end of moving multimer
                        unmergedMovingEndIdx=2*numMultimerLeft+endCared+1
                        
                        self.ald['end'][unmergedMovingEndIdx]=(self.ald['end'][unmergedMovingEndIdx]+leftRight)%self.params['N']
                        
                        #endCared is -1 if right, 0 if left. we want 1 if right, 0 if left.
                        # for the idx_delete for nonmoving multimer, we want +1 if left, 0 if right.
                        numEnds=len(self.ald['end'])
                        mergedEndIdx=(unmergedMovingEndIdx+2*leftRight)%numEnds
                        idx_delete=(2*numMultimerLeft-endCared,mergedEndIdx)
                        self.ald['end']=np.delete(self.ald['end'],idx_delete)
                        self.al['end']=np.delete(self.al['end'],idx_delete)
                        
                        if idx_delete[0]%(numEnds-1)==idx_delete[1]%(numEnds-1):
                            # this is when two oligos at two ends of the space are connected, bring the last one to front
                            self.ald['end']=np.roll(self.ald['end'],1)
                            
                        ###### HYDRO ######
                        # shift the positions for ['hydro']['1hot'] and then reconstruct ['hydro']['1hot']
                        self.ald['hydro']['1hot'][newOligoPos[endCared]]=1
                        self.ald['hydro']['1hot'][oldOligoPos[oppEnd]]=0
                        self.ald['hydro']['pos']=np.where(self.ald['hydro']['1hot'])[0]

                    else:
                        # multimer + monomer
                        ###### END ######
                        # substitute interacting end with new end energy
                        idx_sub=2*numMultimerLeft-endCared
                        #substitute position, using knowledge they must be next to each other
                        #times 2 because we need to diffuse first, then merge.
                        mergingPos=(self.ald['end'][idx_sub]+leftRight*2)%self.params['N']
                        self.ald['end'][idx_sub]=mergingPos
                        self.ald['end'][idx_sub-leftRight]=(self.ald['end'][idx_sub-leftRight]+leftRight)%self.params['N']
                        
                        # substitute rate                        
                        newAlpha, _ =self.aEndCalc([endCalcOligoOldPos],self.state,self.params)
                        self.al['end'][idx_sub]=newAlpha[-endCared]
                        
                        ###### HYDRO ######
                        # add the new monomer in. It's definitely ATP
                        self.ald['hydro']['1hot'][newOligoPos[endCared]]=1
                        self.ald['hydro']['1hot'][oldOligoPos[oppEnd]]=0
                        self.ald['hydro']['1hot'][mergingOligoPos]=1
                        self.ald['hydro']['pos']=np.where(self.ald['hydro']['1hot'])[0]
                        self.al['hydro']=np.append(self.al['hydro'],self.params['rateHydro'])

                else:
                    # note this numMultimerLeft is left to Multimer that does not move.
                    numMultimerLeft=sum(self.multimerIdx<exist)
                    if len(mergingOligoPos)>1:
                        
                        # monomer + multimer
                        ###### END ######
                        # substitute interacting end with new end energy
                        
                        # change endCared to be wrt the non-moving multimer, so right gives 0, left gives -1
                        endCaredTemp=-(endCared+1)
                        # if left, here should be +1, if right, here should be 0. So -endCaredTemp
                        idx_sub=2*numMultimerLeft-endCaredTemp
                        #substitute position, using knowledge they must be next to each other
                        self.ald['end'][idx_sub]=(self.ald['end'][idx_sub]-leftRight)%self.params['N']
                        # substitute rate                        
                        newAlpha, _ =self.aEndCalc([endCalcOligoOldPos],self.state,self.params)
                        self.al['end'][idx_sub]=newAlpha[endCaredTemp]
                        
                        ###### HYDRO ######
                        # add the new monomer in. It's definitely ATP
                        # for position, nothing moves except adding the new one.
                        self.ald['hydro']['1hot'][mergingOligoPos]=1
                        self.ald['hydro']['pos']=np.where(self.ald['hydro']['1hot'])[0]
                        self.al['hydro']=np.append(self.al['hydro'],self.params['rateHydro'])
                        
                    
                    else:
                        # monomer + monomer
                        ###### END ######
                        # create a new pair of ends emerging at 2*numMultimerLeft
                        self.ald['end']=np.insert(self.ald['end'],2*numMultimerLeft,endCalcOligo)
                        
                        newAlpha, _ = self.aEndCalc([endCalcOligoOldPos],self.state,self.params)
                        self.al['end']=np.insert(self.al['end'],2*numMultimerLeft,newAlpha)
                        
                        ###### HYDRO ######
                        # Add both merging and moving oligos. 
                        # for position, nothing moves except adding the new one.
                        self.ald['hydro']['1hot'][mergingOligoPos]=1
                        self.ald['hydro']['1hot'][newOligoPos]=1
                        self.ald['hydro']['pos']=np.where(self.ald['hydro']['1hot'])[0]
                        self.al['hydro']=np.append(self.al['hydro'],[self.params['rateHydro'],self.params['rateHydro']])
                        

                ##### DIFFUSION UPDATE #####
                
                newDiff=self.aDiffCalc([endCalcOligo])
                self.al['diff'][eventIdx*2]=newDiff
                self.al['diff'][eventIdx*2+1]=newDiff
                self.al['diff']=np.delete(self.al['diff'],(exist*2,exist*2+1))
                # ald                
                self.ald['diff']['pos'][eventIdx]=endCalcOligo
                del(self.ald['diff']['pos'][exist])
                
            else:
                # no new linking, only diffusion.
                self.ald['diff']['pos'][eventIdx]=newOligoPos
                ###### END ######
                # also update 'end'
                numMultimerLeft=sum(self.multimerIdx<event['diff']/2)

                if len(newOligoPos)>1:
                    self.ald['end'][numMultimerLeft*2]=newOligoPos[0]
                    self.ald['end'][numMultimerLeft*2+1]=newOligoPos[-1]
                
                ###### HYDRO ######
                # shift position
                self.ald['hydro']['1hot'][newOligoPos[endCared]]=1
                self.ald['hydro']['1hot'][oldOligoPos[oppEnd]]=0
                self.ald['hydro']['pos']=np.where(self.ald['hydro']['1hot'])[0]
            
            # ROLLING ################
            # going to the left, while last subunit is 0, this case we need to take the new oligo to the last index
#            pdb.set_trace()
            
            if exist.size!=0:
                # new linking
                movingEndPos=endCalcOligo[endCared]
                mergingEndPos=mergingOligoPos[-1-endCared]
                if len(newOligoPos)>1:
                    if len(mergingOligoPos)>1:
                        if endCared and mergingEndPos<2:
                            self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],1)
                            self.al['diff']=np.roll(self.al['diff'],2)
                            self.ald['end']=np.roll(self.ald['end'],1)
                            self.al['end']=np.roll(self.al['end'],1)
                        elif not endCared and self.params['N']-merging<3:
                            self.ald['end']=np.roll(self.ald['end'],1)
                            self.al['end']=np.roll(self.al['end'],1)
                    else:
                        if endCared and mergingEndPos<2: #right
                            self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],1)
                            self.al['diff']=np.roll(self.al['diff'],2)
                            self.ald['end']=np.roll(self.ald['end'],2)
                            self.al['end']=np.roll(self.al['end'],2)
                else: #movingOligo is monomer
                    if endCared:
                        #right
                        if mergingEndPos<2:
                            self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],1)
                            self.al['diff']=np.roll(self.al['diff'],2)
                    else:
                        if mergingEndPos==self.params['N']-1:
                            self.ald['end']=np.roll(self.ald['end'],2)
                            self.al['end']=np.roll(self.al['end'],2)
                        elif movingEndPos==self.params['N']-1:
                            self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],-1)
                            self.al['diff']=np.roll(self.al['diff'],-2)
            else:
                # no linking
                if endCared and pos[endCared]==self.params['N']-1: #right and diffuse to 0
                    self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],1)
                    self.al['diff']=np.roll(self.al['diff'],2)
                    self.ald['end']=np.roll(self.ald['end'],2)
                    self.al['end']=np.roll(self.al['end'],2)
                    
                
            
            ###### STATE ######
            # update the state from ald and al's
#            self.state=[x.copy() for x in self.newState]
                
            self.state[0]=self.newState[0].copy()
#            self.newState[0]=self.newState[0][self.newStateIdx]
#            self.newStateIdx=np.arange(self.params['m'])
        
        elif event.keys()[0]=='hydro':
            # hydrolysis event, change the selected subunit/unit into ADP
            pos=self.ald['hydro']['pos'][event['hydro']]
            self.ald['hydro']['1hot'][pos]=1
            self.ald['hydro']['pos']=np.delete(self.ald['hydro']['pos'],event['hydro'])
            self.al['hydro']=np.delete(self.al['hydro'],event['hydro'])
            # update state and newState
            idx=self.getIdx(pos)
            self.state[2][idx]=1
            self.newState[2][idx]=1

        
        
        elif event.keys()[0]=='end':
            
            # end-breaking event
            endIdx=event['end']
            leftRight=((endIdx%2)*2-1)
            endCared=(leftRight+1)/2 # 1 if right, 0 if left
            pos=self.ald['end'][endIdx]
            movingPos=(pos+leftRight)%self.params['N']
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
                else: 
                    newOligoNonMoving=oligoNonMoving[1:]
                newAlpha, _ = self.aEndCalc([newOligoNonMoving],self.state,self.params)
                self.al['end'][endIdx]=newAlpha[endCared] # update newly exposed end energy
                movingEnd_idx=(endIdx+endCared)%len(self.ald['end'])

            else:
                # oligoNonMoving becomes a monomer 
                # delete nonmoving end AND moving end
                # if moving end gets attached later, substitute multimer merged end, or create merged 2-mer's 2 ends.
                nonMovingEndIdx=(endIdx-leftRight)%len(self.ald['end'])
                self.ald['end']=np.delete(self.ald['end'],nonMovingEndIdx) #delete non-moving new monomer
                self.al['end']=np.delete(self.al['end'],nonMovingEndIdx)
                
                # if going right, left one is deleted, so -1 for movingEnd_idx;
                # if going left, right one is deleted does not affect movingEnd_idx
                movingEnd_idx=(endIdx-endCared)%len(self.ald['end']) # right, -1; left, 0
                # update state[2] for new monomer from non-moving end
                self.newState[2][self.getIdx(oligoNonMoving)]=0
                self.ald['end']=np.delete(self.ald['end'],movingEnd_idx)
                self.al['end']=np.delete(self.al['end'],movingEnd_idx)
            
            
            mergingPos=(pos+leftRight*2)%self.params['N']
            oligoMerging_idx=np.where([mergingPos in x for x in self.ald['diff']['pos']])[0]

            if oligoMerging_idx.size==0:
                # going into empty space
                
                ##### HYDRO UPDATE #####
                if len(oligoNonMoving)<3:
                    # 2-mer case. delete both nonmoving and moving
                    self.ald['hydro']['1hot'][pos]=0
                    self.ald['hydro']['1hot'][newOligoNonMoving]=0
                    self.ald['hydro']['pos']=np.where(self.ald['hydro']['1hot'])[0]
                    self.al['hydro']=np.delete(self.al['hydro'],(0,1))
                    
#                    mask=self.ald['hydro']==pos
#                    if any(mask):
#                        self.ald['hydro']=np.delete(self.ald['hydro'],np.where(mask)[0])
#                    if endCared:
#                        newOligoNonMoving=oligoNonMoving[:-1]
#                    else:
#                        newOligoNonMoving=oligoNonMoving[1:]
#                    mask=self.ald['hydro']==newOligoNonMoving
#                    if any(mask):
#                        self.ald['hydro']=np.delete(self.ald['hydro'],np.where(mask)[0])
                        
                    # update state[2]
                    self.newState[2][idx]=0
                    self.newState[2][(idx+leftRight)%self.params['m']]=0
                    
                else:
                    # multimer case. delete moving end
                    self.ald['hydro']['1hot'][pos]=0
                    self.ald['hydro']['pos']=np.where(self.ald['hydro']['1hot'])[0]
                    self.al['hydro']=np.delete(self.al['hydro'],0)
#                    mask=self.ald['hydro']==pos
#                    if mask.any():
#                        self.ald['hydro']=np.delete(self.ald['hydro'],np.where(mask)[0])


                ##### END UPDATE #####
                # do nothing
                
                ##### DIFFUSION UPDATE #####
                # empty, substitute non-moving, insert moving
                self.ald['diff']['pos'][oligoNonMoving_idx]=np.delete(self.ald['diff']['pos'][oligoNonMoving_idx],-endCared)
                # insert pos +1 if to the right, +0 if to the left. which is +endCared
                insPosIdx=oligoNonMoving_idx+endCared
                self.ald['diff']['pos'].insert(insPosIdx,np.array([(mergingPos-leftRight)%self.params['N']]))
                
            else:
                # merging with some oligo
                oligoMerging=self.ald['diff']['pos'][oligoMerging_idx]
                
                ##### END UPDATE #####
                if len(oligoMerging)==1:
                    # oligoMerging is monomer
                
                    # insert two mergedEnds and their energy
                    # in this case, mergingEnd = movingEnd_idx.
                    mergingEnd_idx=movingEnd_idx
#                    if mergingEnd_idx==0 or mergingEnd_idx==len(self.ald['end']):
#                        mergingEnd_idx=(mergingEnd_idx+endCared)%(len(self.ald['end'])+1)
                    
                    # find newly formed 2-mer's position
#                    pdb.set_trace()
                    if endCared: #right
                        newOligoMerging=np.array([mergingPos-1,mergingPos])%self.params['N']
                        #old position is used to get idx from state[0] and calculate rate
                        newOligoMergingOldPos=np.array([mergingPos-2,mergingPos])%self.params['N']
                        # insert new end position emerged from merged two monomers forming 2-mer.
                        self.ald['end']=np.insert(self.ald['end'],mergingEnd_idx,mergingPos)
                        self.ald['end']=np.insert(self.ald['end'],mergingEnd_idx,(mergingPos-leftRight)%self.params['N'])
                    else: 
                        if oligoNonMoving_idx<oligoMerging_idx:
                            # in this case, our oligoNonMoving is spanning cross circular boundary
                            mergingEnd_idx = len(self.ald['end'])
                        newOligoMerging=np.array([mergingPos,mergingPos+1])%self.params['N']
                        newOligoMergingOldPos=np.array([mergingPos,mergingPos+2])%self.params['N']
                        # insert new end position emerged from merged two monomers forming 2-mer.
                        self.ald['end']=np.insert(self.ald['end'],mergingEnd_idx,(mergingPos-leftRight)%self.params['N'])
                        self.ald['end']=np.insert(self.ald['end'],mergingEnd_idx,mergingPos)
                        
                    # calculate new alpha
                    newAlpha,_=self.aEndCalc([newOligoMergingOldPos],self.state,self.params)
                    
                    #they should be the same. same pair.
                    self.al['end']=np.insert(self.al['end'],mergingEnd_idx,newAlpha[1-endCared])
                    self.al['end']=np.insert(self.al['end'],mergingEnd_idx,newAlpha[endCared])
                    
                else:
                    # oligoMerging is multimer
                    # substitute mergingEnd
                
                    # mergingEnd_idx = movingEnd_idx -1 if left, +0 if right.
                    mergingEnd_idx = (movingEnd_idx+endCared-1)%len(self.ald['end'])
                    self.ald['end'][mergingEnd_idx]=(self.ald['end'][mergingEnd_idx]-leftRight)%self.params['N']
                    
                    if endCared: # right
                        newOligoMerging=np.concatenate(([pos+1],oligoMerging))
                        newOligoMergingOldPos=np.concatenate(([pos],oligoMerging))
                    else: 
                        newOligoMerging=np.concatenate((oligoMerging,[pos-1]))
                        newOligoMergingOldPos=np.concatenate((oligoMerging,[pos]))
                    newAlpha,_=self.aEndCalc([newOligoMergingOldPos],self.state,self.params)
                    self.al['end'][mergingEnd_idx]=newAlpha[1-endCared]
                    
                    
                ##### DIFFUSION UPDATE #####
                # substitute nonmoving, sustitute merging
                self.ald['diff']['pos'][oligoNonMoving_idx]=np.delete(self.ald['diff']['pos'][oligoNonMoving_idx],-endCared)
                if endCared: #right
                    self.ald['diff']['pos'][oligoMerging_idx]=np.append([(mergingPos-leftRight)%self.params['N']],self.ald['diff']['pos'][oligoMerging_idx])
                else:
                    self.ald['diff']['pos'][oligoMerging_idx]=np.append(self.ald['diff']['pos'][oligoMerging_idx],[(mergingPos-leftRight)%self.params['N']])

            # update states

            self.state=[x.copy() for x in self.newState]
            
            
            
            # this is when right end broke and it is 0. This makes the nonmoving oligo the largest index oligo.
            if oligoMerging_idx.size==0 and not endCared and oligoNonMoving[-1]<oligoNonMoving[0]:
                self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],-1)
                self.al['diff']['pos']=np.roll(self.al['diff']['pos'],-2)
            
            
            
            ##### ROLLING ######
            if pos==0:
                if endCared:
                    self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],-1)
                    self.al['diff']['pos']=np.roll(self.al['diff']['pos'],-2)
                    if len(oligoNonMoving)>2:
                        # we also need to roll ald['end']
                        self.ald['end']=np.roll(self.ald['end'],-2)
                        self.al['end']=np.roll(self.al['end'],-2)
                else:#move to the left
                    if oligoMerging_idx.size==0:
                        self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],-1)
                        self.al['diff']['pos']=np.roll(self.al['diff']['pos'],-2)
                    elif len(oligoMerging)==1:
                        self.ald['end']=np.roll(self.ald['end'],-2)
                        self.al['end']=np.roll(self.al['end'],-2)
                    else:
                        # merging with multimer
                        pass

            if movingPos==0:
                if endCared:
                    if oligoMerging_idx.size==0:
                        self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],1)
                        self.al['diff']['pos']=np.roll(self.al['diff']['pos'],2)
                    elif len(oligoMerging)==1:
                        if len(oligoNonMoving)<3:
                            self.ald['end']=np.roll(self.ald['end'],2)
                            self.al['end']=np.roll(self.al['end'],2)
                    else:
                        pass
                else: # LEFT
                    if oligoMerging_idx.size==0:
                        pass
                    elif len(oligoMerging)==1:
                        self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],1)
                        self.al['diff']['pos']=np.roll(self.al['diff']['pos'],2)
                    else:
                        self.ald['diff']['pos']=np.roll(self.ald['diff']['pos'],1)
                        self.al['diff']['pos']=np.roll(self.al['diff']['pos'],2)
                        self.ald['end']=np.roll(self.ald['end'],2)
                        self.al['end']=np.roll(self.al['end'],2)
            
            if mergingPos==0:
                if endCared:
                    if oligoMerging_idx.size==0:
                        pass
                    elif len(oligoMerging)==1:
                        if len(oligoNonMoving)<3:
                            self.ald['end']=np.roll(self.ald['end'],2)
                            self.al['end']=np.roll(self.al['end'],2)
                    else:
                        pass
                else: # LEFT
                    pass
#            self.state[0][idx]=(self.state[0][idx]+leftRight)%self.params['N']
#        #make ald['diff']['pos'] correctly ordered
#        order=np.argsort([x[-1] for x in self.ald['diff']['pos']])
#        self.ald['diff']['pos']=[self.ald['diff']['pos'][x] for x in order]
        #update ald['diff']['leftRight']
        self.ald['diff']['leftRight']=np.tile([0,1],len(self.ald['diff']['pos']))
        if type(np.unique(self.state[0])==np.sort(self.state[0]))=='bool':
            pdb.set_trace()
#        if 0 in np.diff(self.ald['hydro']):
#            pdb.set_trace()
#            aa=1
        return event
    
    def step(self):
        event,time=self.drawEvent()
        print(event)
        print(self.state[0])
        print(self.ald)
        print(self.al)
        event=self.execEvent(event,time)
        self.multimerIdx=np.where(np.array(map(len,self.ald['diff']['pos']))>1)[0]
        self.event=event
        self.time+=time
        self.stateSort()
        print(self.al)
        print(self.ald)
        print(self.state[0])
    
    def out(self):
        dic={'pos':self.state[0],'domains':self.state[1],'form':self.state[2],'event':self.event,'time':self.time,'ald_diff_pos':self.ald['diff']['pos'],'ald_end':self.ald['end'],'ald_hydro':self.ald['hydro'],'al_diff':self.al['diff'],'al_end':self.al['end'],'al_hydro':self.al['hydro']}
        return dic
    
    def stateSort(self):
        sortIdx=np.argsort(self.state[0])
        formPos=np.where(self.state[2])[0]
        tempFunc=lambda x: x[sortIdx]
        self.state[:2]=map(tempFunc,self.state[:2])
        self.newState=[x.copy() for x in self.state]
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
#            start=self.getIdx([segments[-1][0]])[0]
            
#            permFunc=lambda x:np.concatenate((x[start:],x[:start]))
#            newState=map(permFunc,state)
#            self.newState=newState
#            reverseStart=self.params['m']-start
#            self.newStateIdx=np.concatenate((self.newStateIdx[reverseStart:],self.newStateIdx[:reverseStart]))
            #did not modify state[0] for future getIdx usage, such as in constructing aldEnds
            
            joined_ends=np.concatenate((segments[-1],segments[0])) 
            #[0,1,2] and [N-2,N-1] becomes [B-2,N-1,0,1,2]
            segments[0]=joined_ends
            del(segments[-1])
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
        posInMultimer=np.concatenate(np.array(aldDiff['pos'])[self.multimerIdx])
        
        posATP=np.where(self.proteinForm_1hot==0)[0] #0 is ATP
        togetherMask=[x in posATP for x in posInMultimer]
        posAllowed=posInMultimer[np.where(togetherMask)]
        aldHydro={}
        aldHydro['pos']=posAllowed
        temp=np.zeros(self.params['N'])
        temp[posAllowed]=1
        aldHydro['1hot']=temp
        alHydro=np.ones(sum(aldHydro['1hot']))*k_h
        
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
        
    
    def checkATP(self,state,form_1hot):
        #turn all ADP form monomers into ATP
        # 0 is ATP state. 1 is ADP state
        hydro=form_1hot
        for t in xrange(self.params['N']):
            if hydro[t]==1:
                if (t+1)%self.params['N'] not in state[0] and (t-1)%self.params['N'] not in state[0]:
                    hydro[t]=0
        hydro_pos=hydro[state[0]]
        return [state[0],state[1],hydro_pos],hydro


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
    NN=1e5
    firstDic=aa.out()
    trace=pd.DataFrame(columns=firstDic.keys())
    trace.append([firstDic],ignore_index=True)
    for n in xrange(int(NN)):
        aa.step()
        trace=trace.append([aa.out()],ignore_index=True)
        trace.to_csv('tempTest.csv')
    pdb.set_trace()
    
#    aa.drawEvent()
#    aa.execEvent({'diff':0},2.)
#    aa.execEvent({'diff':1},2.)
    
    
    print(aa.state[0])
    print(aa.ald)
    print(aa.al)
#    aa.ald['hydro']=np.array([15, 16, 17, 18, 19,  0, 10, 11])
    aa.execEvent({'diff':3},2.) 
    aa.stateSort()
    print(aa.al)
    print(aa.ald)
    print(aa.state[0])
