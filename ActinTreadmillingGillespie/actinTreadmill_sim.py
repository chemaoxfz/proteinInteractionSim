# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 21:17:18 2015

@author: xfz

Simulation of actin treadmilling in a general 1D protein interaction setting.
"""

import numpy as np
import pdb
import time


class unit:
    # class for one unit, could be subunit or monomer
    def __init__(self,pos,domain,form,params):
        self.pos=pos
        self.domain=domain
        self.form=form
        self.is_subunit=-1 #whether this unit is subunit of an oligo or not. It will refer to the oligo object if it is a subunit.
        self.params=params
        
    def nb(self,leftRight,units):
        #find left or right neightboring unit among units
        pos=np.array([x.pos for x in units])
        idx=np.where(pos==(self.pos+leftRight)%self.params['N'])[0]
        if len(idx)==0:
            # no neighboring unit
            return -1
        else:
            return units[idx[0]]
            
    def hydro(self):
        #hydrolyze this unit. Note that the dynamics is set up such that only unit inside an oligo that is also ATP will be chosen to hydrolyze.
#        assert self.form==0
#        assert self.is_subunit!=-1
        self.form=1
        motherOligo=self.is_subunit
        # in case this unit is close to end, so update end rate.
        motherOligo.ends_rate=motherOligo.ends_rate_calc(motherOligo.subunits)
        
    def diffuse(self,leftRight,units):
        # diffusion action of this unit
        form_oligo=False
        # if form oligo or join oligo, then form_oligo becomes a dictionary
        self.pos=(self.pos+leftRight)%self.params['N']
        nb=self.nb(leftRight,units)
        event_code=0 #event code is used to capture statistics of what kind of interaction happened
        if nb!=-1:
            if nb.is_subunit==-1:
                #monomer+monomer
                event_code=1
                form_oligo={}
                form_oligo['monomer']=nb
                # is used for if form_oligo['monomer']
                form_oligo['oligo']=oligo([self,nb],self.params)
            else:
                #monomer+multimer
                event_code=2
                form_oligo={}
                form_oligo['monomer']=0
                form_oligo['oligo']=nb.is_subunit
                nb.is_subunit.attach([self])
        return form_oligo,event_code
    
    def energize(self):
        # energize corresponds to the unit becomes a monomer, and therefore instantly becomes ATP if it was ADP.
        # as we energize everytime a monomer is created, asserting no longer make sense.
#        assert self.form==1 
        self.form=0.
    
class oligo:
    # class for oligos, which correspond to >=2 units next to each other.
    def __init__(self,subunits,params):
        # params should be the same as params from intSpace
        # params contains stuff that don't change, such as energy function and parameters
        self.params=params
        self.subunits=subunits
        self.enslave(self.subunits)
        self.ends=self.findEnd(self.subunits)
        self.diff_rate=self.diff_rate_calc()
        self.ends_rate=self.ends_rate_calc(self.subunits)
        self.hash=np.random.rand()
        
    def snapshot(self,time):
        # this function is to facilitate output of useful data for recording.
        dic={}
        lE=self.ends[0].pos
        rE=self.ends[1].pos
        dic['pos_leftEnd']=lE
        dic['form_leftEnd']=self.ends[0].form
        dic['form_leftEnd2']=self.ends[0].nb(1,self.subunits).form
        dic['pos_rightEnd']=rE
        dic['form_rightEnd']=self.ends[1].form
        dic['form_rightEnd2']=self.ends[1].nb(-1,self.subunits).form
        dic['len']=len(self.subunits)
        forms=[x.form for x in self.subunits]
        dic['form_D']=sum(forms)
        dic['form_T']=dic['len']-dic['form_D']
        dic['time']=time
        dic['hash']=self.hash
        un=self.ends[1]
        form=[]
        for i in xrange(len(self.subunits)):
            form.append(un.form)
            un=un.nb(-1,self.subunits)
        return dic,form
        
    def attach(self,attach_units):
        # attach a unit to this oligo. Could be attaching another oligo (i.e. a list of subunits), or just one unit (a monomer).
        # oligo_units should be a list, even if it's a monomer
        self.subunits=np.concatenate((self.subunits,attach_units))
        self.enslave(attach_units)
        self.ends=self.findEnd(self.subunits)
        self.diff_rate=self.diff_rate_calc()
        self.ends_rate=self.ends_rate_calc(self.subunits)
        
    def end_break(self,leftRight,units):
        # end_break dynamics. or depolymerization.
        # left end break or right end break
        
        oligo_vanish=False
        
        lr=-(leftRight+1)/2
        movingUnit=self.ends[lr]
        if len(self.subunits)==2:
            self.setFree(self.subunits)
            oligo_vanish=True
        else:
            self.setFree([movingUnit])
            endIdx=np.where([movingUnit is x for x in self.subunits])[0]
            self.subunits=np.delete(self.subunits,endIdx)
            self.ends=self.findEnd(self.subunits)
            self.diff_rate=self.diff_rate_calc()
            self.ends_rate=self.ends_rate_calc(self.subunits)
        
        # above update oligo, below update movingUnit
        form_oligo,event_code=movingUnit.diffuse(leftRight,units)
        event_code=-3*lr+event_code
            
        return oligo_vanish,form_oligo,event_code


    def diffuse(self,leftRight,units):
        #diffusion dynamics
        oligo_vanish=False
        lr=-(leftRight+1)/2
        movEnd=self.ends[lr]
        for x in self.subunits:
            x.pos=(x.pos+leftRight)%self.params['N']
        
        nb=movEnd.nb(leftRight,units)
        event_code=4
        if nb!=-1:
            if nb.is_subunit==-1:
                #multimer+monomer
                oligo_vanish=nb
                self.attach([nb])
                event_code=5
            else:
                #multimer+multimer
                # let the longer one survive
                mergingOligo=nb.is_subunit
                if len(self.subunits)>=len(mergingOligo.subunits):
                    oligo_vanish=self
                    mergingOligo.attach(self.subunits)
                else:
                    oligo_vanish=mergingOligo
                    self.attach(mergingOligo.subunits)
                event_code=7
        
        return oligo_vanish,event_code

    def enslave(self,subunits):
        # a romantic way of saying making one unit a subunit of this oligo.
        for x in subunits:
            x.is_subunit=self
        
    def setFree(self,subunits):
        # romance is great.
        for x in subunits:
            x.is_subunit=-1
    
    def findEnd(self,subunits):
        # find, among its own subunits, what are the end subunits.
        leftEnd=-1
        rightEnd=-1
        for x in subunits:
            if x.nb(1,subunits)==-1:
                rightEnd=x
            elif x.nb(-1,subunits)==-1:
                leftEnd=x
#        if leftEnd==-1 or rightEnd==-1:
#            pdb.set_trace()
        return [leftEnd,rightEnd]
        
    def diff_rate_calc(self):
        # calculate the diffusion rate. Note that here I made the diffusion rate decrease by square-root of the length.
        alDiff=1/np.sqrt(len(self.subunits))*self.params['D_mono']*self.params['T']
        return alDiff 
    
    def ends_rate_calc(self,units):
        # calculate the rate of end-breaking events.
        leftEnd=self.ends[0]
        leftEndNb=leftEnd.nb(1,units)
        rightEnd=self.ends[1]
        rightEndNb=rightEnd.nb(-1,units)
#        if rightEndNb==-1 or leftEndNb==-1:
#            pdb.set_trace()
        endDomains=[[leftEnd.domain[1],leftEndNb.domain[0]],[rightEndNb.domain[1],rightEnd.domain[0]]]
        endForms=[[leftEnd.form,leftEndNb.form],[rightEndNb.form,rightEnd.form]]
        delE=self.params['H'](endDomains,endForms)
        return np.exp(delE/self.params['T'])
        
        
def H_actin(endDomainVec,endFormVec,eps=-10.,xi=-1.):
    # energy function that return the interaction energy between two proteins.
    # form: 1 is ADP, 0 is ATP.
    # domain: 0 is barbed, 1 is pointed.

    h=k=2
    H_array=np.zeros([h,h,k,k])
    H_array[0][0]=np.array([[0,eps],
                           [eps,0]])
    H_array[0][1]=np.array([[0,xi],
                           [eps,0]])
    H_array[1][0]=H_array[0][1].T
    H_array[1][1]=np.array([[0,xi],
                           [xi,0]])
    getE=lambda idx:H_array[endFormVec[idx][0]][endFormVec[idx][1]][endDomainVec[idx][0]][endDomainVec[idx][1]]
    return np.array(map(getE,xrange(len(endDomainVec))))
        

class InteractionSpace:
    # the space that we are simulating.
    def __init__(self,init):
        params=init['params']
        self.params=params        
        if init['mode']=='unit':
            # this corresponds to initialize with units objects provided.
            self.units=init['units']
        else:
            # this is the more common way to initialize, which is by providing position, domain, and form of the proteins.
            pos=init['pos']
            domain=init['domain']
            form=init['form']
            units=[]
            for idx in xrange(len(pos)):
                aa=unit(pos[idx],domain[idx],form[idx],params)
                units.append(aa)
            self.units=np.array(units)
        
        
        self.oligos,self.monomers=self.gen_oligos(self.units)
#        self.check_monomer_atp()
        self.al,self.ald=self.alConstruct()    
        self.event=None
        self.current_time=0.
        self.time=0.
        self.event_code=None
        # event code: for diff, 0 is mono no merging, 1 is mono + mono, 2 is mono+oligo (L end), 3 is mono+oligo (R end), 4 is oligo no merging, 5 is oligo+mono(L), 6 is oligo+mono(R) 7 is oligo+oligo.
        # for end, 0 is left end breaking, no merging, 1 if left merge mono, 2 is left merge oligo, 3,4,5 are corresponding right ones.

    def check_monomer_atp(self):
        # to check whether all the monomers are ATP.
        # this function is not really used. It's for debugging purpose.
        for x in self.monomers:
            x.form=0.

    def gen_oligos(self,units):
        # given a list of units, find what are the oligos in their and generate these oligo objects, as well as a list of left-over monomers.
        def consecutive(pos):
            # from a position array find which are consecutive to each other.
            endI=np.where(np.diff(pos)!=1)[0]
            segments=np.split(pos, endI+1)
            
            if (segments[0][0]-segments[-1][-1])%self.params['N']==1:
                #i.e. first and last segment is one segment
                joined_ends=np.concatenate((segments[-1],segments[0])) 
                #[0,1,2] and [N-2,N-1] becomes [B-2,N-1,0,1,2]
                segments[0]=joined_ends
                del(segments[-1])
            return segments
            
        pos=np.array([x.pos for x in units])
        pos_sorted=np.sort(pos)
        segments=consecutive(pos_sorted)
        oligos=[]
        monomers=[]
        for seg in segments:
            oligo_units=[units[np.where(pos==p)[0][0]] for p in seg]
            if len(oligo_units)>1:
                oligos.append(oligo(oligo_units,self.params))
            else:
                monomers.append(oligo_units[0])
        
        return np.array(oligos),np.array(monomers)
    
    def alConstruct(self):
        # al, is list of alpha's, the rate of each reaction used in gillespie algorithm.
        # ald is the description of everything that is needed to execute each reaction.
        # this function constructs al and ald.
        # ald's contents: diff, ald is oligos; end, ald is oligos; hydro, ald is units
        al={}
        ald={}
        aa=np.array([x.diff_rate for x in self.oligos]+[self.params['D_mono'] for x in self.monomers])
        al['diff']=np.tile(aa,[2,1]).T.flatten()
        ald['diff']=np.concatenate((self.oligos,self.monomers))
        al['end']=np.array([x.ends_rate for x in self.oligos]).flatten()
        ald['end']=self.oligos
        assert sum(map(len,[x.subunits for x in self.oligos]))==self.params['m']-len(self.monomers)
        ald['hydro']=[x for x in self.units if x not in self.monomers and x.form==0]
        al['hydro']=np.ones(len(ald['hydro']))*self.params['rateHydro']
        return al,ald
    
    def drawEvent(self):
        # using al, the rates, to draw the next event and corresponding time.
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
        return eventDic,time
    
    def execEvent(self,event,time):
        # execute the event.
        if event.keys()[0]=='diff':
            self.execDiff(event['diff'])
        elif event.keys()[0]=='end':
            self.execEnd(event['end'])
        elif event.keys()[0]=='hydro':
            self.execHydro(event['hydro'])
        else:
            raise ValueError('BAD:the event type is not one of diff, hydro, end')
        self.time+=time

    def execDiff(self,eventIdx):
        # execute a diffusion event.
        # diff idx is along list of oligos
        oligoMoving=self.ald['diff'][eventIdx/2] #the oligo (oligo/monomer) that is chosen
        leftRight=eventIdx%2*2-1
        
        if isinstance(oligoMoving,oligo):
            #this is an oligo
            oligo_vanish,event_code=oligoMoving.diffuse(leftRight,self.units)
            if event_code==5:
                self.event_code=(leftRight+1)/2+event_code
            if oligo_vanish:
                if isinstance(oligo_vanish,oligo):
                    idx=np.where([oligo_vanish is x for x in self.oligos])[0]
                    self.oligos=np.delete(self.oligos,idx)
                else:
                    idx=np.where([oligo_vanish is x for x in self.monomers])[0]
                    self.monomers=np.delete(self.monomers,idx)
            # if no merging, then no modification needed.
        else:
            # this is an monomer
            form_oligo,event_code=oligoMoving.diffuse(leftRight,self.units)
            if event_code==2:
                self.event_code=1-(leftRight+1)/2+event_code
            if form_oligo:
                # not empty
                mono=form_oligo['monomer']
                if mono:
                    # monomer+monomer
                    idx=np.where([x in [mono,oligoMoving] for x in self.monomers])[0]
                    self.monomers=np.delete(self.monomers,idx)
                    self.oligos=np.insert(self.oligos,0,form_oligo['oligo'])
                else:
                    # monomer+multimer
                    idx=np.where([oligoMoving is x for x in self.monomers])[0]
                    self.monomers=np.delete(self.monomers,idx)
                #empty, do nothing
        
       
    
    def execHydro(self,eventIdx):
        # execute hydrolysis event.
        self.ald['hydro'][eventIdx].hydro()
    
    def execEnd(self,eventIdx):
        # execute an end-breaking or depolymerization event.
        oligoEndBreak=self.ald['end'][eventIdx/2]
        leftRight=eventIdx%2*2-1
        lr=-(leftRight+1)/2
        unitMoving=oligoEndBreak.ends[lr]
        oligo_vanish,form_oligo,self.event_code=oligoEndBreak.end_break(leftRight,self.units)
        if form_oligo:
            # not empty
            mono=form_oligo['monomer']
            if mono:
                # monomer + monomer (mergeOligo)
                idx=np.where([x in [mono,unitMoving] for x in self.monomers])[0]
                self.monomers=np.delete(self.monomers,idx)
                self.oligos=np.insert(self.oligos,0,form_oligo['oligo'])
            else:
                # monomer + multimer (mergeOligo)
                idx=np.where([unitMoving is x for x in self.monomers])[0]
                self.monomers=np.delete(self.monomers,idx)
        else:
            #empty, add the end to monomers
            self.monomers=np.insert(self.monomers,0,unitMoving)
            unitMoving.energize()
        
        if oligo_vanish:
            idx=np.where([oligoEndBreak is x for x in self.oligos])[0]
            self.oligos=np.delete(self.oligos,idx)
            
            idx=np.where([unitMoving is not x for x in oligoEndBreak.subunits])[0]
            nonmoving_unit=oligoEndBreak.subunits[idx[0]]
            self.monomers=np.insert(self.monomers,0,nonmoving_unit)
            nonmoving_unit.energize()
                
    
    def step(self):
        # take a step in simulation
        event,time=self.drawEvent()
        self.execEvent(event,time)
#        if self.check_is_subunit():
#            pdb.set_trace()
        self.al,self.ald=self.alConstruct()
#        if not self.check_ends_rate():pdb.set_trace()
#        if self.current_time>6000: pdb.set_trace()
        self.event=event
        self.time=time
        self.current_time+=time

    def check_ends_rate(self):
        # makde sure the end rates is right. For debugging purpose.
        return all([(x.ends_rate==x.ends_rate_calc(x.subunits)).all() for x in self.oligos])
        
    def check_is_subunit(self):
        # for debugging purpose. Check that each subunit of an oligo have its is_subunit pointing to this oligo
        tf=[]
        for x in self.units:
            if x.is_subunit!=-1:
                temp=[x is y for y in x.is_subunit.subunits]
                tf.append(sum(temp)==1)
                if sum(temp)!=1:
                    pdb.set_trace()
            
        return not all(tf)
        
    # the following are some convenient functions for debugging or printing purposes.
    def disp_pos(self):
        pos=np.sort([x.pos for x in self.units])
        return pos
        
    def disp_mono(self):
        mono_pos=np.sort([x.pos for x in self.monomers])
        return mono_pos
        
    def disp_oligo(self):
        oligo_pos=[np.sort([x.pos for x in y.subunits]) for y in self.oligos]
        return oligo_pos
    
    def disp(self):
        pos=self.disp_pos
        mono=self.disp_mono
        oligo=self.disp_oligo
        return pos,mono,oligo
        
    def timeSteps(self,NN=100000):
        st=time.time()
        for n in xrange(NN):
            self.step()
        ed=time.time()
        return ed-st




