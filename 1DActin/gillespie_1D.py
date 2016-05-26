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
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pdb
import pickle
#import cProfile
import sys
from multiprocessing import Pool
import pandas as pd
import time
import matplotlib.pyplot as plt

class unit:
    def __init__(self,pos,domain,form,params):
        self.pos=pos
        self.domain=domain
        self.form=form
        self.is_subunit=-1
        self.params=params
        
    def nb(self,leftRight,units):
        pos=np.array([x.pos for x in units])
        idx=np.where(pos==(self.pos+leftRight)%self.params['N'])[0]
        if len(idx)==0:
            return -1
        else:
            return units[idx[0]]
            
    def hydro(self):
        assert self.form==0
        assert self.is_subunit!=-1
        self.form=1
        motherOligo=self.is_subunit
        motherOligo.ends_rate=motherOligo.ends_rate_calc(motherOligo.subunits)
        
    def diffuse(self,leftRight,units):
        form_oligo=False
        # if form oligo or join oligo, then form_oligo becomes a dictionary
        self.pos=(self.pos+leftRight)%self.params['N']
        nb=self.nb(leftRight,units)
        if nb!=-1:
            if nb.is_subunit==-1:
                #monomer+monomer
                form_oligo={}
                form_oligo['monomer']=nb
                # is used for if form_oligo['monomer']
                form_oligo['oligo']=oligo([self,nb],self.params)
            else:
                #monomer+multimer
                form_oligo={}
                form_oligo['monomer']=0
                form_oligo['oligo']=nb.is_subunit
                nb.is_subunit.attach([self])
        return form_oligo
    
    def energize(self):
        # as we energize everytime a monomer is created, asserting no longer make sense.
#        assert self.form==1 
        self.form=0.
    
class oligo:
    def __init__(self,subunits,params):
        # params should be the same as params from intSpace
        # params contains stuff that don't change, such as energy function and parameters
    
        self.params=params
        self.subunits=subunits
        self.enslave(self.subunits)
        self.ends=self.findEnd(self.subunits)
        self.diff_rate=self.diff_rate_calc()
        self.ends_rate=self.ends_rate_calc(self.subunits)
#        self.history=np.empty(0,dtype=object)
        self.hash=np.random.rand()
        
    def snapshot(self,time):
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
        return dic
        
    def history_entry(self):
        self.history=np.append(self.history,self.snapshot)
        
    def attach(self,attach_units):
        # oligo_units should be a list, even if it's a monomer
        self.subunits=np.concatenate((self.subunits,attach_units))
        self.enslave(attach_units)
        self.ends=self.findEnd(self.subunits)
        self.diff_rate=self.diff_rate_calc()
        self.ends_rate=self.ends_rate_calc(self.subunits)
        
    def end_break(self,leftRight,units):
        #left end break or right end break
        
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
        form_oligo=movingUnit.diffuse(leftRight,units)
        return oligo_vanish,form_oligo


    def diffuse(self,leftRight,units):
        oligo_vanish=False
        lr=-(leftRight+1)/2
        movEnd=self.ends[lr]
        for x in self.subunits:
            x.pos=(x.pos+leftRight)%self.params['N']
        
        nb=movEnd.nb(leftRight,units)
        
        if nb!=-1:
            if nb.is_subunit==-1:
                #multimer+monomer
                oligo_vanish=nb
                self.attach([nb])
            else:
                #multimer+multimer
                oligo_vanish=self
                nb.is_subunit.attach(self.subunits)
        
        return oligo_vanish

    def enslave(self,subunits):
        for x in subunits:
            x.is_subunit=self
        
    def setFree(self,subunits):
        for x in subunits:
            x.is_subunit=-1
    
    def findEnd(self,subunits):
        leftEnd=-1
        rightEnd=-1
        for x in subunits:
            if x.nb(1,subunits)==-1:
                rightEnd=x
            elif x.nb(-1,subunits)==-1:
                leftEnd=x
        if leftEnd==-1 or rightEnd==-1:
            pdb.set_trace()
        return [leftEnd,rightEnd]
        
    def diff_rate_calc(self):
        alDiff=1/np.sqrt(len(self.subunits))*self.params['D_mono']*self.params['T']
        return alDiff 
    
    def ends_rate_calc(self,units):
        leftEnd=self.ends[0]
        leftEndNb=leftEnd.nb(1,units)
        rightEnd=self.ends[1]
        rightEndNb=rightEnd.nb(-1,units)
        if rightEndNb==-1 or leftEndNb==-1:
            pdb.set_trace()
        endDomains=[[leftEnd.domain[1],leftEndNb.domain[0]],[rightEndNb.domain[1],rightEnd.domain[0]]]
        endForms=[[leftEnd.form,leftEndNb.form],[rightEndNb.form,rightEnd.form]]
        delE=self.params['H'](endDomains,endForms)
        return np.exp(delE/self.params['T'])
        
        
def H_actin(endDomainVec,endFormVec,eps=-10.,xi=-1.):
    # a=(pos,form, domain type, protein type)
    # form: 1 is ADP, 0 is ATP.
    # ORDER MATTERS. (a,b) is an ORDERED pair, meaning protein a is left, b is right.
    #   side of interaction. This is important when distance is calculated.

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
    def __init__(self,init):
        params=init['params']
        self.params=params        
        if init['mode']=='unit':
            self.units=init['units']
        else:
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

    def check_monomer_atp(self):
        for x in self.monomers:
            x.form=0.

    def gen_oligos(self,units):
        def consecutive(pos):
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
        # diff, ald is oligos; end, ald is oligos; hydro, ald is units
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
        # diff idx is along list of oligos
        oligoMoving=self.ald['diff'][eventIdx/2]
        leftRight=eventIdx%2*2-1
        
        if isinstance(oligoMoving,oligo):
            oligo_vanish=oligoMoving.diffuse(leftRight,self.units)
            if oligo_vanish:
                if isinstance(oligo_vanish,oligo):
                    idx=np.where([oligo_vanish is x for x in self.oligos])[0]
                    self.oligos=np.delete(self.oligos,idx)
                else:
                    idx=np.where([oligo_vanish is x for x in self.monomers])[0]
                    self.monomers=np.delete(self.monomers,idx)
            # if no merging, then no modification needed.
        else:
            form_oligo=oligoMoving.diffuse(leftRight,self.units)
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
#        pdb.set_trace()
    
    def execHydro(self,eventIdx):
        try:
            self.ald['hydro'][eventIdx].hydro()
        except IndexError: 
            pdb.set_trace()
    
    def execEnd(self,eventIdx):
        oligoEndBreak=self.ald['end'][eventIdx/2]
        leftRight=eventIdx%2*2-1
        lr=-(leftRight+1)/2
        unitMoving=oligoEndBreak.ends[lr]
#        pdb.set_trace()
        oligo_vanish,form_oligo=oligoEndBreak.end_break(leftRight,self.units)
#        pdb.set_trace()
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
        return all([(x.ends_rate==x.ends_rate_calc(x.subunits)).all() for x in self.oligos])
        
    def check_is_subunit(self):
        tf=[]
        for x in self.units:
            if x.is_subunit!=-1:
                temp=[x is y for y in x.is_subunit.subunits]
                tf.append(sum(temp)==1)
                if sum(temp)!=1:
                    pdb.set_trace()
            
        return not all(tf)
                
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


#################################################################


class intSpaceSim():
    def __init__(self,intSpace,params={'mode':'nstep','time':1e3,'nstep':10000}):
        self.intSpace=intSpace
        self.params=params
    
    @staticmethod
    def event_type_func(event_type):
            if event_type=='diff':
                return 0
            elif event_type=='end':
                return 1
            elif event_type=='hydro':
                return 2

    def sim(self):
        trace=pd.DataFrame(0,index=np.arange(self.params['nstep']),columns = ['event_type','time','current_time'])
        oligoData=pd.DataFrame(np.nan,index=np.arange(5*self.params['nstep']),columns = ['hash', 'time', 'pos_leftEnd','pos_rightEnd','len','form_leftEnd','form_rightEnd','form_D','form_T','form_leftEnd2','form_rightEnd2'])
        monoData=pd.DataFrame(0,index=np.arange(self.params['nstep']),columns = ['len','form_D','time'])
#        oligoData.loc[0,:]={'hash':np.nan, 'time':np.nan, 'pos_leftEnd':np.nan,'pos_rightEnd':np.nan,'len':np.nan,'form_leftEnd':np.nan,'form_rightEnd':np.nan,'form_D':np.nan,'form_T':np.nan,'form_leftEnd2':np.nan,'form_rightEn2':np.nan}
#        monoData.loc[0,:]={'len':np.nan,'form_D':np.nan,'time':np.nan}
#        trace.loc[0,:]={'event_type':np.nan,'time':np.nan,'current_time':np.nan}
        oligoData=[]
        monoData=[]
        trace=[]
        if self.params['mode']=='nstep':
            for t in xrange(self.params['nstep']):
                self.intSpace.step()
                
                for oligo in self.intSpace.oligos:
#                    oligoData=oligoData.append(oligo.snapshot(self.intSpace.current_time),ignore_index=True)
#                    oligoData.loc[oligoData_counter]=oligo.snapshot(self.intSpace.current_time)
                    oligoData.append(oligo.snapshot(self.intSpace.current_time))
                mono_form=[x.form for x in self.intSpace.monomers]
#                monoData.loc[t]=               
                monoData.append({'len':len(self.intSpace.monomers),'form_D':sum(mono_form),'time':self.intSpace.current_time})
#                trace.loc[t]={'event_type':intSpaceSim.event_type_func(self.intSpace.event.keys()[0]),'time':self.intSpace.time,'current_time':self.intSpace.current_time}
                trace.append({'event_type':intSpaceSim.event_type_func(self.intSpace.event.keys()[0]),'time':self.intSpace.time,'current_time':self.intSpace.current_time})
        else: raise ValueError('mode in params of simulation is not valid')
        oligoData=pd.DataFrame(oligoData)
        monoData=pd.DataFrame(monoData)
        trace=pd.DataFrame(trace)
#        oligoData=oligoData.dropna()
        return oligoData,monoData,trace
    
    @staticmethod
    def save(fN,oligoData):
        pd.DataFrame.to_csv(oligoData,fN+'.csv')
    
    def posCalc(self,oligoData):
        if isinstance(oligoData,str):
            oligoData=pd.DataFrame.from_csv(oligoData+'.csv')
        pos=[]
        for i,row in oligoData.iterrows():
            N=self.intSpace.params['N']
            lE=row['pos_leftEnd']
            rE=row['pos_rightEnd']
            if lE>rE:
                pos.append(((lE+rE-N)/2.)%N)
            else:
                pos.append((lE+rE)/2.)
        
        oligoData['pos']=pos
        return oligoData
    
    @staticmethod
    def oligoGraph(N,fN,oligoData,graphParams={'truncate':False,'cutoffTime':1.,'cutoffNStep':10,'mode':'centered'}):
        # mode: centered, not_centered, average
        # oligo position as a function of time, each oligo is a trace
        if isinstance(oligoData,str):
            oligoData=pd.DataFrame.from_csv(oligoData+'.csv')
        
        fig=plt.subplot()
        ax=fig.axes
        oligos_h=np.unique(oligoData['hash'])
        for h in oligos_h:
            tempDic=oligoData.loc[oligoData['hash']==h]
            pos=tempDic['pos']
            tEnd=tempDic['time'][tempDic.index[-1]]
            tStart=tempDic['time'][tempDic.index[0]]
            dur=tEnd-tStart
            if graphParams['truncate']:
                if (dur<graphParams['cutoffTime']) or len(pos)<graphParams['cutoffNStep']:
                    continue
            diff=np.diff(pos)
            pos_diff=[]
            
            # pos_diff is the difference in average position of the oligo between time steps
            # here the circular space is linearized for displacement calculations.
            for d in diff:
                if abs(d)>N-2:
                    # it goes from one side of circular boundary to another side.
                    d-=np.sign(d)*N
                pos_diff.append(d)
            
            #pos_cum is cumulated displacement from position of creation
            pos_cum=np.zeros(len(pos_diff)+1)
            pos_cum[1:]=np.cumsum(pos_diff)
            if graphParams['mode']=='centered':
                xAxis=tempDic['time']-tStart
                plt.plot(xAxis,pos_cum,'k')
                leftEndPos=pos_cum-tempDic['len']/2.
                rightEndPos=pos_cum+tempDic['len']/2.
                plt.fill_between(xAxis,leftEndPos,rightEndPos,color='k',alpha=0.3)
                leftDMask=tempDic['form_leftEnd']==1
                leftTMask=np.logical_not(leftDMask)
                rightDMask=tempDic['form_rightEnd']==1
                rightTMask=np.logical_not(rightDMask)
                leftD2Mask=tempDic['form_leftEnd2']==1
                leftT2Mask=np.logical_not(leftDMask)
                rightD2Mask=tempDic['form_rightEnd2']==1
                rightT2Mask=np.logical_not(rightDMask)
                color_T='r'
                color_D='b'
                markerSize=0.6
                plt.scatter(xAxis[leftDMask],leftEndPos[leftDMask],color=color_D,s=markerSize,alpha=0.3)
                plt.scatter(xAxis[leftTMask],leftEndPos[leftTMask],color=color_T,s=2*markerSize,alpha=1)
                plt.scatter(xAxis[rightDMask],rightEndPos[rightDMask],color=color_D,s=markerSize,alpha=0.3)
                plt.scatter(xAxis[rightTMask],rightEndPos[rightTMask],color=color_T,s=2*markerSize,alpha=1)
                plt.scatter(xAxis[leftD2Mask],leftEndPos[leftD2Mask]+1,color=color_D,s=markerSize,alpha=0.3)
                plt.scatter(xAxis[leftT2Mask],leftEndPos[leftT2Mask]+1,color=color_T,s=2*markerSize,alpha=1)
                plt.scatter(xAxis[rightD2Mask],rightEndPos[rightD2Mask]-1,color=color_D,s=markerSize,alpha=0.3)
                plt.scatter(xAxis[rightT2Mask],rightEndPos[rightT2Mask]-1,color=color_T,s=2*markerSize,alpha=1)
                ax.set_xlabel('time since born (AU)')
                ax.set_ylabel('displacement since born')
            elif graphParams['mode']=='not_centered':
                plt.plot(tempDic['time'],pos_cum+pos[pos.index[0]],'k')
                ax.set_xlabel('time')
                ax.set_ylabel('pos')
            else:
                raise ValueError('not valid mode option in graphParams')
        plt.savefig(fN+'_oligoGraph'+'_'+graphParams['mode']+'.pdf')
    
    @staticmethod
    def runningCalc(fN,summaryData,npts=1000):
        if isinstance(summaryData,str):
            summaryData=pd.DataFrame.from_csv(summaryData+'.csv')
        running={'v':[],'vl':[],'l':[],'v(l)':[],'atp':[],'atp(l)':[]}
        
        t=10.
        selected=[]
        while len(selected)==0:
            t+=1
            selected=[(v,min(t-t_start,t_total),l,adp) for v,t_start,t_total,l,adp in zip(summaryData['v'],summaryData['t_start'],summaryData['t_total'],summaryData['l'],summaryData['form_D']) if t>t_start]
        
        t_init=t
        times=np.linspace(t_init,max(summaryData['t_end']),npts)
        for t in times:
            selected=[(v,min(t-t_start,t_total),l,adp) for v,t_start,t_total,l,adp in zip(summaryData['v'],summaryData['t_start'],summaryData['t_total'],summaryData['l'],summaryData['form_D']) if t>t_start]
            running['v'].append(np.sum([x[0]*x[1] for x in selected])/t)
            lt=np.sum([x[1]*x[2] for x in selected])
            running['l'].append(lt/t)
            vl=np.sum([x[0]*x[1]*x[2] for x in selected])
            running['vl'].append(vl/t)
            running['v(l)'].append(vl/lt) # length weighted velocity time-course
            atp=np.sum([x[1]*(1-x[3])*x[2] for x in selected])
            running['atp(l)'].append( atp/lt)
            running['atp'].append(atp/t)

        for label in ['v','vl','l','v(l)','atp','atp(l)']:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.plot(times,running[label],'-b',lw=2,label=label)
            ax.set_xlabel('time')
            ax.set_ylabel(label)
            plt.savefig(fN+'_running_'+label+'.pdf')
            
    
    @staticmethod
    def averageCalc(oligoData,N_intSpace,take_latter_half=True):
        if isinstance(oligoData,str):
            oligoData=pd.DataFrame.from_csv(oligoData+'.csv')
        
        if take_latter_half:
            l=len(oligoData)
            data=oligoData.loc[oligoData.index[l/2:]]
        oligos_h=np.unique(data['hash'])
        oligoSummary=pd.DataFrame(index=oligos_h,columns = ['v','v_sd','l','l_sd','t_start','t_end','t_total','form_left_D','form_right_D','form_D','form_left_D_2','form_right_D_2'])
        for h in oligos_h:
            tempDic=oligoData.loc[oligoData['hash']==h]
            if len(tempDic)<2:
                continue
            t_start=tempDic['time'][tempDic.index[0]]
            t_end=tempDic['time'][tempDic.index[-1]]
            pos=tempDic['pos']
            leng=tempDic['len']
            diff=np.diff(pos)
            pos_diff=[]
            N=N_intSpace
            for d in diff:
                if abs(d)>N-2:
                    d-=np.sign(d)*N
                pos_diff.append(d)
            l=len(pos_diff)+1
            pos_cum=np.zeros(l)
            pos_cum[1:]=np.cumsum(pos_diff)
            t_total=t_end-t_start
#            if t_total<1: # Get rid of noise
#                continue
            timeDiff=np.diff(tempDic['time'])
            v=pos_diff/timeDiff
            v_mean=pos_cum[-1]/t_total
            leng_mean=np.sum(leng[:-1]*timeDiff)/t_total
            
#            if leng_mean<5: # get rid of too short oligos
#                continue
            v_sd=np.std(v)
            leng_sd=np.std(leng)
            #form_left is time-weighted prob of left end to be ATP
            form_left_D=np.sum(tempDic['form_leftEnd'][:-1]*timeDiff)/t_total
            form_right_D=np.sum(tempDic['form_rightEnd'][:-1]*timeDiff)/t_total
            form_left_D_2=np.sum(tempDic['form_leftEnd2'][:-1]*timeDiff)/t_total
            form_right_D_2=np.sum(tempDic['form_rightEnd2'][:-1]*timeDiff)/t_total
            form_D = np.sum((tempDic['form_D']/leng)[:-1]*timeDiff)/t_total
            oligoSummary.loc[h]={'v':v_mean,'l':leng_mean,'v_sd':v_sd,'l_sd':leng_sd,'t_start':t_start,'t_end':t_end,'t_total':t_total,'form_left_D':form_left_D,'form_right_D':form_right_D,'form_D':form_D,'form_left_D_2':form_left_D_2,'form_right_D_2':form_right_D_2}
        return oligoSummary
    
    @staticmethod
    def averageGraph(fN):
        oligoSummary=pd.DataFrame.from_csv(fN+'.csv')
        leng=oligoSummary['l']
        v=oligoSummary['v']
        t=oligoSummary['t_total']
        v_mean_weighted=np.sum(v*leng*t)/np.sum(leng*t)
        
        nCol=100
        v_count,v_cenc=np.histogram(v,nCol)
        l_count,l_cenc=np.histogram(leng,nCol)
#        t_count,t_cenc=np.histogram(t,nCol,range=(1,120))
        t_count,t_cenc=np.histogram(t,nCol)
        
        func_v=lambda idx:(v_cenc[idx]+v_cenc[idx+1])/2.
        func_l=lambda idx:(l_cenc[idx]+l_cenc[idx+1])/2.
        func_t=lambda idx:(t_cenc[idx]+t_cenc[idx+1])/2.
        v_xaxis=map(func_v,range(nCol))
        l_xaxis=map(func_l,range(nCol))
        t_xaxis=map(func_t,range(nCol))
        fig=plt.figure()
#        ax=fig.add_subplot(111)
#        ax.set_ylabel('probability')
        ax1=fig.add_subplot(311)
        plt.plot(v_xaxis,v_count/float(sum(v_count)),'-k',lw=2)
#        ax1.set_xlabel('mean velocity')
        ax1.legend(['velocity'])
        ax1.set_title('Distribution, v_wm='+str(v_mean_weighted))
#        plt.savefig(fN+'_oligoSummary'+'_v.pdf')
        
        ax2=fig.add_subplot(312)
        plt.plot(l_xaxis,l_count/float(sum(l_count)),'-k',lw=2)
        ax2.legend(['length'])
        ax2.set_ylabel('probability')
        
        ax3=fig.add_subplot(313)
        ax=fig.axes
        plt.plot(t_xaxis,t_count/float(sum(t_count)),'-k',lw=2)
        ax3.legend(['lifetime'])
        ax3.set_xlabel('lifetime')
        
        plt.savefig(fN+'_graph.pdf')
    
    @staticmethod
    def monoGraph(fN_mono,fN_trace):
        monoData=pd.DataFrame.from_csv(fN_mono+'.csv')
        fig=plt.figure()
        ax=fig.add_subplot(111)
        plt.plot(monoData['time'],monoData['len'],'-k',lw=1)
        ll=np.array(monoData['len'])[:-1]
        runningMean=np.cumsum(ll*np.diff(monoData['time']))/monoData['time'][1:]
        plt.plot(monoData['time'][1:],runningMean,'-r',lw=2)
        ax.legend(['number of monomers','running average'])
        ax.set_title('Time Trajectory of Number of Monomers')
        ax.set_ylabel('number')
        ax.set_xlabel('time (AU)')
        ax.set_ylim([0.,10.])
        plt.savefig(fN+'_monoGraph.pdf')
    
    @staticmethod
    def statistics(trace):
        if isinstance(trace,str):
            trace=pd.DataFrame.from_csv(trace+'.csv')
        total_num=float(len(trace))
        diff_ratio=sum(trace['event_type']==intSpaceSim.event_type_func('diff'))
        end_ratio=sum(trace['event_type']==intSpaceSim.event_type_func('end'))
        hydro_ratio=len(trace)-diff_ratio-end_ratio
        print('diff:'+str(diff_ratio)+'/'+str(total_num)+', end:'+str(end_ratio)+'/'+str(total_num)+', hydro:'+str(hydro_ratio)+'/'+str(total_num))
#        pdb.set_trace()
    
    @staticmethod
    def timeWindowVelocity(oligoSummary):
        if isinstance(oligoSummary,str):
            oligoSummary=pd.DataFrame.from_csv(oligoSummary+'.csv')
    
    @staticmethod
    def treadmillingFilaments(oligoSummary):
        if isinstance(oligoSummary,str):
            oligoSummary=pd.DataFrame.from_csv(oligoSummary+'.csv')
        threshold=0.2
        dur_thr=100.
        left_mask=oligoSummary['form_left_D']>1-threshold
        right_mask=oligoSummary['form_right_D']<threshold
        dur_mask=oligoSummary['t_total']>dur_thr
        pdb.set_trace()
        all_mask=dur_mask & left_mask & right_mask
        temp=oligoSummary[all_mask]
        v_tread=np.sum(temp['v']*temp['l']*temp['t_total'])/np.sum(temp['l']*temp['t_total'])
        print v_tread
        return v_tread

def createSim(fN,NN=10000,xi=-1.,rateHydro=1e-3,N=50,m=30,eps=-100.):
    N=N
    k=2
    h=2
    m_array=np.array([[0,m],[0,0]])
#    xi=-0.1
    eps=eps # to make it essentially impossible
#    rateHydro=1e-1
    m=np.sum(m_array)
    params={'N':N,'m':m,'k':k,'h':h,
                'T':1.,'xi':None,'eps':None,
                'isCircular':True,
                'transRotRatio':1.,
                'rateHydro':rateHydro,
                'm_array':m_array,
                'dist_cutoff':2.,
                'D_mono':1.,
                'H':lambda x,y:H_actin(x,y,eps=eps,xi=xi)}


    position=np.random.choice(range(N),m,replace=False)
    proteinOrientation=np.repeat(np.reshape(np.ravel(np.indices((k,k)),1),[k**2,2]),params['m_array'].flatten(),axis=0)
    # [1,0]
    # for orientation, barbed end is 0, so barbed end is on the positive side. So velocity should be positive.
#    proteinForm=np.floor(np.random.rand(params['N'])*params['k'])
    proteinForm=np.ones(params['N'])
    mask=np.ones(params['N'])
    mask[position]=0
    proteinForm[np.where(mask)]=-1 #-1 indicate no protein, 0 is ATP, 1 is ADP
    proteinForm_pos=proteinForm[position]

    init={}
    init['pos']=position
    init['domain']=proteinOrientation
    init['form']=proteinForm_pos
    init['params']=params
    init['mode']='list'
    aa=InteractionSpace(init)
#    t=aa.timeSteps(NN=NN)
#    print t
    
    sim_params={}
    sim_params['time']=1e3
    sim_params['nstep']=NN
    sim_params['mode']='nstep'
    bb=intSpaceSim(aa,params=sim_params)

    oligoData,monoData,trace=bb.sim()
    oligoData=bb.posCalc(oligoData)
    intSpaceSim.save(fN+'_oligo',oligoData)
    intSpaceSim.save(fN+'_mono',monoData)
    intSpaceSim.save(fN+'_trace',trace)
    oligoSummary=intSpaceSim.averageCalc(oligoData,bb.intSpace.params['N'])
    intSpaceSim.save(fN+'_summary',oligoSummary)
    return fN


if __name__ == "__main__":
    N=50
    m=30
    NN=100000
    xi=-1.
    rateHydro=1e-1
    
    fN='N-'+str(N)+'_m-'+str(m)+'_NN-'+str(NN)+'_xi-'+str(xi)+'_rh-'+str(rateHydro)
#    createSim(fN,NN=NN,xi=xi,rateHydro=rateHydro)
    intSpaceSim.runningCalc(fN,fN+'_summary')
    intSpaceSim.oligoGraph(N,fN,fN+'_oligo',graphParams={'truncate':True,'cutoffTime':20.,'cutoffNStep':NN/2,'mode':'centered'})
    intSpaceSim.averageGraph(fN+'_summary')
    intSpaceSim.monoGraph(fN+'_mono',fN+'_trace')
    intSpaceSim.statistics(fN+'_trace')

#    v_tread=[]
#    for i in xrange(10):
#        createSim(fN,NN=NN,xi=-1.,rateHydro=1e-2)
#        v_tread.append(intSpaceSim.treadmillingFilaments(fN+'_summary'))
#    pdb.set_trace()
