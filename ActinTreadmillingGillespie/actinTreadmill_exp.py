# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 11:45:42 2016

@author: chemaoxfz

Doing experiments and collect data from runs of actinTreadmill_sim.py.
"""

#################################################################
import matplotlib
matplotlib.use('Agg')
from actinTreadmill_sim import InteractionSpace
import pickle
from multiprocessing import Pool
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

class intSpaceSim():
    # a simulation class that gets all the methods needed for simulation and plotting in one place.
    def __init__(self,intSpace,params={'nstep':10000}):
        self.intSpace=intSpace
        self.params=params
    
    @staticmethod
    def event_type_func(event_type):
        # for easy encoding of event types.
        if event_type=='diff':
            return 0
        elif event_type=='end':
            return 1
        elif event_type=='hydro':
            return 2

    def sim(self):
        # execute the simulation and record data.
        oligoData=[]
        monoData=[]
        trace=[]
        form=[]
        for t in xrange(self.params['nstep']):
            self.intSpace.step()
            for oligo in self.intSpace.oligos:
                sp,fm=oligo.snapshot(self.intSpace.current_time)
                oligoData.append(sp)
                form.append(fm)
            mono_form=[x.form for x in self.intSpace.monomers]
            monoData.append({'len':len(self.intSpace.monomers),'form_D':sum(mono_form),'time':self.intSpace.current_time})
            trace.append({'event_type':intSpaceSim.event_type_func(self.intSpace.event.keys()[0]),'event_code':self.intSpace.event_code,'time':self.intSpace.time,'current_time':self.intSpace.current_time})
        oligoData=pd.DataFrame(oligoData)
        monoData=pd.DataFrame(monoData)
        trace=pd.DataFrame(trace)
        return oligoData,monoData,trace,form
    
    @staticmethod
    def save(fN,oligoData):
        pd.DataFrame.to_csv(oligoData,fN+'.csv')
    
    def posCalc(self,oligoData):
        # calculate the position of oligos and add it to oligoData
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
        # this is to plot time trajectory of oligos in this simulation.
        # to make sure the simulation has converged, I usually set cutoffNStep to be half of the total number of steps.
        # mode: centered, not_centered. 
        # centered means plotting the time=0 to be time when an oligo is born, and pos=0 to be the position the oligo is born.
        # not_centered means plotting time to be simulation time, position to be interaction space position, but the oligo's relative motion to its position of born is not circular, but linearize.
        # oligo position as a function of time, each oligo is a trace
        if isinstance(oligoData,str):
            # in this case, a file name is given.
            oligoData=pd.DataFrame.from_csv(oligoData+'.csv')
        
        fig=plt.subplot()
        ax=fig.axes
        oligos_h=np.unique(oligoData['hash'])
        # each oligo has a distinct hash, so we can distinguish each entry is for which oligo using the hash.
        for h in oligos_h:
            tempDic=oligoData.loc[oligoData['hash']==h]
            # tempDic contains all entries of the oligo with hash h.
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
                
                ax.set_xlabel('time since born (AU)')
                ax.set_ylabel('displacement since born')
            elif graphParams['mode']=='not_centered':
                xAxis=tempDic['time']
                plt.plot(xAxis,pos_cum+pos[pos.index[0]],'k')
                pos_cum=pos_cum+pos[pos.index[0]]
                ax.set_xlabel('time')
                ax.set_ylabel('pos')
            else:
                raise ValueError('not valid mode option in graphParams')
            
            # black line is mean position of an oligo
            # grey area is length of oligo
            # at each time, the two ends of an oligo and the next-to-end subunits' hydrolysis state are also plotted.
            # with red means ATP (hotter, more energy) and blue means ADP (colder, lower energy)
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
        plt.savefig(fN+'_oligoGraph'+'_'+graphParams['mode']+'.pdf')
    
    @staticmethod
    def runningCalc(fN,summaryData,npts=1000):
        # to calculate some running averages. Looking at these may help determine whether the simulation has converged to a steady state.
        if isinstance(summaryData,str):
            summaryData=pd.DataFrame.from_csv(summaryData+'.csv')
        
        # v is treadmilling velocity
        # vl, = v*l is treadmilling velocity times length of the oligo
        # l is oligo length
        # v(l) is length-weighted treadmilling velocity
        # atp is average number of subunits being atp
        # l_atp is average probability of left end (pointed end, - end) of an oligo being ATP.
        # r_atp is the same for right end (barbed end, + end)
        # l2_atp is the same for left next-to-end subunit being ATP
        running={'v':[],'vl':[],'l':[],'v(l)':[],'atp':[],'atp(l)':[],'l_atp':[],'r_atp':[],'l2_atp':[],'r2_atp':[]}
        
        # a lowest time to start counting averages
        t=10.
        selected=[]
        while len(selected)==0:
            # find the first non-empty entry
            t+=1
            selected=[(v,min(t-t_start,t_total),l,adp,l_adp,r_adp,l2_adp,r2_adp) for v,t_start,t_total,l,adp,l_adp,r_adp,l2_adp,r2_adp in zip(summaryData['v'],summaryData['t_start'],summaryData['t_total'],summaryData['l'],summaryData['form_D'],summaryData['form_left_D'],summaryData['form_left_D_2'],summaryData['form_right_D'],summaryData['form_right_D_2']) if t>t_start]
        
        t_init=t
        times=np.linspace(t_init,max(summaryData['t_end']),npts)
        for t in times:
            # one issue, t is not necessarily all time, as multiple filaments may coexist. i.e. sum(t for each filament till time t)>t
            selected=[(v,min(t-t_start,t_total),l,adp,l_adp,r_adp,l2_adp,r2_adp) for v,t_start,t_total,l,adp,l_adp,r_adp,l2_adp,r2_adp in zip(summaryData['v'],summaryData['t_start'],summaryData['t_total'],summaryData['l'],summaryData['form_D'],summaryData['form_left_D'],summaryData['form_left_D_2'],summaryData['form_right_D'],summaryData['form_right_D_2']) if t>t_start]
            running['v'].append(np.sum([x[0]*x[1] for x in selected])/t)
            lt=np.sum([x[1]*x[2] for x in selected])
            running['l'].append(lt/t)
            vl=np.sum([x[0]*x[1]*x[2] for x in selected])
            running['vl'].append(vl/t)
            running['v(l)'].append(vl/lt) # length weighted velocity time-course
            atp=np.sum([x[1]*(1-x[3])*x[2] for x in selected])
            running['atp(l)'].append( atp/lt)
            running['atp'].append(atp/t)
            l_atp,l2_atp,r_atp,r2_atp=[np.sum([(1-x[i])*x[1] for x in selected]) for i in [4,5,6,7]]
            running['l_atp'].append(l_atp/t)
            running['l2_atp'].append(l2_atp/t)
            running['r_atp'].append(r_atp/t)
            running['r2_atp'].append(r2_atp/t)

        for label in running.keys():
            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.plot(times,running[label],'-b',lw=2,label=label)
            ax.set_xlabel('time')
            ax.set_ylabel(label)
            plt.savefig(fN+'_running_'+label+'.pdf')
    
    @staticmethod
    def traceRunningPlot(fN,trace,npts=1000):
        # same as running plot, but for some quantities contained in trace file.
        if isinstance(trace,str):
            trace=pd.DataFrame.from_csv(trace+'.csv')
        t_end=max(trace['current_time'])
        # here I directly hard-codedly discarded data points before half of the duration of simulation
        times=np.linspace(t_end/2,t_end,npts)
        
        # kl_on is left (pointed, barbed) on-rate, i.e. rate of diffusion event that attach a monomer on left end of an oligo.
        # similar for kr_on, kl_off, kr_off
        trrun={'kl_on':np.zeros(npts),'kr_on':np.zeros(npts),'kl_off':np.zeros(npts),'kr_off':np.zeros(npts)}
        diff=intSpaceSim.event_type_func('diff')
        end=intSpaceSim.event_type_func('end')
        # dics is what pair of event_type and event_code would correspond to which rate's calculation.
        dics={(diff,3):'kl_on',(diff,6):'kl_on',(end,2):'kl_on',(diff,2):'kr_on',(diff,5):'kr_on',(end,5):'kr_on',(end,0):'kl_off',(end,1):'kl_off',(end,2):'kl_off',(end,3):'kr_off',(end,4):'kr_off',(end,5):'kr_off'}
        for _,row in trace.iterrows():
            try:
                trrun[dics[(row['event_type'],row['event_code'])]][times>row['current_time']]+=1
            except KeyError:
                continue

        for label in trrun.keys():
            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.plot(times,trrun[label]/times,'-b',lw=2,label=label)
            ax.set_xlabel('time')
            ax.set_ylabel(label)
            plt.savefig(fN+'_trrun_'+label+'.pdf')
                
    @staticmethod
    def oneValCalc(summary,trace,m,time_cutoff='half'):
        # overall summary of the whole simulation, essentially the last point of running plot.
        # used for varying parameter comparison.
        # two parts. one from the trace file. one from the summary data file.
        if isinstance(trace,str):
            trace=pd.DataFrame.from_csv(trace+'.csv')
        
        t_end=trace['current_time'][trace.index[-1]]
        if isinstance(time_cutoff,str):
            if time_cutoff=='half':
                time_cutoff=t_end/2
            else:
                raise ValueError('Invalid time_cutoff option')

        t_total=t_end-time_cutoff
        tr={'kl_on':0.,'kr_on':0.,'kl_off':0.,'kr_off':0.}
        diff=intSpaceSim.event_type_func('diff')
        end=intSpaceSim.event_type_func('end')
        dics={(diff,3):'kl_on',(diff,6):'kl_on',(end,2):'kl_on',(diff,2):'kr_on',(diff,5):'kr_on',(end,5):'kr_on',(end,0):'kl_off',(end,1):'kl_off',(end,2):'kl_off',(end,3):'kr_off',(end,4):'kr_off',(end,5):'kr_off'}
        trace_truncated=trace.loc[trace['current_time']>time_cutoff]
        for et,ec in zip(trace_truncated['event_type'],trace_truncated['event_code']):
            try:
                tr[dics[(et,ec)]]+=1
            except KeyError:
                continue
        for key in tr.keys():
            tr[key]=tr[key]/t_total
        
        
        if isinstance(summary,str):
            summary=pd.DataFrame.from_csv(summary+'.csv')
        
        # v(lt) is average treadmilling velocity both life-time and length weighted.
#        su={'v(lt)':0.,'v(t)':0.,'atp':0.,'l_atp':0.,'r_atp':0.,'l2_atp':0.,'r2_atp':0.}
        su={}
        v,t,l,adp,l_adp,r_adp,l2_adp,r2_adp = [summary[key] for key in ['v','t_total','l','form_D','form_left_D','form_left_D_2','form_right_D','form_right_D_2']]

        t_total=np.sum(t)
        lt=np.sum(t*l)
        su['v(t)']=np.sum(v*t)/t_total
        su['v(lt)']=np.sum([v*t*l])/t_total
        su['l']=lt/t_total
        su['atp']=np.sum(t*(1-adp)*l)/lt
        su['l_atp'],su['l2_atp'],su['r_atp'],su['r2_atp']=[np.sum((1-x)*t)/t_total for x in [l_adp,r_adp,l2_adp,r2_adp]]

        # calculate form and length. form is atp probability
        form_keys=['form_'+str(i+1) for i in xrange(m)]
        l_keys=['l_'+str(i+1) for i in xrange(m)]
        l_time_total=dict(zip(l_keys,[0.]*m))
        l_prob=l_time_total.copy()
        form_time_total=dict(zip(form_keys,[0.]*m))
        form_prob=form_time_total.copy()

        for lkey in l_keys:
            l_time_total[lkey]=np.sum(summary[lkey])
            l_prob[lkey]=l_time_total[lkey]/t_total
        i=0
        for fkey in form_keys:
            form_time_total[fkey]=np.sum(summary[fkey])
            form_prob[fkey]=1-form_time_total[fkey]/max(np.sum([l_time_total[x] for x in l_keys[i:]]),1e-7)
            i+=1
        
        su.update(form_prob)
        su.update(l_prob)
        su.update(tr)
        return pd.DataFrame(su,index=[0])
    
    @staticmethod
    def averageCalc(oligoData,form,N_intSpace,m_intSpace,take_latter_half=True):
        # calculate summary data for each oligo. (averaged over life time)
        # contains treadmilling velocity, start - end time, length, form at ends and next-to-ends
        # also contains probability for ATP or ADP form at each subunit position (relative to right end)
        # and probability for being in each length
        if isinstance(oligoData,str):
            oligoData=pd.DataFrame.from_csv(oligoData+'.csv')
        
        if take_latter_half:
            l=len(oligoData)
            data=oligoData.loc[oligoData.index[l/2:]]
#            form=form[l/2:]
        oligos_h=np.unique(data['hash'])
        oligoSummary=[]
        
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
#            v_sd=np.std(v)
#            leng_sd=np.std(leng)
            
            #form_left is time-weighted prob of left end to be ATP
            form_left_D=np.sum(tempDic['form_leftEnd'][:-1]*timeDiff)/t_total
            form_right_D=np.sum(tempDic['form_rightEnd'][:-1]*timeDiff)/t_total
            form_left_D_2=np.sum(tempDic['form_leftEnd2'][:-1]*timeDiff)/t_total
            form_right_D_2=np.sum(tempDic['form_rightEnd2'][:-1]*timeDiff)/t_total
            form_D = np.sum((tempDic['form_D']/leng)[:-1]*timeDiff)/t_total
            

            # form prob calculation, for relative distance from oligo's right end (barbed)
            form_time=np.zeros(m_intSpace)
            form_oligo=[]
            for i in tempDic.index[:-1]:
                temp=len(form[i])
                form_oligo.append(form[i]+[0]*(m_intSpace-temp))
            form_time=np.sum((np.array(form_oligo).T*timeDiff),axis=1) #amount of time spend in ADP for each position from right end. Divide by cumsum(length time) to get probability
            form_dic=dict(zip(['form_'+str(i+1) for i in xrange(m_intSpace)],form_time))

            # calculate length distribution
            l_dic=dict(zip(['l_'+str(i+1) for i in xrange(m_intSpace)],[0]*m_intSpace)) # amount of time spent in each length. Divide by total time to get probability
            for ln,t in zip(leng,timeDiff):
                l_dic['l_'+str(ln)]+=t
            
            dic={'v':v_mean,'l':leng_mean,'t_start':t_start,'t_end':t_end,'t_total':t_total,'form_left_D':form_left_D,'form_right_D':form_right_D,'form_D':form_D,'form_left_D_2':form_left_D_2,'form_right_D_2':form_right_D_2}
            dic.update(form_dic)
            dic.update(l_dic)
            oligoSummary.append(dic)
            
        oligoSummary=pd.DataFrame(oligoSummary)
        return oligoSummary
    
    @staticmethod
    def averageGraph(fN):
        # plot distribution graphs of average velocity, length, and lifetime
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
    def monoGraph(fN):
        # plot monomer data over time.
        monoData=pd.DataFrame.from_csv(fN+'.csv')
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
        # calculate statistics of types of event
        if isinstance(trace,str):
            trace=pd.DataFrame.from_csv(trace+'.csv')
        total_num=float(len(trace))
        diff_ratio=sum(trace['event_type']==intSpaceSim.event_type_func('diff'))
        end_ratio=sum(trace['event_type']==intSpaceSim.event_type_func('end'))
        hydro_ratio=total_num-diff_ratio-end_ratio
        print('diff:'+str(diff_ratio)+'/'+str(total_num)+', end:'+str(end_ratio)+'/'+str(total_num)+', hydro:'+str(hydro_ratio)+'/'+str(total_num))
        
    @staticmethod
    def treadmillingFilaments(oligoSummary):
        # an attempt to focus only on the filaments that should treadmill
        # not useful in general
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
        
     

class createSim(object):
    # a class for really creating a simulation. Also callable to execute the simulation.
#(fN,NN=10000,xi=-1.,rateHydro=1e-3,N=50,m=30,eps=-100.):
    def __init__(self,params={'NN':10000,'xi':-1.,'rateHydro':1e-3,'N':50,'m':30,'eps':-100.}):
        self.params=params
        
    def __call__(self,fN):
        NN,xi,rateHydro,N,m,eps=[self.params[x] for x in ['NN','xi','rateHydro','N','m','eps']]
        k=2
        h=2
        m_array=np.array([[0,m],[0,0]])
        eps=eps # to make it essentially impossible
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
    
        # RANDOM INITIALIZATION    
        position=np.random.choice(range(N),m,replace=False)
        proteinOrientation=np.repeat(np.reshape(np.ravel(np.indices((k,k)),1),[k**2,2]),params['m_array'].flatten(),axis=0)
        # [1,0]
        # for orientation, barbed end is 0, so barbed end is on the positive side. So velocity should be positive.
        proteinForm=np.ones(params['N'])
        mask=np.ones(params['N'])
        mask[position]=0
        proteinForm[np.where(mask)]=-1 #-1 indicate no protein, 0 is ATP, 1 is ADP
        proteinForm_pos=proteinForm[position]
    
        #initilize the space
        init={}
        init['pos']=position
        init['domain']=proteinOrientation
        init['form']=proteinForm_pos
        init['params']=params
        init['mode']='list'
        aa=InteractionSpace(init)
        
        #initialize the simulation object
        sim_params={'nstep':NN}
        bb=intSpaceSim(aa,params=sim_params)
    
        # simulate, then calculate data and save them.
        oligoData,monoData,trace,form=bb.sim()
        oligoData=bb.posCalc(oligoData)
        intSpaceSim.save(fN+'_oligo',oligoData)
        intSpaceSim.save(fN+'_mono',monoData)
        intSpaceSim.save(fN+'_trace',trace)
        oligoSummary=intSpaceSim.averageCalc(oligoData,form,bb.intSpace.params['N'],bb.intSpace.params['m'])
        intSpaceSim.save(fN+'_summary',oligoSummary)
        oneVal=intSpaceSim.oneValCalc(oligoSummary,trace,bb.intSpace.params['m'])
        intSpaceSim.save(fN+'_oneVal',oneVal)
        return fN


def createSim_star(ar):
    # a wrapper for parallel computation
    np.random.seed()
    return createSim(ar[1])(ar[0])

def repeatSim(params_list,nRep,maxNCore=24):
    # using parallel computation to simulation several instances
    # nRep is number of repeats for each parametre
    # params_list is list of parameters to simulate, for each parameter nRep instances will be simulated.
    args_list=params_list*nRep
    ll=len(params_list)
    fN_list=[x+'_REP-'+str(i/ll) for x,i in zip([fN_func(param) for param in params_list]*nRep,xrange(nRep*len(params_list)))]
    args_list=[(fN,x) for fN,x in zip(fN_list,args_list)]
    pool=Pool(min(len(params_list)*nRep,maxNCore))
    pool.map(createSim_star,args_list)
    
    # for debug purpose: run things non-parallel.
#    [createSim_star(args_list[i]) for i in xrange(nRep)]
    return args_list



def resultPlot(varName,varScale,params_list,nRep):
    # plot results of a repeatSim
    x_list=[param[varName] for param in params_list]
    if varName=='xi':
        i=0
        for x in x_list:
            x_list[i]=-x
            i+=1
    fN_list=[fN_func(param) for param in params_list]
    meanList=[]
    stdList=[]
    for fN in fN_list:
        repRslt=[]
        for i in xrange(nRep):
            temp=fN+'_REP-'+str(i)
            repRslt.append(pd.DataFrame.from_csv(temp+'_oneVal'+'.csv').loc[0])
        repRslt=pd.DataFrame(repRslt)
        
        meanDic=dict.fromkeys(repRslt.keys())
        stdDic=dict.fromkeys(repRslt.keys())
        for key in repRslt.keys():
            meanDic[key]=np.mean(repRslt[key])
            stdDic[key]=np.std(repRslt[key])
        meanList.append(meanDic)
        stdList.append(stdDic)
    
    meanList=pd.DataFrame(meanList)
    stdList=pd.DataFrame(stdList)
    
    labels=['v(lt)','v(t)','atp','l_atp','r_atp','l2_atp','r2_atp']
    for label in labels:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(x_list,meanList[label],'-b',lw=2,label=label+'_mean')
        ax.fill_between(x_list,meanList[label]-stdList[label],meanList[label]+stdList[label],color='purple',alpha=0.3,label=label+'_std')
        ax.set_xlabel(varName)
        ax.set_ylabel(label)
        ax.set_xscale(varScale)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc=3, bbox_to_anchor=(0.,1.02,1.,.102),ncol=2,mode='expand',borderaxespad=0.)
        plt.savefig(fN+'_trrun_'+label+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()
    
    for fN in fN_list:
        l_mean,l_std=distrPlot(fN,nRep=nRep,keyStub='l',m=params_list[0]['m'])
        f_mean,f_std=distrPlot(fN,nRep=nRep,keyStub='form',m=params_list[0]['m'])
    
    pickle.dump({'mean':meanList,'std':stdList,'var':x_list,'f_mean':f_mean,'f_std':f_std,'l_mean':l_mean,'l_std':l_std},open(varName+'_test.p','wr'))
    
 
def distrPlot(fN,nRep=30,keyStub='l',m=30):
    # plot the distribution of length or form for each subunit location
    keys=[keyStub+'_'+str(i+1) for i in xrange(m)]
    ys=[]
    for i in xrange(nRep):
        temp=fN+'_REP-'+str(i)
        su=pd.DataFrame.from_csv(temp+'_oneVal.csv').loc[0]
        y=[su[key] for key in keys]
        ys.append(y)
    ys=np.array(ys)
    y_mean=np.mean(ys,axis=0)
    y_std=np.std(ys,axis=0)
    x=range(m)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(x,y_mean,'-b',lw=2,label=keyStub+'--'+fN)
    ax.fill_between(x,y_mean-y_std,y_mean+y_std,color='purple',alpha=0.3,label='std')
    if keyStub=='l':
        ax.set_xlabel('length')
    elif keyStub=='form':
        ax.set_xlabel('pos from + end')
    else:
        raise ValueError('keyStub not valid')
    ax.set_ylabel('prob')
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc=3, bbox_to_anchor=(0.,1.02,1.,.102),ncol=1,mode='expand',borderaxespad=0.)
    plt.savefig(fN+'_DISTR-'+keyStub+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')   
    plt.close()
    return y_mean,y_std


def fN_func(param):
    # generate file name from parameters.
    NN,xi,rateHydro,N,m,eps=[param[x] for x in ['NN','xi','rateHydro','N','m','eps']]
    fN=dataFolder+'N-'+str(N)+'_m-'+str(m)+'_NN-'+str(NN)+'_xi-'+str(xi)+'_rh-'+str(rateHydro)+'_eps-'+str(eps)
    return fN          

def defaultParam(changeVar={}):
    # degenrate a parameter, with certains variables changed if indicated in changeVar.
    dft={'NN':1000000,'xi':-1.,'rateHydro':1e-1,'N':50,'m':30,'eps':-100.}
    for key,val in changeVar.iteritems():
        dft[key]=val
    return dft

if __name__ == "__main__":
    global dataFolder
#    dataFolder='/data/fangzhou/actin'
#    nRep=10
#    maxNCore=4
#    rH_list=np.logspace(-3,2,10)
#    params_list=[defaultParam({'rateHydro':rH}) for rH in rH_list]
#    args_list=repeatSim(params_list,nRep,maxNCore=maxNCore)
#    resultPlot('rateHydro','log',params_list,nRep)
#    
#    xi_list=-np.logspace(-3,1,10)
#    params_list=[defaultParam({'xi':xi}) for xi in xi_list]
#    args_list=repeatSim(params_list,nRep,maxNCore=maxNCore)
#    resultPlot('xi','log',params_list,nRep)
    
    parser=argparse.ArgumentParser(description='Simulation of actin treadmilling in 1D space.')
#    parser.add_argument('-m','--mode',type=str,nargs='?',default='var',
#                        help='What mode you would like to run. multiple is varying parameters: supply []; single is running one instance, supply []')
    parser.add_argument('-r','--rep',type=int,nargs='?',default=10,
                        help='number of repeats for var case')
    parser.add_argument('-c','--core',type=int,nargs='?', default=24,
                        help='number of cores to run in parallel')
    parser.add_argument('-o','--output',type=str,nargs='?', default='',
                        help='Which folder you would like the results and graphs to be outputed to.')
    parser.add_argument('-N','--nstep',type=int,nargs='?', default=100000,
                        help='number of steps for each simulation to run.')
    parser.add_argument('-v','--var',type=str,nargs='?',default='rateHydro',
                        help='the variable to vary. rateHydro, xi, N, m, eps')
    parser.add_argument('-l','--lim',type=int,nargs=2,default=[-3,1],
                        help='min and max of the variable value')
    parser.add_argument('-s','--scale',type=str,nargs='?',default='log',
                        help='scale for the variable values, log or linear')
    parser.add_argument('-n','--npts',type=int,nargs='?',default=10,
                        help='number of pts for the variable values between min and max')
    
    args=parser.parse_args()    
    
    
    nRep=args.rep
    maxNCore=args.core
    dataFolder=args.output
    NN=args.nstep
    varName=args.var
    varLim=args.lim
    varScale=args.scale
    varNPts=args.npts
    
    if varScale=='log':
        var_list=np.logspace(varLim[0],varLim[1],varNPts)
    elif varScale=='linear':
        var_list=np.linearspace(varLim[0],varLim[1],varNPts)
    else:
        raise ValueError('invalid variable scale option')
    
    params_list=[defaultParam({varName:var,'NN':NN}) for var in var_list]
    args_list=repeatSim(params_list,nRep,maxNCore=maxNCore)
    resultPlot(varName,varScale,params_list,nRep)
    
    
    
#    aa=defaultParam()
#    fN=fN_func(aa)
#    createSim(aa)(fN)
#    intSpaceSim.oligoGraph(aa['N'],fN,fN+'_oligo',graphParams={'truncate':True,'cutoffTime':20.,'cutoffNStep':aa['NN']/200,'mode':'centered'})

#    createSim(defaultParam())(fN_func(defaultParam()))
#    intSpaceSim.runningCalc(fN,fN+'_summary')
#    intSpaceSim.traceRunningPlot(fN,fN+'_trace')
#    
#    intSpaceSim.averageGraph(fN+'_summary')
#    intSpaceSim.monoGraph(fN+'_mono')
#    intSpaceSim.statistics(fN+'_trace')