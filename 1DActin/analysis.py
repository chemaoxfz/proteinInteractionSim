# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 02:36:08 2015

@author: xfz
"""
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import pdb
import matplotlib.animation as animation
from OneD_New import InteractionSpace
from OneD_New import H_actin

def actinAnalysis(fN):
    stats=pickle.load(open(fN,'rb'))
    energies=stats['energy']
    nRepeat=len(stats['energy'][0])
    nParam=len(stats['stepSize'])

    nRepeat=len(stats['energy'][0])
    nParam=len(stats['HParamStack'])
    filler = np.frompyfunc(lambda x: list(), 1, 1)
    pos=np.empty([nParam,nRepeat],dtype=np.object)
    filler(pos,pos)
    time=np.empty([nParam,nRepeat],dtype=np.object)
    filler(time,time)
    lEndForm=np.empty([nParam,nRepeat],dtype=np.object)
    filler(lEndForm,lEndForm)
    rEndForm=np.empty([nParam,nRepeat],dtype=np.object)
    filler(rEndForm,rEndForm)
    monoNum=np.empty([nParam,nRepeat],dtype=np.object)
    filler(monoNum,monoNum)
    
    
    monoMean=np.zeros([nParam,nRepeat])
    
    lEndMean=np.zeros([nParam,nRepeat])
    rEndMean=np.zeros([nParam,nRepeat])
    
    defaultPos=-stats['initParams']['N']/5    
    
    for paramIdx in xrange(nParam):    
        for repeat in xrange(nRepeat):
            nPts=len(stats['energy'][paramIdx][repeat])
            for ptIdx in xrange(nPts):
                state=stats['states'][paramIdx][repeat][ptIdx]
                fila=filament(stats['states'][paramIdx][repeat][ptIdx],stats['initParams'],defaultPos)
                filapos=filamentPos(fila,stats['initParams']['N'])
                pos[paramIdx][repeat].append(filapos)
#                pdb.set_trace()
                lEndForm[paramIdx][repeat].append([state[2][np.where(state[0]==x[0])] for x in fila])
                rEndForm[paramIdx][repeat].append([state[2][np.where(state[0]==x[-1])] for x in fila])
                monoNum[paramIdx][repeat].append(freeMonoNum(state,fila,stats['initParams']))
                time[paramIdx][repeat].append(stats['obsStart']+stats['stepSize'][paramIdx]*ptIdx)
                
                # Detect whether monomer ADP exists
#                for t in xrange(20):
#                    if state[2][t]!=0 and (state[0][t]+1)%30 not in state[0] and (state[0][t]-1)%30 not in state[0]:
#                            print('shooooooooooooooooooooooooot!')
#                            print(state)
#                            pdb.set_trace()
#                            print(t)
#                            print(ptIdx)
#                            aa=InteractionSpace(stats['initParams'],lambda a,b:H_actin(a,b,1,1,1,1,1))
#                            print(aa.checkATP(state))
#                            pdb.set_trace()
                            
#    pdb.set_trace()
#        lEndMean[paramIdx][repeat]=np.mean((lEndForm[paramIdx][repeat]))
#        rEndMean[paramIdx][repeat]=np.mean((rEndForm[paramIdx][repeat]))
#        monoMend[paramIdx][repeat]=np.mean((monoNum[paramIdx][repeat]))
    pickle.dump({'pos':pos,'time':time,'lEndForm':lEndForm,'rEndForm':rEndForm,'monoNum':monoNum},open(fN+'_anal.p','wr'))

    fig, ax = plt.subplots(1)
    colors=plt.cm.coolwarm(np.linspace(0,1,nParam))
    for figIdx in xrange(nParam):
        maxLen=max(map(len,pos[0][0]))
        if maxLen==1:
            ax.plot(time[0][0], pos[0][0],color=[figIdx], lw=0.5,  label=str(stats['HParamStack'][figIdx]))
        else: 
            pos=map(lambda x: padPos(x,maxLen,defaultPos),pos[0][0])
            ax.plot(time[0][0], pos,color=[figIdx], lw=0.5,  label=str(stats['HParamStack'][figIdx]))
        ax.set_title('position over time,(T,eps,xi)= (%.4f' %stats['HParamStack'][figIdx][0]+',%.4f)'%stats['HParamStack'][figIdx][1]+',%.4f'%stats['HParamStack'][figIdx][2])
        
    #    ax.legend(loc='upper left')
        ax.set_xlabel('time')
        ax.set_ylabel('pos')
        plt.savefig(fN+str(figIdx)+'.pdf')
        plt.show()
    
#    for figIdx in xrange(nParam):
#        fig, ax = plt.subplots(1)
#        ax.plot(time[0][0], monoNum[0][0],'bo', lw=0.5,  label='monoNum')
#        ax.set_title('monoNum over time,(xi,T)= (%.4f' %stats['xiTStack'][figIdx][0]+',%.4f)'%stats['xiTStack'][figIdx][1])
#        
#    #    ax.legend(loc='upper left')
#        ax.set_xlabel('time')
#        ax.set_ylabel('monoNum')
#        plt.savefig(fN+str(figIdx)+'_mono.pdf')
#        plt.show()
        


def plotFromFile(fileName):  
    stats=pickle.load(open(fileName,'rb'))
    energies=stats['energy']
    nRepeat=len(stats['energy'][0])
    nParam=len(stats['stepSize'])
    
    filler = np.frompyfunc(lambda x: list(), 1, 1)
    pos=np.empty([nParam,nRepeat],dtype=np.object)
    filler(pos,pos)
    time=np.empty([nParam,nRepeat],dtype=np.object)
    time=filler(time,time)

    
    defaultPos=-stats['initParams']['N']/5    
    
    for paramIdx in xrange(nParam):    
        for repeat in xrange(nRepeat):
            nPts=len(stats['energy'][paramIdx][repeat])
            for ptIdx in xrange(nPts):
#                pdb.set_trace()
                fila=filament(stats['states'][paramIdx][repeat][ptIdx],stats['initParams'],defaultPos)
                filapos=filamentPos(fila,stats['initParams']['N'])
#                if paramIdx==0 and ptIdx>22: pdb.set_trace()
                pos[paramIdx][repeat].append(filapos)
                time[paramIdx][repeat].append(stats['obsStart']+stats['stepSize'][paramIdx]*ptIdx)
#    time=[]
#    for idx in xrange(len(stats['stepSize'])):
#        time.append(stats['obsStart']+np.arange(np.floor(stats['obsDur']/stats['stepSize'][idx])*stats['stepSize'][idx]))
#    pdb.set_trace()
#    data={'pos':pos,'time':time}
#    pickle.dump(data,open('asdf','wr'))

    

def padPos(pos,maxLen,defaultPos):
    return np.pad(pos,(0,maxLen-len(pos)),'constant',constant_values=(defaultPos))

def filament(state,params,defaultPos):
    segments=consecutive(state,params)
    fila=[seg for seg in segments if len(seg) > params['m']/2.]
    if not fila:
        fila=[np.array([defaultPos])]
    return fila

def filamentPos(fila,N):
    return [filaMedian(fil)+1e-5 for fil in fila] # adding 1e-5 to make sure it's not zero.

def filaMedian(fil):
    l=len(fil)
    if l%2==0:
        return (fil[l/2]+fil[l/2-1])/2.
    else: return fil[(l-1)/2]

def consecutive(state,params):
    pos=np.sort(state[0])
    segments=np.array(np.split(pos, np.where(np.diff(pos) != 1)[0]+1))
    if (segments[0][0]-segments[-1][-1])%params['N']==1:
        ends=np.concatenate((segments[-1],segments[0])) #[0,1,2] and [N-2,N-1] becomes [-2,-1,0,1,2]
        segments=np.array([ends,segments[1:-1]])
    return segments

def freeMonoNum(state,fila,params):
#    connectedSlots=np.union1d(np.union1d(state[0],(state[0]+1)%params['N']),(state[0]-1)%params['N'])
    aa=[]
    occupied=[np.union1d(aa,x) for x in fila]
    N=params['N']
#    pdb.set_trace()
    mono=[x for x in state[0] if x not in occupied[0]]
    monoNum=len(mono)
    return monoNum

def formPercentage(state,params):
    pass

def update_plot(i):
    global states
#    pdb.set_trace()
    x = np.tile(states[0][i][0],2)+np.concatenate((np.zeros(m),np.ones(m)*deviation))
    x = np.ravel(np.reshape(x,[2,m]),1)
    y = np.zeros(2*m)
    z = np.vstack((x,y)).transpose()
    points.set_offsets(z)
    colors = np.lib.pad(np.repeat(np.reshape(states[0][i][1].flatten()*1/float(k-1),[2*m,1]),2,axis=1),(0,1),'constant',constant_values=(0))
    points.set_facecolor(colors)
    
    xcenters=states[0][i][0]+deviation/2
    ycenters=np.zeros(m)
    zcenters = np.vstack((xcenters,ycenters)).transpose()
    centers.set_offsets(zcenters)
    colorscenters = np.zeros([3,m])
    colorscenters[2]=states[0][i][2]*1/float(h-1)
    colorscenters[0]=1-colorscenters[2]
    colorscenters=colorscenters.T
    centers.set_facecolor(colorscenters)        
    
    aa=stepSize*i
    
    time_text.set_text('time=%.1f' % aa)
    energy_text.set_text('energy=%.1f' % energy[0][i])
    
    defaultPos=-20
    fila=filament(states[0][i],params,defaultPos)
    filapos=filamentPos(fila,params['N'])
    pos_text.set_text('pos=%.1f'% filapos[0])
    return points,time_text,energy_text

def TE_animate(fN,elt):
    global states
    global deviation
    global points
    global centers
    global time_text
    global energy_text
    global pos_text
    global m
    global k
    global h
    global stepSize
    global energy
    global params
    simDict=pickle.load(open(fN,'r'))
    states=simDict['states'][elt]
    energy=simDict['energy'][elt]
    m=simDict['initParams']['m']
    k=simDict['initParams']['k']
    N=simDict['initParams']['N']
    h=simDict['initParams']['h']
    params=simDict['initParams']
    stepSize=simDict['stepSize'][0]
    
    interval=500
    frames=len(states[0])
    
    fig=plt.figure(figsize=(27,2))
    ax=fig.add_subplot(111,autoscale_on=False,xlim=(-1,N),ylim=(-N/5.,N/5.))
    deviation=0.3 #deviation is the distance between the center of right domain from left domain of the same protein
#    pdb.set_trace()
    x = np.tile(states[0][0][0],2)+np.concatenate((np.zeros(m),np.ones(m)*deviation))
    x = np.ravel(np.reshape(x,[2,m]),1)
    y = np.zeros(2*m)
    z = np.vstack((x,y)).transpose()
    #currently 1 is yellow, 0 is black. so barbed is yellow. pointed is black.
    colors = np.lib.pad(np.repeat(np.reshape(states[0][0][1].flatten()*1/float(k-1),[2*m,1]),2,axis=1),(0,1),'constant',constant_values=(0))
    points=ax.scatter(x,y,s=np.pi*30,c=colors)   
    
    xcenters=states[0][0][0]+deviation/2
    ycenters=np.zeros(m)
    colorscenters = np.zeros([3,m])
    colorscenters[2]=states[0][0][2]*1/float(h-1)
    colorscenters[0]=1-colorscenters[2]
    colorscenters=colorscenters.T
    centers=ax.scatter(xcenters,ycenters,s=np.pi*30,c=colorscenters,marker=(5,0))    
    
    time_text=ax.text(0.02,0.95,'',transform=ax.transAxes)
    energy_text=ax.text(0.02,0.90,'',transform=ax.transAxes)
    pos_text=ax.text(0.02,0.85,'',transform=ax.transAxes)
    ani=animation.FuncAnimation(fig,update_plot,np.arange(1,frames),
                                interval=interval, blit=False)
#    ani.save(fN+'.gif',writer='imagemagick',fps=8)
    plt.show(ax)

if __name__ == "__main__":
#    fN='./data/pleaseSuceed.p3'
    fN='./noAlt'
#    plotFromFile(fN)
    actinAnalysis(fN)
#    elt=0
#    TE_animate(fN,elt)