#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import matplotlib.pyplot as plt


# In[9]:


def plotgE(seqList,labelList,s=False):
    
    plt.figure(figsize=(10,6))
    ax=plt.gca()
        
    for counter,seq in enumerate(seqList):

        filename='extractednE_'+seq+'.npy'

        data=np.load(path+filename)
        
        #number of s bonds
        if s:
            x=np.arange(0,data.shape[0],1)
            data=np.flip(data)
        #energy =-s*epsilon 
        else:
            x=np.arange(-data.shape[0]+1,1,1)
        

        plt.plot(x,data,label=labelList[counter])
    
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
        ax = plt.gca()
        ax.yaxis.get_offset_text().set_fontsize(18+5)
#     handles, labels = ax.get_legend_handles_labels()
#     ax.legend(handles[::-1], labels[::-1], title='Domain size',title_fontsize=15,fontsize=15,loc=9)
    
def plotgELog(seqList,labelList,s=False):
    
    plt.figure(figsize=(10,6))
    ax=plt.gca()
    
    for counter,seq in enumerate(seqList):

        filename='extractednE_'+seq+'.npy'

        data=np.load(path+filename)
        
        #number of s bonds
        if s:
            x=np.arange(0,data.shape[0],1)
            data=np.flip(data)
        #energy =-s*epsilon 
        else:
            x=np.arange(-data.shape[0]+1,1,1)

        plt.semilogy(x,data,label=labelList[counter])
        plt.xticks(ticks=[0,2,4,6,8,10,12])


# In[11]:


path='densityOfStatesData/'


blockSizes=[1,2,3,4,6,12]
seqList=['L24_block'+str(i) for i in blockSizes]
labelList=[str(i) for i in blockSizes]

saveBool=True

labelSize=24
plotgE(seqList,labelList,True)
plt.xlabel('Self-bonds $s$',fontsize=labelSize)
plt.ylabel('Density of states $g(s)$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=18)

if saveBool==True:
    plt.savefig('fig2_dos_linear.svg')

  
plotgELog(seqList,labelList,True)
plt.xlabel('Self-bonds $s$',fontsize=22)
plt.ylabel('$g(s)$ (log scale)',fontsize=22)
plt.tick_params(axis='both',labelsize=18)
# if saveBool==True:
#     plt.savefig('fig2_dos_log.svg')


# In[ ]:




