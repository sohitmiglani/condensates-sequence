#!/usr/bin/env python
# coding: utf-8

# analyze structure of the dense phase in a multicanonical simulation, where we compare different seqs at same tau

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.optimize
from matplotlib import cm


# In[2]:


#set up param space
blockSizeList=[1,2,3,4,6,12]
seqList=['L24_b'+str(i) for i in blockSizeList]
linearD=30
j=0.05
length=24


saveBool=True


tauList=[-0.0025,-0.005,-0.008,-0.011,-0.013,-0.015,-0.019]
#step 2: -0.011
#step3: -0.013
#step 4: -0.015
#5: -0.017
#6: -0.019
#7: -0.021
#8: -0.023

path='densePhase/'
tcPath='criticalPoint/'

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']


# bin data by density and average

# In[3]:


binWidth=0.02
binEdges=np.arange(0,1,binWidth)
binMiddles=binEdges[1:]-binWidth/2

def binData(densitySeries,dataSeries):
    
    densitySeries=densitySeries.to_numpy()
    dataSeries=dataSeries.to_numpy()

    binMeans=np.full(len(binEdges)-1,-1,dtype=float)
    binError=np.full(len(binEdges)-1,-1,dtype=float)
    
    binIndices=np.digitize(densitySeries,binEdges)
    
    for i in range(len(binEdges)-1):
        
        data_thisBin=dataSeries[binIndices==i]
        
        if len(data_thisBin)>0:
            binMeans[i]=np.mean(data_thisBin)
            binError[i]=np.std(data_thisBin)
    
    #remove empty bins
    x=binMiddles[binMeans>-1]
    binMeans=binMeans[binMeans>-1]
    binError=binError[binError>-1]
    
    return x,binMeans,binError


# In[4]:


def getWeightingData(seq,beta1):
    
    #determine filename
    searchSeq=seq+'_'
    allFiles=os.listdir(path)
    for thisFile in allFiles:
        preweightingCondition=thisFile.find('preweighting')!=-1
        seqCondition=thisFile.find(searchSeq)!=-1
        if preweightingCondition and seqCondition:
            fileBeta1=thisFile[thisFile.find('beta1')+6:thisFile.find('.txt')]
            fileBeta1=float(fileBeta1)
            roundedFileBeta1=round(fileBeta1,6)
            if abs(roundedFileBeta1-beta1)<0.0001:
                weightingFilename=thisFile
    
    
    #load data
    multicanonicalWeightingArray=np.loadtxt(path+weightingFilename,skiprows=9,delimiter=' ')
    multicanonicalWeighting={}
    for i in range(multicanonicalWeightingArray.shape[0]):
        multicanonicalWeighting[int(multicanonicalWeightingArray[i,0])]=multicanonicalWeightingArray[i,1]
        
    return multicanonicalWeighting


# In[5]:


#load data, reweight, find phi1 and phi2
def getPhaseBoundary(data,preweighting):
    
    #get weighting as new column
    
    dataCopy=data.copy()
    
    dataCopy['N']=dataCopy['density']*(linearD**3)/length
    dataCopy['N']=dataCopy['N'].astype('int32')
    weightArray=np.array([preweighting[i] for i in dataCopy['N']])
    weights=np.log(weightArray)
    
    dataCopy['w']=np.exp(weights)

    
#     plt.hist(dataCopy['density'],weights=dataCopy['w'],density=True,bins=np.arange(0,1,0.02))
    
    #get means
        
    binEdges=np.arange(np.min(dataCopy['N'])-0.5,np.max(dataCopy['N'])+0.5,1)
    hist,_ = np.histogram(dataCopy['N'],weights=dataCopy['w'], density=True,bins=binEdges)
    
    prelimThreshold=400
    
    #find modes
    countVector=binEdges[1:]-0.5
    thresholdIndex=np.argwhere(abs(countVector-prelimThreshold)<10**-5)[0,0]

    
    #in histogram bin index
    mode1Index=np.argmax(hist[0:thresholdIndex])
    mode2Index=np.argmax(hist[thresholdIndex:])+thresholdIndex
    

    realThresholdIndex=int(round((mode1Index+mode2Index)/2))
    
    diluteData=dataCopy[dataCopy['N']<=countVector[realThresholdIndex]]
    denseData=dataCopy[dataCopy['N']>countVector[realThresholdIndex]]
    
    diluteMean=np.average(diluteData['N'],weights=diluteData['w'])
    denseMean=np.average(denseData['N'],weights=denseData['w'])
    
    phi1=diluteMean*length/(linearD**3)
    phi2=denseMean*length/(linearD**3)

    
    return phi1,phi2
    


# now plot only at phi1, phi2, but as a function of tau

# In[6]:


#take points at +/- phiRadius from phi1, phi2
phiRadius=0.01

#return mean and std of points within phiRadius of testPhi
def getCoexistenceData(densitySeries,dataSeries,testPhi):
    
    densitySeries=densitySeries.to_numpy()
    dataSeries=dataSeries.to_numpy()
    
    dataSeries=dataSeries[abs(densitySeries-testPhi)<=phiRadius]
    
    
    
    
    return np.mean(dataSeries),np.std(dataSeries)


# In[7]:


#make dictionary of critical points for tau calculation, and array for plotting
betaCDict={}

#for each sequence, mean betaC, std betaC, mean phiC, std phiC
criticalData=np.zeros([len(seqList),4])

for counter,seq in enumerate(seqList):
        
    tcFilename='tcFit_'+seq+'_j0.05_rep'
    betaList=[]
    phiList=[]
    for rep in range(3):
        thisData=np.load(tcPath+tcFilename+str(rep)+'.npy')
        betaList.append(thisData[0])
        phiList.append(thisData[2])

    betaC=np.mean(betaList)
    phiC=np.mean(phiList)
    phiCStd=np.std(phiList)
        
    tList=[1/i for i in betaList]
    tMean=np.mean(tList)
    tStd=np.std(tList)
    
    criticalData[counter]=np.array([tMean,tStd,phiC,phiCStd])
    
    betaCDict[seq]=betaC


# In[8]:


#get phi2 values so you can normalize color map
phi2List=[]
for seq in seqList:
    
    
    betaC=betaCDict[seq]
    

    for counter,thisTau in enumerate(tauList):    
      
    
        beta1=betaC/(thisTau+1)
        beta1=round(beta1,6)

        
        filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
        thisData=pd.read_csv(path+filename)
        
        #weighting data
        weighting=getWeightingData(seq,beta1)
        
        #get phi1 and phi2
        phi1,phi2=getPhaseBoundary(thisData,weighting)
        phi2List.append(phi2)

phi2List.sort()
norm = cm.colors.Normalize(vmax=phi2List[-1], vmin=phi2List[0])
norm_phi2 = cm.colors.Normalize(vmax=phi2List[-1], vmin=phi2List[0])


# In[10]:


#set overall figure properties and make legends

markerset=['.','_','^','s','*','o']
cmap = plt.cm.get_cmap('viridis_r')
markersizes=[18,14,12,10,16,12]

figureSizes=[8,6]
labelSize=24
tickSize=18


#some dummy plots for legends
plt.figure()
for i in range(len(seqList)):
    if seqList[i]=='L24_b12':
        plt.plot(1,1,marker=markerset[i],linewidth=0,color='k',markerfacecolor='none',markersize=markersizes[i],label=blockSizeList[i])
    elif seqList[i]=='L24_b2':
        plt.plot(1,1,marker=markerset[i],linewidth=0,color='k',mew=3,markerfacecolor='none',markersize=markersizes[i],label=blockSizeList[i])

    else:
        plt.plot(1,1,marker=markerset[i],linewidth=0,color='k',markersize=markersizes[i],label=blockSizeList[i])
plt.legend(title='Block size',title_fontsize=15,fontsize=15)

if saveBool:
    plt.savefig('fig4_domainLabels.svg')


plt.figure()
Z,_=np.meshgrid(phi2List, phi2List)
plt.imshow(Z, cmap=cmap, norm=norm)
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label('$\\phi_{\mathrm{dense}}$',fontsize=labelSize)
cbar.ax.tick_params(labelsize=14)

if saveBool:
    plt.savefig('fig4_phi2_colorbar.svg')


# In[10]:



#for properties measured directly in simulation
def plotDensePhaseProperty(statisticName,yAxisLabel):

    

    plt.figure(figsize=(figureSizes[0],figureSizes[1]))


    for counter1,seq in enumerate(seqList):

        thisMarker=markerset[counter1]
        thisMarkerSize=markersizes[counter1]

        betaC=betaCDict[seq]

        #tau, mean of means, std of means
        phi1Points=np.zeros([len(tauList),3])
        phi2Points=np.zeros([len(tauList),3])


        for counter,thisTau in enumerate(tauList):    


            beta1=betaC/(thisTau+1)
            beta1=round(beta1,6)


            filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
            thisData=pd.read_csv(path+filename)

            #weighting data
            weighting=getWeightingData(seq,beta1)

            #get phi1 and phi2
            phi1,phi2=getPhaseBoundary(thisData,weighting)


            thisColor=cmap(norm(phi2))

            #get data around phi1 and phi2
            phi2Data=getCoexistenceData(thisData['density'],thisData[statisticName],phi2)

            phi2Points[counter]=np.array([abs(thisTau),phi2Data[0],phi2Data[1]])



            if seq=='L24_b12':
                plt.errorbar(abs(thisTau),phi2Data[0],yerr=phi2Data[1],marker=thisMarker,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)
            
            elif seq=='L24_b2':
                plt.errorbar(abs(thisTau),phi2Data[0],yerr=phi2Data[1],marker=thisMarker,mew=3,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)

            else:
                plt.errorbar(abs(thisTau),phi2Data[0],yerr=phi2Data[1],marker=thisMarker,color=thisColor,markerSize=thisMarkerSize)

    plt.xlabel('Scaled temperature $|T-T_\mathrm{c}|$ / $T_\mathrm{c}$',fontsize=labelSize)
    plt.ylabel(yAxisLabel,fontsize=labelSize)
    plt.tick_params(axis='both',labelsize=tickSize)
    
    if statisticName=='s_mean':
        plt.margins(y=0.5)
        x0, x1, y0, y1 = plt.axis()
        plt.ylim(ymin=y0+(y1-y0)*0.2-0.3)

    
    plt.xlim(xmin=0)

    plt.xticks(ticks=np.arange(0,0.02025,0.0025),labels=['0','','0.005','','0.01','','0.015','','0.02'])
    
    if statisticName=='rg_mean':
        plt.ylim(1.8,2.13)
    

    
    if saveBool:
        plt.savefig('fig4_phi2_'+str(statisticName)+'.svg', bbox_inches = "tight")


# In[12]:


plotDensePhaseProperty('s_mean','Self-bonds $s$')
plotDensePhaseProperty('t_mean','Trans-bonds $t$')
plotDensePhaseProperty('rg_mean','Radius of gyration $R_g$')


# In[13]:


#viscosity at phi2 with density color scheme


plt.figure(figsize=(figureSizes[0],figureSizes[1]))


for counter1,seq in enumerate(seqList):
    
    thisMarker=markerset[counter1]
    thisMarkerSize=markersizes[counter1]


    betaC=betaCDict[seq]
    


    for counter,thisTau in enumerate(tauList):    
      
    
        beta1=betaC/(thisTau+1)
        beta1=round(beta1,6)

        
        filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
        thisData=pd.read_csv(path+filename)
        
        #weighting data
        weighting=getWeightingData(seq,beta1)
        
        #get phi1 and phi2
        phi1,phi2=getPhaseBoundary(thisData,weighting)
        

        thisColor=cmap(norm(phi2))
        
        #get data around phi1 and phi2
        phi2Data=getCoexistenceData(thisData['density'],thisData['t_mean'],phi2)

        
        thisViscosity=np.exp(beta1)*(1/beta1)*(phi2/length)*phi2Data[0]**2
        thisError=2*np.exp(beta1)*(1/beta1)*(phi2/length)*phi2Data[1]*phi2Data[0]

        
        
        if seq=='L24_b12':
            plt.errorbar(abs(thisTau),thisViscosity,yerr=thisError,marker=thisMarker,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)
        
        elif seq=='L24_b2':
            plt.errorbar(abs(thisTau),thisViscosity,yerr=thisError,marker=thisMarker,mew=3,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)
            
        else:
            plt.errorbar(abs(thisTau),thisViscosity,yerr=thisError,marker=thisMarker,color=thisColor,markerSize=thisMarkerSize)


plt.ylabel('"Viscosity" $\\eta$ ($\\tau_0/m^3$)',fontsize=labelSize)
plt.xlabel('Scaled temperature $|T-T_\mathrm{c}|$ / $T_\mathrm{c}$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)

plt.xlim(xmin=0)
plt.ylim(ymax=9)
plt.xticks(ticks=np.arange(0,0.02025,0.0025),labels=['0','','0.005','','0.01','','0.015','','0.02'])

if saveBool:
    plt.savefig('fig4_phi2_viscosity.svg', bbox_inches = "tight")


# In[14]:


#diffusion at phi2 with density color scheme

plt.figure(figsize=(figureSizes[0]+1,figureSizes[1]))


for counter1,seq in enumerate(seqList):
    
    thisMarker=markerset[counter1]
    thisMarkerSize=markersizes[counter1]


    betaC=betaCDict[seq]
    


    for counter,thisTau in enumerate(tauList):    
      
    
        beta1=betaC/(thisTau+1)
        beta1=round(beta1,6)

        
        filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
        thisData=pd.read_csv(path+filename)
        
        #weighting data
        weighting=getWeightingData(seq,beta1)
        
        #get phi1 and phi2
        phi1,phi2=getPhaseBoundary(thisData,weighting)
        

        thisColor=cmap(norm(phi2))
        
        #get data around phi1 and phi2
        phi2Data=getCoexistenceData(thisData['density'],thisData['t_mean'],phi2)

        
        thisD=(1/np.exp(beta1))/phi2Data[0]
        thisError=(1/np.exp(beta1))*(1/(phi2Data[0]**2))*phi2Data[1]

        
        
        if seq=='L24_b12':
            plt.errorbar(abs(thisTau),thisD,yerr=thisError,marker=thisMarker,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)
        
        elif seq=='L24_b2':
            plt.errorbar(abs(thisTau),thisD,yerr=thisError,marker=thisMarker,mew=3,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)
            
            
        else:
            plt.errorbar(abs(thisTau),thisD,yerr=thisError,marker=thisMarker,color=thisColor,markerSize=thisMarkerSize)


plt.ylabel('Diffusivity $D$ ($m^2/\\tau_0$)',fontsize=labelSize)
plt.xlabel('Scaled temperature $|T-T_\mathrm{c}|$ / $T_\mathrm{c}$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)


plt.xlim(xmin=0)
plt.ylim(ymax=0.065)
plt.xticks(ticks=np.arange(0,0.02025,0.0025),labels=['0','','0.005','','0.01','','0.015','','0.02'])
plt.tight_layout()
#plt.savefig('fig3_phi2_diffusion.svg')

if saveBool:
    plt.savefig('fig4_phi2_diffusion.svg')


# In[15]:


#Rg at phi1


#get phi1 values so you can normalize color map
phi1List=[]
for seq in seqList:
    
    
    betaC=betaCDict[seq]
    

    for counter,thisTau in enumerate(tauList):    
      
    
        beta1=betaC/(thisTau+1)
        beta1=round(beta1,6)

        
        filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
        thisData=pd.read_csv(path+filename)
        
        #weighting data
        weighting=getWeightingData(seq,beta1)
        
        #get phi1 and phi2
        phi1,phi2=getPhaseBoundary(thisData,weighting)
        phi1List.append(phi1)

phi1List.sort()
norm_phi1 = cm.colors.Normalize(vmax=phi1List[-1], vmin=phi1List[0])


plt.figure(figsize=(figureSizes[0],figureSizes[1]))



for counter1,seq in enumerate(seqList):

    thisMarker=markerset[counter1]
    thisMarkerSize=markersizes[counter1]

    betaC=betaCDict[seq]

    #tau, mean of means, std of means
    phi1Points=np.zeros([len(tauList),3])
    phi2Points=np.zeros([len(tauList),3])


    for counter,thisTau in enumerate(tauList):    


        beta1=betaC/(thisTau+1)
        beta1=round(beta1,6)


        filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
        thisData=pd.read_csv(path+filename)

        #weighting data
        weighting=getWeightingData(seq,beta1)

        #get phi1 and phi2
        phi1,phi2=getPhaseBoundary(thisData,weighting)


        thisColor=cmap(norm_phi1(phi1))

        #get data around phi1 and phi2
        phi1Data=getCoexistenceData(thisData['density'],thisData['rg_mean'],phi1)

        phi1Points[counter]=np.array([abs(thisTau),phi1Data[0],phi1Data[1]])



        if seq=='L24_b12':
            plt.errorbar(abs(thisTau),phi1Data[0],yerr=phi1Data[1],marker=thisMarker,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)

        elif seq=='L24_b2':
            plt.errorbar(abs(thisTau),phi1Data[0],yerr=phi1Data[1],marker=thisMarker,mew=3,color=thisColor,markerfacecolor='none',markersize=thisMarkerSize)

        else:
            plt.errorbar(abs(thisTau),phi1Data[0],yerr=phi1Data[1],marker=thisMarker,color=thisColor,markerSize=thisMarkerSize)

plt.xlabel('Scaled temperature $|T-T_\mathrm{c}|$ /$ T_\mathrm{c}$',fontsize=labelSize)
plt.ylabel('Radius of gyration $R_g$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)


plt.xlim(xmin=0)
plt.xticks(ticks=np.arange(0,0.02025,0.0025),labels=['0','','0.005','','0.01','','0.015','','0.02'])

plt.margins(y=0.5)
x0, x1, y0, y1 = plt.axis()
plt.ylim(1.8,2.13)

plt.tight_layout()


# if saveBool:
#     plt.savefig('fig4_phi1_rg.svg')
    
plt.figure()
Z,_=np.meshgrid(phi2List, phi2List)
plt.imshow(Z, cmap=cmap, norm=norm_phi1)
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label('$\\phi_{\mathrm{dilute}}$',fontsize=labelSize)
cbar.ax.tick_params(labelsize=14)

# if saveBool:
#     plt.savefig('fig4_phi1_colorbar.svg')


# In[16]:


#plot Rg dense and dilute together
offset=0.0006
cmap2 = plt.cm.get_cmap('viridis_r')
cmap1 = plt.cm.get_cmap('cool')


#get phi1 values so you can normalize color map
phi1List=[]
for seq in seqList:
    
    
    betaC=betaCDict[seq]
    

    for counter,thisTau in enumerate(tauList):    
      
    
        beta1=betaC/(thisTau+1)
        beta1=round(beta1,6)

        
        filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
        thisData=pd.read_csv(path+filename)
        
        #weighting data
        weighting=getWeightingData(seq,beta1)
        
        #get phi1 and phi2
        phi1,phi2=getPhaseBoundary(thisData,weighting)
        phi1List.append(phi1)

phi1List.sort()
norm_phi1 = cm.colors.Normalize(vmax=phi1List[-1], vmin=phi1List[0])


plt.figure(figsize=(figureSizes[0],figureSizes[1]))



for counter1,seq in enumerate(seqList):

    thisMarker=markerset[counter1]
    thisMarkerSize=markersizes[counter1]

    betaC=betaCDict[seq]

    #tau, mean of means, std of means
    phi1Points=np.zeros([len(tauList),3])
    phi2Points=np.zeros([len(tauList),3])


    for counter,thisTau in enumerate(tauList):    


        beta1=betaC/(thisTau+1)
        beta1=round(beta1,6)


        filename='densePhaseData_'+seq+'_beta_'+str(beta1)+'_j'+str(j)+'.dat'
        thisData=pd.read_csv(path+filename)

        #weighting data
        weighting=getWeightingData(seq,beta1)

        #get phi1 and phi2
        phi1,phi2=getPhaseBoundary(thisData,weighting)


        thisColor_phi1=cmap1(norm_phi1(phi1))
        thisColor_phi2=cmap2(norm_phi2(phi2))


        #get data around phi1 and phi2
        phi1Data=getCoexistenceData(thisData['density'],thisData['rg_mean'],phi1)

        phi1Points[counter]=np.array([abs(thisTau),phi1Data[0],phi1Data[1]])

        phi2Data=getCoexistenceData(thisData['density'],thisData['rg_mean'],phi2)
        phi2Points[counter]=np.array([abs(thisTau),phi2Data[0],phi2Data[1]])


        if seq=='L24_b12':
            plt.errorbar(abs(thisTau),phi2Data[0],yerr=phi2Data[1],marker=thisMarker,color=thisColor_phi2,markerfacecolor='none',markersize=thisMarkerSize)
            plt.errorbar(abs(thisTau)+offset,phi1Data[0],yerr=phi1Data[1],marker=thisMarker,color=thisColor_phi1,markerfacecolor='none',markersize=thisMarkerSize)

            
        elif seq=='L24_b2':
            plt.errorbar(abs(thisTau),phi2Data[0],yerr=phi2Data[1],marker=thisMarker,mew=3,color=thisColor_phi2,markerfacecolor='none',markersize=thisMarkerSize)
            plt.errorbar(abs(thisTau)+offset,phi1Data[0],yerr=phi1Data[1],marker=thisMarker,mew=3,color=thisColor_phi1,markerfacecolor='none',markersize=thisMarkerSize)

        else:
            plt.errorbar(abs(thisTau),phi2Data[0],yerr=phi2Data[1],marker=thisMarker,color=thisColor_phi2,markerSize=thisMarkerSize)
            plt.errorbar(abs(thisTau)+offset,phi1Data[0],yerr=phi1Data[1],marker=thisMarker,color=thisColor_phi1,markerSize=thisMarkerSize)

            
plt.xlabel('Scaled temperature $|T-T_\mathrm{c}|$ /$ T_\mathrm{c}$',fontsize=labelSize)
plt.ylabel('Radius of gyration $R_g$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)


plt.xlim(xmin=0)
plt.xticks(ticks=np.arange(0,0.02025,0.0025),labels=['0','','0.005','','0.01','','0.015','','0.02'])

plt.margins(y=0.5)
x0, x1, y0, y1 = plt.axis()
plt.ylim(ymin=y0+(y1-y0)*0.2)


# plt.tight_layout()

# plt.savefig('fig3_rg_both.svg')


if saveBool:
    plt.savefig('fig4_rg_both.svg')
    
plt.figure()
Z,_=np.meshgrid(phi2List, phi2List)
plt.imshow(Z, cmap=cmap1, norm=norm_phi1)
cbar=plt.colorbar(orientation='horizontal')
cbar.set_label('$\\phi_{\mathrm{dilute}}$',fontsize=labelSize)
cbar.ax.tick_params(labelsize=14)

# plt.savefig('fig4_phi1_colorbar.svg')

# if saveBool:
#     plt.savefig('fig3_phi1_colorbar.svg')


# In[17]:
