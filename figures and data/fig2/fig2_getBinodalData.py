#!/usr/bin/env python
# coding: utf-8

# load data and preweighting iteratively and return:
# 1) phi1 and phi2
# 2) dilute fraction (check that the fitting worked)
# 3) surface tension (test for separation)
# 4) mu*
# 5) generate picture of reweighted histogram as a check

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.optimize
from scipy import stats
from sklearn import mixture


# convention:
#     beta0=temp of the simulation from which the current weighting was extracted (not meaningful for uniform)
#     beta1=temp of the actual simulation being analyzed
#     betaNew=temp of the simulation we want to do next, after reweighting

# In[2]:


def getReweightedData(data,beta0,mu0,beta1,mu1,nWeighting):
    
    
    dataCopy=data.copy()
    
    weightArray=np.array([nWeighting[i] for i in data['N']])
    
    
        
    weights=-(beta1-beta0)*data['E']+(beta1*mu1-beta0*mu0)*data['N']+np.log(weightArray)
    
    dataCopy['w']=np.exp(weights)
    
    return dataCopy

#generate data following distribution but without carrying around weighting
def getUnweightedData(data):
    
    binWidth=5
    binEdges=np.arange(np.min(data['N'])-binWidth/2,np.max(data['N'])+binWidth/2,binWidth)
    hist = np.histogram(data['N'],weights=data['w'], density=True,bins=binEdges)    
    countVector=binEdges[1:]-binWidth/2
    
# #use date directly
#     unweightedPMF=stats.rv_discrete( values=(countVector, hist))
#     newData=unweightedPMF.rvs(size=30000)

#fit PDF
    
    hist_dist = scipy.stats.rv_histogram(hist)
    newData=hist_dist.rvs(size=200000)


    
    return newData


# In[3]:


# #given an empirical histogram, find the deltaMu which gives equal weights to peaks. Return the deltaMu and shifted histogram

#basic form of PDF
def gaussian(x,mu,sigma):
    prefactor=1/(sigma*np.sqrt(2*np.pi))
    exponential=np.exp(-(((x-mu)/sigma)**2)/2)
    return prefactor*exponential

#a and b are weights
def doubleGaussian(x,mu1,sigma1,mu2,sigma2,a):
    b=1-a
    return a*gaussian(x,mu1,sigma1)+b*gaussian(x,mu2,sigma2)

#given empirical pdf, fit to double gaussian. returns optimized parameters
def fitPDF(empiricalPDF):
    guess=[300,100,700,100,0.5]
    xData=list(empiricalPDF.keys())
    yData=list(empiricalPDF.values())
    params,cov=scipy.optimize.curve_fit(doubleGaussian,xData,yData,guess)
    return params

def getGMM(data):
    
    unweightedData=getUnweightedData(data)
    unweightedData=unweightedData.reshape(-1, 1)
    # fit a Gaussian Mixture Model with two components
    clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf.fit(unweightedData)
    
    return clf

def getDiluteFraction(data):
    
    #get pdf from data
#     binEdges=np.arange(np.min(data['N'])-0.5,np.max(data['N'])+0.5,1)
#     hist,_ = np.histogram(data['N'],weights=data['w'], density=True,bins=binEdges)    
    
#     countVector=binEdges[1:]-0.5
    
#     dataPDF=dict(zip(countVector, hist))
    
#     fitParams=fitPDF(dataPDF)
#     diluteFraction=fitParams[-1]


#Gaussian mixture model

    thisGMM=getGMM(data)
    means=thisGMM.means_
    means=[means[0,0],means[1,0]]
    
    weights=thisGMM.weights_
    diluteFraction=weights[np.argmin(means)]
    

    return diluteFraction


# In[4]:


#define a function mapping mu1 to dense fraction
def reweightedDiluteFraction(trialMu1,*paramTuple):

    #unpack the tuple containing parameters
    data,beta0,mu0,beta1=paramTuple
    
    reweightedData=getReweightedData(data,beta0,mu0,beta1,trialMu1,multicanonicalWeighting)
    diluteFraction=getDiluteFraction(reweightedData)
    
    
    return diluteFraction-0.5

def getEqualWeights(data,beta0,mu0,beta1):
    
#     solTuple=mu1=scipy.optimize.fsolve(reweightedDiluteFraction,mu0,args=(data,beta0,mu0,beta1),maxfev=1000,full_output=1)
#     mu1=solTuple[0]
#     mu1=mu1[0]
    
    
#     print(solTuple)
    leftEdge=mu0-0.01
    rightEdge=mu0+0.01
    brentCondition=False

    while(brentCondition==False):
        leftValue=reweightedDiluteFraction(leftEdge,data,beta0,mu0,beta1)
        rightValue=reweightedDiluteFraction(rightEdge,data,beta0,mu0,beta1)
        
        if np.sign(leftValue)==np.sign(rightValue):
            
            if leftValue<0:
                leftEdge=leftEdge-0.01     
                
            
            if rightValue>0:
                rightEdge=rightEdge+0.01
        
        else:
            brentCondition=True


    mu1=scipy.optimize.brentq(reweightedDiluteFraction,leftEdge,rightEdge,args=(data,beta0,mu0,beta1),xtol=0.0001,maxiter=5000)
    
    solutionData=getReweightedData(data,beta0,mu0,beta1,mu1,multicanonicalWeighting)
    return (mu1,solutionData)


def getMeans(data):
    
    #get pdf from data
#     binEdges=np.arange(np.min(data['N'])-0.5,np.max(data['N'])+0.5,1)
#     hist,_ = np.histogram(data['N'],weights=data['w'], density=True,bins=binEdges)    
    
#     countVector=binEdges[1:]-0.5
    
#     dataPDF=dict(zip(countVector, hist))
    
#     fitParams=fitPDF(dataPDF)
    
#     phi1,phi2=fitParams[1],fitParams[3]


#Gaussian mixture model

    thisGMM=getGMM(data)
    means=thisGMM.means_ 
    
    means=np.array([means[0,0],means[1,0]])
    means=np.sort(means)
        
    phi1=means[0]*length/volume
    phi2=means[1]*length/volume
    
    return phi1,phi2
    


# In[5]:


def getSurfaceTension(data,binSize):
    
    tensionHist,edges=np.histogram(data['N'],weights=data['w'],bins=np.arange(0,1000,binSize),density=True)
    
    prelimThreshold=500
    
    #find modes
    thresholdIndex=np.argwhere(edges>prelimThreshold)[0,0]
    
    
    mode1Index=np.argmax(tensionHist[0:thresholdIndex])
    mode2Index=np.argmax(tensionHist[thresholdIndex:])+thresholdIndex
    
    mode1=tensionHist[mode1Index]
    mode2=tensionHist[mode2Index]
    
    pMax=(mode1+mode2)/2
    
    #find min between modes
    pMin=np.min(tensionHist[mode1Index:mode2Index])
    
    prefactor=1/(2*beta1*(linearD**2))
    
    surfaceTension=prefactor*np.log(pMax/pMin)
    
    return surfaceTension


# In[6]:


def getThermalizedDataframe(data,finalBeta):
    #find step where we reach final temp
    cooledStep=0
    cooledTuple=np.nonzero((data.iloc[:,1]==finalBeta).values)
    if len(cooledTuple[0]>0):
        cooledStep=cooledTuple[0][0]

    #thermalize
    thermalizedData=data.iloc[cooledStep:]
    thermalizedData=thermalizedData.iloc[int(len(thermalizedData)/5):]
    
    return(thermalizedData)


# In[7]:


def getBetaPairs(seqArg):
    newSeqArg=seqArg+'_'
    betaPairs=[]
    allFiles=os.listdir(path)
    for thisFile in allFiles:
        if thisFile.find(newSeqArg)!=-1 and thisFile.find('preweighting')!=-1:
            loc1=thisFile.find('beta0')+6
            loc2=thisFile.find('beta1')-1
            loc3=thisFile.find('beta1')+6
            loc4=thisFile.find('.txt')
            
            thisBeta0=thisFile[loc1:loc2]
            thisBeta1=thisFile[loc3:loc4]
            
            betaPairs.append([float(thisBeta0),float(thisBeta1)])
            
    
    return betaPairs


# In[8]:


#specify sequence and temperatures
linearD=30
volume=linearD**3
j=0.05
# numReps=3
numReps=3

seq='L24_b2'
length=24


path='/rawData/'
betaPairList=getBetaPairs(seq)



#beta1, phi1 mean, phi1 SD, phi2 mean, phi2 SD
dataSummary=np.zeros([len(betaPairList),5])
plotting=False


# In[24]:


#load data and process

phaseLine=np.zeros([len(betaPairList),3]) #T mu list

for counter,betaPair in enumerate(betaPairList):
    
    beta0,beta1=betaPair

    phi1Reps=[]
    phi2Reps=[]
    muStarReps=[]
    
    for repNumber in range(numReps):
        filename='writeRecord_multicanonical_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'_j'+str(j)+'_rep'+str(repNumber)+'.dat'
        weightingFilename='preweighting_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'.txt'

        #load
        thisData=pd.read_csv(path+filename,sep=' ',header=None,names=['steps','beta','N','E'],index_col=False)

        #thermalize
        thermalizedData=getThermalizedDataframe(thisData,beta1)

        #load multicanonical weighting data. It will be just -1: -1 for uniform weighting

        multicanonicalWeightingArray=np.loadtxt(path+weightingFilename,skiprows=9,delimiter=' ')
        multicanonicalWeighting={}
        for i in range(multicanonicalWeightingArray.shape[0]):
            multicanonicalWeighting[int(multicanonicalWeightingArray[i,0])]=multicanonicalWeightingArray[i,1]

         #load mu1 from preweighting file
        weightingFileStream=open(path+weightingFilename,"r")
        weightingLines=weightingFileStream.readlines()
        mu1=float(weightingLines[7])

        #remove eta
        dataGCE=getReweightedData(thermalizedData,beta1,mu1,beta1,mu1,multicanonicalWeighting)
        

    #     #reweight to equal weights
        beta1MuCritical,beta1Sol=getEqualWeights(dataGCE,beta1,mu1,beta1)
        
        muStarReps.append(beta1MuCritical)

        #get phi overlap, surface tension
        phi1,phi2=getMeans(beta1Sol)
        diluteFraction=getDiluteFraction(beta1Sol)
        
        surfaceTension=getSurfaceTension(beta1Sol,10)
        
        phi1Reps.append(phi1)
        phi2Reps.append(phi2)
    
    phi1_mean=np.mean(phi1Reps)
    phi1_sd=np.std(phi1Reps)
    
    phi2_mean=np.mean(phi2Reps)
    phi2_sd=np.std(phi2Reps)
    
    dataSummary[counter]=np.array([beta1,phi1_mean,phi1_sd,phi2_mean,phi2_sd])
    
    
    muStar_mean=np.mean(muStarReps)
    muStar_std=np.std(muStarReps)
    phaseLine[counter]=np.array([muStar_mean,muStar_std,1/beta1])
    
#     print(surfaceTension)
    
    
    #export plots
    bins=np.arange(0,1000,5)
    if plotting:
        plt.figure()
        plt.hist(thermalizedData['N'],bins=bins,density=True)
        plt.xlabel('$N$')
        plt.ylabel('$P(N)$')
        plt.title(seq+', $\\beta=$'+str(beta1)+', $\\tilde{H}$')
#         plt.savefig(seq+'_beta1_'+str(beta1)+'_eta.png')

        plt.figure()
        plt.hist(dataGCE['N'],weights=dataGCE['w'],bins=bins,density=True)
        plt.xlabel('$N$')
        plt.ylabel('$P(N)$')
        plt.title(seq+', $\\beta=$'+str(beta1)+', $H$')
#         plt.savefig(seq+'_beta1_'+str(beta1)+'_noEta.png')

        plt.figure()
        plt.hist(beta1Sol['N'],weights=beta1Sol['w'],bins=bins,density=True)
        plt.xlabel('$N$')
        plt.ylabel('$P(N)$')
        plt.title(seq+', $\\beta=$'+str(beta1)+', equal weights'+'\n'+'$\\phi_1=$'+str(round(phi1,3))+', $\\phi_2=$'+str(round(phi2,3)))
#         plt.savefig(seq+'_beta1_'+str(beta1)+'_reweighted.png')


# In[44]:


#for ell=2:
criticalPoint=np.array([-10.488682202333335,1.0888366012595907])

plt.figure(figsize=(8,6))
plt.errorbar(phaseLine[:,0],phaseLine[:,2],yerr=None,xerr=phaseLine[:,1],markersize=12,marker='.',linestyle='-')
plt.plot(criticalPoint[0],criticalPoint[1],marker='*',markersize=15)
plt.xlim(-10.65,-10.43)
plt.ylim(1.061,1.095)
plt.xlabel('Chemical potential $\\mu$ $(\\epsilon)$',fontsize=22)
plt.ylabel('Temperature $T$ $(\\epsilon/k_\mathrm{B}$)',fontsize=22)
plt.tick_params(axis='both',labelsize=18)
# plt.savefig('phaseLine.svg',bbox_inches = "tight")


# In[ ]:


np.save('phaseData_'+seq+'_j'+str(j)+'.npy',dataSummary)


# In[ ]:


#example of multicanonical weighting for methods
if seq=='L24_b2':

    beta0,beta1=betaPairList[-1]

    phi1Reps=[]
    phi2Reps=[]

    for repNumber in range(numReps):
        filename='writeRecord_multicanonical_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'_j'+str(j)+'_rep'+str(repNumber)+'.dat'
        weightingFilename='preweighting_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'.txt'

        #load
        thisData=pd.read_csv(path+filename,sep=' ',header=None,names=['steps','beta','N','E'],index_col=False)

        #thermalize
        thermalizedData=getThermalizedDataframe(thisData,beta1)

        #load multicanonical weighting data. It will be just -1: -1 for uniform weighting

        multicanonicalWeightingArray=np.loadtxt(path+weightingFilename,skiprows=9,delimiter=' ')
        multicanonicalWeighting={}
        for i in range(multicanonicalWeightingArray.shape[0]):
            multicanonicalWeighting[int(multicanonicalWeightingArray[i,0])]=multicanonicalWeightingArray[i,1]

         #load mu1 from preweighting file
        weightingFileStream=open(path+weightingFilename,"r")
        weightingLines=weightingFileStream.readlines()
        mu1=float(weightingLines[7])

        #remove eta
        dataGCE=getReweightedData(thermalizedData,beta1,mu1,beta1,mu1,multicanonicalWeighting)


    #     #reweight to equal weights
        beta1MuCritical,beta1Sol=getEqualWeights(dataGCE,beta1,mu1,beta1)

        #get phi overlap, surface tension
        phi1,phi2=getMeans(beta1Sol)
        diluteFraction=getDiluteFraction(beta1Sol)

        surfaceTension=getSurfaceTension(beta1Sol,10)

        phi1Reps.append(phi1)
        phi2Reps.append(phi2)

    phi1_mean=np.mean(phi1Reps)
    phi1_sd=np.std(phi1Reps)

    phi2_mean=np.mean(phi2Reps)
    phi2_sd=np.std(phi2Reps)

    dataSummary[counter]=np.array([beta1,phi1_mean,phi1_sd,phi2_mean,phi2_sd])

#     print(surfaceTension)

    yMax=0.00925
    figSize=[8,6]
    labelSize=22
    tickSize=18
    
    
    #export plots
    bins=np.arange(0,1100,5)
    if plotting:
        plt.figure(figsize=(figSize[0],figSize[1]))
        plt.hist(thermalizedData['N'],bins=bins,density=True)
        plt.xlabel('Polymer number $N$',fontsize=labelSize)
        plt.ylabel('Probability $\\tilde{P}(N)$',fontsize=labelSize)
        plt.ylim(ymax=yMax)
        plt.tick_params(axis='both',labelsize=tickSize)
        plt.tight_layout()


        plt.savefig('methods_multicanonical_'+seq+'_beta1_'+str(beta1)+'_raw.svg')

        plt.figure(figsize=(figSize[0],figSize[1]))
        plt.hist(dataGCE['N'],weights=dataGCE['w'],bins=bins,density=True)
        plt.xlabel('Polymer number $N$',fontsize=labelSize)
        plt.ylabel('Probability $P(N)$',fontsize=labelSize)
        plt.ylim(ymax=yMax)
        plt.tick_params(axis='both',labelsize=tickSize)
        plt.tight_layout()

#         plt.savefig('methods_multicanonical_'+seq+'_beta1_'+str(beta1)+'_noEta.svg')

        plt.figure(figsize=(figSize[0],figSize[1]))
        plt.hist(beta1Sol['N'],weights=beta1Sol['w'],bins=bins,density=True)
        plt.xlabel('Polymer number $N$',fontsize=labelSize)
        plt.ylabel('Probability $P(N)$',fontsize=labelSize)
        plt.ylim(ymax=yMax)
        plt.tick_params(axis='both',labelsize=tickSize)
        plt.tight_layout()

#         plt.savefig('methods_multicanonical_'+seq+'_beta1_'+str(beta1)+'_reweighted.svg')


# In[ ]:




