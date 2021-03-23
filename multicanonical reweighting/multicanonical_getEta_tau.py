

# given a histogram of p(E,N) at mu0,beta0, reweight it and return an estimate of p(N) at mu1,beta1
#however, we want (T-Tc)/Tc to be the same for all sequences, so we provide betaC=1/Tc

import sys
import numpy as np
import pandas as pd
import scipy.optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# convention:
#     beta0=temp of the simulation from which the current weighting was extracted (not meaningful for uniform)
#     beta1=temp of the actual simulation being analyzed
#     betaNew=temp of the simulation we want to do next, after reweighting
#     betaNew=beta1+deltaBeta



#begin loading data
seq=sys.argv[1]
betaC=float(sys.argv[2])

unweighted=int(sys.argv[3])
beta0=float(sys.argv[4])
beta1=float(sys.argv[5])

#tau=(T-Tc)/Tc
tau=float(sys.argv[6])

tc=1/betaC
tNew=tc*(tau+1)
betaNew=1/tNew


linearD=30 #lattice size
j=0.05
length=24


if unweighted==1:    #we're starting from an unweighted simulation
	
	path='/tigress/bweiner/idp/grandCanonical/multicanonical_6.29.20/block_uniform/'
	filename='writeRecord_multicanonical_etaUniform_'+seq+'_e'+str(beta1)+'_j'+str(j)+'_rep0.dat'
	weightingFilename='preweighting_'+seq+'_uniform_beta1_'+str(beta1)+'.txt'


	
else: #we're loading a simulation where we've weighted H according to higher T histogram

	path='/tigress/bweiner/idp/grandCanonical/multicanonical_6.29.20/block_tauStep1/'
	filename='writeRecord_multicanonical_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'_j'+str(j)+'.dat'
	weightingFilename='preweighting_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'.txt'
	



#load data


thisData=pd.read_csv(path+filename,sep=' ',header=None,names=['steps','beta','N','E'])

#find step where we reach final temp
cooledStep=0
cooledTuple=np.nonzero(thisData.iloc[:,1]==beta1)
if len(cooledTuple[0]>0):
	cooledStep=cooledTuple[0][0]

#thermalize
thermalizedData=thisData.iloc[cooledStep:]
thermalizedData=thermalizedData.iloc[int(len(thermalizedData)/5):]
thermalizedData['w']=np.full((thermalizedData.shape[0],1),1)

#load multicanonical weighting data. It will be just -1: -1 for uniform weighting

multicanonicalWeightingArray=np.loadtxt(path+weightingFilename,skiprows=9,delimiter=' ')
multicanonicalWeighting={}
for i in range(multicanonicalWeightingArray.shape[0]):
	multicanonicalWeighting[int(multicanonicalWeightingArray[i,0])]=multicanonicalWeightingArray[i,1]

 #load mu1 from preweighting file
weightingFileStream=open(path+weightingFilename,"r")
weightingLines=weightingFileStream.readlines()
mu1=float(weightingLines[7])
	





def getReweightedData(data,beta0,mu0,beta1,mu1,nWeighting):
    
    
    dataCopy=data.copy()
    
    weightArray=np.array([nWeighting[i] for i in data['N']])
    
        
    weights=-(beta1-beta0)*data['E']+(beta1*mu1-beta0*mu0)*data['N']+np.log(weightArray)
    
    dataCopy['w']=np.exp(weights)
    
    return dataCopy
		



#get data with multicanonical weighting removed
dataGCE=getReweightedData(thermalizedData,beta1,mu1,beta1,mu1,multicanonicalWeighting)



#given a pdf which isn't clearly separated, return the (pdf,mu1) at a new temp which has equal weights
#define equal weights using halfway between modes as threshold


# get dilute fraction by setting threshold at halfway between 2 peaks
def overlappingDiluteFraction(data):    
	
	prelimThreshold=500
	
	binEdges=np.arange(np.min(data['N'])-0.5,np.max(data['N'])+0.5,1)
	hist,_ = np.histogram(data['N'],weights=data['w'], density=True,bins=binEdges)

	
	#find modes
	countVector=binEdges[1:]-0.5
	thresholdIndex=np.argwhere(abs(countVector-prelimThreshold)<10**-5)[0,0]

	
	#in histogram bin index
	mode1Index=np.argmax(hist[0:thresholdIndex])
	mode2Index=np.argmax(hist[thresholdIndex:])+thresholdIndex
	

	realThresholdIndex=int(round((mode1Index+mode2Index)/2))
	
	
	dilute=hist[0:realThresholdIndex]    
	diluteFraction=np.sum(dilute)
	
	return diluteFraction

	
#define a function mapping mu1 to dense fraction
def reweightedDiluteFraction(trialMu1,*paramTuple):

	#unpack the tuple containing parameters
	data,beta0,mu0,beta1=paramTuple
	
	reweightedData=getReweightedData(data,beta0,mu0,beta1,trialMu1,multicanonicalWeighting)
	reweightedFraction=overlappingDiluteFraction(reweightedData)
	
	return reweightedFraction-0.5

def getEqualWeights_overlapping(data,beta0,mu0,beta1):
	
	mu1=scipy.optimize.fsolve(reweightedDiluteFraction,mu0-0.005,args=(data,beta0,mu0,beta1))
	
	solutionData=getReweightedData(data,beta0,mu0,beta1,mu1,multicanonicalWeighting)
	return (mu1[0],solutionData)




#if you have data where the peaks overlap, just estimate means using the midway point between modes
def getMeans_overlapping(data):
	
	binEdges=np.arange(np.min(data['N'])-0.5,np.max(data['N'])+0.5,1)
	hist,_ = np.histogram(data['N'],weights=data['w'], density=True,bins=binEdges)
	
	prelimThreshold=400
	
	#find modes
	countVector=binEdges[1:]-0.5
	thresholdIndex=np.argwhere(abs(countVector-prelimThreshold)<10**-5)[0,0]

	
	#in histogram bin index
	mode1Index=np.argmax(hist[0:thresholdIndex])
	mode2Index=np.argmax(hist[thresholdIndex:])+thresholdIndex
	

	realThresholdIndex=int(round((mode1Index+mode2Index)/2))
	
	diluteData=data[data['N']<=countVector[realThresholdIndex]]
	denseData=data[data['N']>countVector[realThresholdIndex]]
	
	diluteMean=np.average(diluteData['N'],weights=diluteData['w'])
	denseMean=np.average(denseData['N'],weights=denseData['w'])
	
	phi1=diluteMean*length/(linearD**3)
	phi2=denseMean*length/(linearD**3)

	
	return phi1,phi2
	




#solve for equal weights at a new temperature, given overlapping distribution

solMu,solData=getEqualWeights_overlapping(dataGCE,beta1,mu1,betaNew)
solDiluteFraction=overlappingDiluteFraction(solData)




plt.hist(dataGCE['N'],weights=dataGCE['w'],bins=np.arange(10,1200,1),density=True,label='old')
plt.hist(solData['N'],weights=solData['w'],bins=np.arange(10,1200,1),density=True,label='new')
plt.legend()
plt.savefig('reweightingSolution_'+seq+'_beta1_'+str(beta1)+'_betaNew_'+str(betaNew)+'.png')



#convert dataframe to pdf
binEdges=np.arange(np.min(solData['N'])-0.5,np.max(solData['N'])+0.5,1)
H,_=np.histogram(solData['N'],weights=solData['w'],bins=binEdges,density=True)

countVector=binEdges[1:]-0.5
solPDF=np.array([countVector,H]).T




#smooth reweighting, truncate noisy tails
from scipy.signal import savgol_filter
smoothedData = savgol_filter(solPDF[:,1], 99, 3) # window size 99, polynomial order 3

pdfBoundaries=np.argwhere(smoothedData>0.001).flatten()

# plt.plot(solPDF[:,0],solPDF[:,1])
# plt.plot(solPDF[pdfBoundaries[0]:pdfBoundaries[-1],0],smoothedData[pdfBoundaries[0]:pdfBoundaries[-1]])

smoothPDF=np.array([solPDF[pdfBoundaries[0]:pdfBoundaries[-1],0],smoothedData[pdfBoundaries[0]:pdfBoundaries[-1]]]).T

# plt.figure()
# plt.plot(smoothPDF[:,0],smoothPDF[:,1])




#add uniform tails
minWeight=smoothPDF[0]
maxWeight=smoothPDF[-1]
maxNp=2*(linearD**3)/length

leftTail=np.array([np.arange(0,minWeight[0]),np.full((int(minWeight[0]),),minWeight[1])]).T
rightTail=np.array([np.arange(maxWeight[0]+1,maxNp+1),np.full((int(maxNp-maxWeight[0]),),maxWeight[1])]).T

finalWeightingPDF=np.concatenate((leftTail,smoothPDF,rightTail),axis=0)






# #write output file
distType='beta0_'+str(beta1)+'_beta1_'+str(betaNew)
filename='preweighting_'+seq+'_'+distType+'.txt'
with open(filename, 'w') as f:
	f.write('beta0:\n')
	f.write(str(beta1)+'\n')
	f.write('mu0:\n')
	f.write(str(mu1)+'\n')
	f.write('beta1: \n')
	f.write(str(betaNew)+'\n')
	f.write('mu1:\n')
	f.write(str(solMu)+'\n')
	f.write('N, p(N):\n')
	for i in range(len(finalWeightingPDF)):
		f.write(str(int(finalWeightingPDF[i,0]))+" "+str(finalWeightingPDF[i,1])+'\n')




#uniform weighting
# maxNp=2*(linearD**3)/length
# uniformWeighting=np.array([np.arange(0,maxNp+1),np.full((int(maxNp)+1,),1)]).T

# distType='uniform'
# filename='preweighting_'+seq+'_'+distType+'.txt'
# with open(filename, 'w') as f:
#     f.write('beta0:\n')
#     f.write(str(0)+'\n')
#     f.write('mu0:\n')
#     f.write(str(0)+'\n')
#     f.write('beta1: \n')
#     f.write(str(0)+'\n')
#     f.write('mu1:\n')
#     f.write(str(solMu)+'\n')
#     f.write('N, p(N):\n')
#     for i in range(len(uniformWeighting)):
#         f.write(str(int(uniformWeighting[i,0]))+" "+str(uniformWeighting[i,1])+'\n')

