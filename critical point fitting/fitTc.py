
# coding: utf-8

# given a histogram of p(E,N) at mu0,beta0, reweight it and fit to Ising universality distribution to find Tc

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.optimize
import scipy.integrate
import scipy.stats
import time
import sys


# In[2]:


#m is the order parameter, m0 is a scale dependent factor. we normalize to unit variance. see https://arxiv.org/pdf/cond-mat/9909343.pdf
a=0.1582
c=0.7762
m0=1.1342

def isingDist(m):
	normalization=1/(m0*2.0546)
	argument=-(((m/m0)**2-1)**2)*(a*(m/m0)**2+c)
	return normalization*np.exp(argument)

def getIsingHistogram(edges):
	histogram=np.zeros([edges.shape[0]-1])
	
	for i in range(len(edges)-1):
		histogram[i]=scipy.integrate.quad(isingDist, edges[i], edges[i+1])[0]
		
	#renormalize to mimic the pdf
	histogram=histogram/(edges[1]-edges[0])
		
	return histogram


# In[3]:


#simulate sampling from ising order parameter distribution
# class your_distribution(scipy.stats.rv_continuous):
#     def _pdf(self, x):
#         return isingDist(x)

# distribution = your_distribution()

# isingSim=np.zeros([1000])
# for i in range(1000):
#     isingSim[i]=distribution.rvs()


# In[4]:


#load data
seq=sys.argv[1]
linearD=30
beta1=float(sys.argv[2])
rep=sys.argv[3]
j=0.05

path='/tigress/bweiner/idp/grandCanonical/multicanonical_8.13.20/run0_multiL_reps/'
filename='writeRecord_multicanonical_etaUniform_'+seq+'_e'+str(beta1)+'_j'+str(j)+'_rep'+str(rep)+'.dat'
weightingFilename='preweighting_'+seq+'_uniform_beta1_'+str(beta1)+'.txt'


thisData=pd.read_csv(path+filename,sep=' ',header=None,names=['steps','beta','N','E'])

#find step where we reach final temp
cooledStep=0
cooledTuple=np.nonzero(thisData.iloc[:,1]==beta1)
if len(cooledTuple[0]>0):
	cooledStep=cooledTuple[0][0]
	
 #load mu1 from preweighting file
weightingFileStream=open(path+weightingFilename,"r")
weightingLines=weightingFileStream.readlines()
mu1=float(weightingLines[7])

#thermalize
thermalizedData=thisData.iloc[cooledStep:]
thermalizedData=thermalizedData.iloc[int(len(thermalizedData)/10):]
thermalizedData['w']=np.full((thermalizedData.shape[0],1),1)




# In[5]:


#using full data set is prohibitively slow, so we subsample
sparseFilter=[]
sparseness=1
for i in range(thermalizedData.shape[0]):
	if i%sparseness==0:
		sparseFilter.append(True)
	else:
		sparseFilter.append(False)
sparseData=thermalizedData[sparseFilter]


# In[6]:


#s is mxing parameter. Return two numpy arrays
def getFieldMixing(data,beta0,mu0,beta1,mu1,s):
	
	
	weights=-(beta1-beta0)*data['E']+(beta1*mu1-beta0*mu0)*data['N']
	weights=np.exp(weights)
	
	mixedField=data['N']-s*data['E']
	
	
	return mixedField,weights
	
#now obtain a histogram of mixed field with unit variance. Returns single numpy array
def getNormalizedPDF(mixedField,weights):
	

	try:
		#subtract (weighted) mean of mixed field
		meanMixedField=np.average(mixedField,weights=weights)
		shiftedField=mixedField-meanMixedField
		
	 
		#rescale to have unit variance
		varianceScale=scipy.optimize.fsolve(varianceFunction,100,args=(shiftedField,weights)) 
		unitVarianceField=shiftedField/abs(varianceScale)
		
		

	except ZeroDivisionError:
		unitVarianceField=np.zeros([mixedField.shape[0]])

	

	return unitVarianceField

	
def varianceFunction(varianceScale,samples,weights):
	
	scaledData=samples/varianceScale
	
	weightedMean=np.average(scaledData,weights=weights)
	variance=np.average((scaledData-weightedMean)**2,weights=weights)
	
	return variance-1
	
	


# In[7]:


#compare histogram to Ising distribution
isingEdges=np.arange(-2.5,2.505,.05)
isingBinMidpoints=isingEdges[1:]-(isingEdges[1]-isingEdges[0])/2
isingHistogram=getIsingHistogram(isingEdges)

def getSquareResidual(fieldData,fieldWeights,edges,isingHist):
	
	try:
		dataHist,_=np.histogram(fieldData,bins=edges,weights=fieldWeights,density=True)
		
		residual=dataHist-isingHist
		squareResidual=np.power(residual,2)
		residualSum=np.sum(squareResidual)
		
		#chi squared:
	#     squareResidual=np.divide(squareResidual,isingHist)
		
	
	except ValueError:
		residualSum=10

	return residualSum


# In[8]:


#define transformations for s, so we can efficiently search over the desired parameter ranges
def sTransform(x,a,b,sign):
	return sign*a*x**b

#x1 and x2 are desired solve parameter bounds, generally 0.045-0.055. s1 and s2 are the desired real parameter bounds
def sTransformParams(x1,x2,s1,s2):
	b=np.log(s2/s1)/np.log(x2/x1)
	
	logA=np.log(s1)-b*np.log(x1)
	
	a=np.exp(logA)
	
	return a,b


# In[9]:


#solve for mu1,beta1,s which gives best match to Ising data
def functionToMinimize(testParams,a,b,sign):
	
	testMu=testParams[0]
	testBeta=testParams[1]
	
	testS=sTransform(testParams[2],a,b,sign)


#     mixedField_unormalized,testWeights=getFieldMixing(thermalizedData,beta1,mu1,testBeta,testMu,testS)
	mixedField_unormalized,testWeights=getFieldMixing(sparseData,beta1,mu1,testBeta,testMu,testS)


	testField=getNormalizedPDF(mixedField_unormalized,testWeights)
	
	thisSquareResidual=getSquareResidual(testField,testWeights,isingEdges,isingHistogram)
	
	return thisSquareResidual

#deltaMu and deltaBeta: deviation of initial condition from simulation value
#a,b: parameters for power function scaling of mixing parameter s
#sign: +1/-1 sign of s
def solveForCriticalPoint(deltaMu,deltaBeta,a,b,sign):
	
	x0=[mu1+deltaMu,beta1+deltaBeta,0.05]

	sol=scipy.optimize.minimize(functionToMinimize,x0,args=(a,b,sign),method='Nelder-Mead')
	
	#get results and process them
	solMu,solBeta,solSParam=sol.x
	solS=sTransform(solSParam,a,b,sign)
	
	residual=sol.fun
	
	return solMu,solBeta,solS,residual


# In[10]:


sRanges=[[0.001,0.01],[0.01,0.1],[0.1,1]]

initialConditions=[]
for deltaMu in np.arange(-0.05,0.06,0.01):
	for deltaBeta in np.arange(-0.002,0.003,0.001):
		for sPair in sRanges:
			for signs in [-1,1]:
				
				theseParams=[round(deltaMu,3),round(deltaBeta,3),sPair[0],sPair[1],signs]
				initialConditions.append(theseParams)



outFile='tcFit_'+seq+'_j'+str(j)+'_rep'+str(rep)+'.dat'


# In[13]:

for theseConditions in initialConditions:
	a,b=sTransformParams(0.045,0.055,theseConditions[2],theseConditions[3])
	solMu,solBeta,solS,residual=solveForCriticalPoint(theseConditions[0],theseConditions[1],a,b,theseConditions[4])
	
	f = open(outFile, "a")
	f.write(str(solMu)+' '+str(solBeta)+' '+str(solS)+' '+str(residual)+'\n')
	f.close()




