#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import operator
import scipy.stats
import scipy.optimize


# In[2]:


#calculate heuristic

def getOrderParameter(g1,Pbond,plusRatio):
    
    Nd=4*Pbond
    minusRatio=1-plusRatio
    
    
    normTerm=4*Pbond
    renormalization=normTerm**(1/2) 


    entropy=g1/renormalization
    

    
    nPlus=length*plusRatio
    nMinus=length*minusRatio
    
    overlapTerm=1/(((minusRatio)**nPlus)*((plusRatio)**nMinus))
    
    return np.log(entropy*overlapTerm)
    


# In[3]:


#calculate g(1)


#auxiliary functions called by g(1), g(2)
#argument: sequence as array
#returns: two lists: contiguous bonds and non-contiguous bonds
def getPairs(seqArray):
   
    seqList=seqArray.tolist()
    
    allBonds=[]
    
    for motif1Position,motif1 in enumerate(seqList):
        if motif1==1:
            for motif2Position,motif2 in enumerate(seqList):
                if motif2==-1:
                    allBonds.append([motif1Position,motif2Position])
                    
    contiguousBonds=[]
    loopBonds=[]
    
    for thisBond in allBonds:
        if abs(thisBond[0]-thisBond[1])==1:
            contiguousBonds.append(thisBond)
        else:
            loopBonds.append(thisBond)
                    
    return [contiguousBonds,loopBonds]


#argument: 2 motifs in loop, the contiguous bond
#returns: number of configs corresponding to that loop, apart from scaling amplitude
# w_walks*w_loop/A1
def getLoopConfigs(motifPair,contiguousBond=[]):
    
    #n=|i-j|
    n=abs(motifPair[0]-motifPair[1])
    
    #no contiguous bonds: the loop is n=|i-j| steps and the walk is L-1-n steps
    if len(contiguousBond)==0:
        
        walkContribution=(mu**(-1))*((length-1-n)/(length-1))**(gamma-1)
        loopContribution=(n-1)**(-3*nu)
        
        totalContribution=walkContribution*loopContribution

        
    #figure out if contigous bond shortens the loop or the walk
    else:  
        motifPair.sort()
        contiguousBond.sort()
        if contiguousBond[0]>motifPair[0] and contiguousBond[1]<motifPair[1]: #bond is inside loop, so loop is n-1 steps
            
            walkContribution=(mu**(-2))*((length-1-n)/(length-1))**(gamma-1)
            loopContribution=(n-2)**(-3*nu)
            
            totalContribution=walkContribution*loopContribution

            
        else:                                                         #bond is outside loop,so walk is shorter
            
            walkContribution=(mu**(-2))*((length-n-2)/(length-1))**(gamma-1)
            loopContribution=(n-1)**(-3*nu)
            
            totalContribution=walkContribution*loopContribution

    
    
        
    return totalContribution





#arg: set of contiguous bonds
#returns: number of ways to choose 2 bonds which don't conflict
def getNumberOfContigBondPairs(contigBonds):
    
    numPairs=0
    
    for bond1 in contigBonds:
        for bond2 in contigBonds:
            
            intersection=list(set(bond1) & set(bond2))
            if len(intersection)==0:
                numPairs+=1
    
    
    #divide by 2 to avoid double-counting
    return numPairs/2

#arg: a particular contiguous bond, the set of loop bonds
#returns: all loop bonds that don't involve contig motifs
def getAllowedLoops(contigBond,loopBonds):
    
    allowedLoops=[]
    
    for thisLoop in loopBonds:
        intersection=list(set(contigBond)&set(thisLoop))
        if len(intersection)==0:
            allowedLoops.append(thisLoop)
    
    return allowedLoops


#argument: pairs, loop parameter
#returns: g(1)/g(0)=contigbonds+a1*loops
def g1Norm(pairs,a1):
    
    contigBonds=pairs[0]
    loopBonds=pairs[1]
    
    contigBondConfigs=len(contigBonds)*(mu**(-1))*((length-2)/(length-1))**(gamma-1)
    
    loopConfigs=0
    
    for thisBond in loopBonds:
        thisLoopConfig=getLoopConfigs(thisBond)
        loopConfigs+=thisLoopConfig

    totalConfigs=contigBondConfigs+a1*loopConfigs
    
    return totalConfigs

def getG1(seq,a1):
    
    thesePairs=getPairs(seq)
    
    g1=g1Norm(thesePairs,a1)
    
#     return g1/g0
    return g1



#fit to data to determine loopWeight


#args: sequence to fit to, whether you want g(1) or g(2)
#returns: a1
def fitLoopWeight(fitSeqBlockSize,fitS):
    
    filename='dos_L24_block'+str(fitSeqBlockSize)+'.out' 
    fittingData=np.loadtxt(dosPath_fitting+filename)
    
    
    datag0=fittingData[0]
    datag1=fittingData[1]/datag0
    datag2=fittingData[2]/datag0
    

    fitSeqArray=np.loadtxt(path_seqs+'polySpecs_L24_b'+str(fitSeqBlockSize)+'.txt',dtype=int,skiprows=2)

    
    pairs=getPairs(fitSeqArray)
    
    if fitS==1:
        sol = scipy.optimize.root_scalar(lambda x:g1Norm(pairs,x)-datag1,x0=0.1,x1=100)
    else:
        sol = scipy.optimize.root_scalar(lambda x:g2Norm(pairs,x)-datag2,x0=0.1,x1=100)
        
    return sol.root


# In[ ]:





# In[4]:


#first, fit to known g(s) to determine loop weight

#exponents and scaling parameters
gamma=1.157
nu=0.588
muDict={'fcc':10.037}
length=24

lattice='fcc'
mu=muDict[lattice]

dosPath_fitting='densityOfStates/'
path_seqs='sequences/'

fitSeqBlockSize=12
fitS=1
loopWeight=fitLoopWeight(fitSeqBlockSize,fitS)


# In[5]:


#calculate correlation metric Pbond

#now, assume a bond has formed
#choose 1 neighbor from each protein
#what is probability that these are compatible?
#argument: seq as array

def getCompatibleNeighborProb(seq):
    
    length=len(seq)

    numMatches=0
    numStates=0

    for i in range(len(seq)):
        for j in range(len(seq)):

            protein1Monomer=seq[i]
            protein2Monomer=seq[j]
            
            if protein1Monomer==-1*protein2Monomer: #proceed if these two monomers lead to bond
                
                #get neighbors who might participate in bonds
                protein1NeighborList=[]
                protein2NeighborList=[]
                
                for protein1NeighborDirection in [-1,1]:
                    if i+protein1NeighborDirection>-1 and i+protein1NeighborDirection<length: #check for end monomers
                        protein1NeighborList.append(seq[i+protein1NeighborDirection])
                        
                for protein2NeighborDirection in [-1,1]:
                    if j+protein2NeighborDirection>-1 and j+protein2NeighborDirection<length: #check for end monomers
                        protein2NeighborList.append(seq[j+protein2NeighborDirection])
            
                for neighbor1 in protein1NeighborList:
                    for neighbor2 in protein2NeighborList:
                        numStates+=1
                        
                        if neighbor1==-1*neighbor2:
                            numMatches+=1
                

    return numMatches/numStates


# In[6]:


#block sequences


blockSizeList=[1,2,3,4,6,12]
j='0.05'

numReps=3
dataPairs_block=np.zeros([len(blockSizeList),4])

path_critical='tcData/tcData_L24_block/'


for counter,blockSize in enumerate(blockSizeList):
    
    
    seqArray=np.loadtxt(path_seqs+'polySpecs_L24_b'+str(blockSize)+'.txt',dtype=int,skiprows=2)
    
    g1=getG1(seqArray,loopWeight)
    thisPbond=getCompatibleNeighborProb(seqArray)
    
    
    repList=[]

    #get mean Tc data
    for rep in range(numReps):

        criticalFilename='tcFit_L24_b'+str(blockSize)+'_j'+str(j)+'_rep'+str(rep)+'.npy'
        thisCriticalData=np.load(path_critical+criticalFilename)
        

        betac=thisCriticalData[0]
        repList.append(betac)
    
    meanTc=np.mean(repList)
    sdTc=np.std(repList)

    thisOrderParameter=getOrderParameter(g1,thisPbond,0.5)

    dataPairs_block[counter]=np.array([-thisOrderParameter,1/meanTc,sdTc,thisPbond])


plt.figure(figsize=(8,6))

for i in range(len(blockSizeList)):
    plt.plot(dataPairs_block[i,0],dataPairs_block[i,1],'.',markersize=10,label=str(blockSizeList[i]))

plt.legend(title='Block size',title_fontsize=15,fontsize=15)
plt.ylabel('$T_\\mathrm{c}$',fontsize=18)

plt.xlabel('Condensation parameter $\Psi$',fontsize=18)
plt.tick_params(axis='both',labelsize=15)

#fit
slope,intercept,rValue,_,_=scipy.stats.linregress(dataPairs_block[:,0],dataPairs_block[:,1])

fitY=intercept+slope*dataPairs_block[:,0]
plt.plot(dataPairs_block[:,0],fitY,label='$R^2=$'+str(round(rValue,8)))

plt.savefig('orderParam_block_g1Calc_fit.svg')
print(slope,intercept)


# In[7]:


#generate random sequences and calculate order parameter

numSeqs=20000
orderParams=np.zeros(numSeqs)
pBonds=np.zeros(numSeqs)
baseSeq=np.concatenate((np.full(int(length/2),1,dtype=int),np.full(int(length/2),-1,dtype=int)))

for i in range(numSeqs):
    np.random.shuffle(baseSeq)
    
    thisG1=getG1(baseSeq,loopWeight)
    thisPbond=getCompatibleNeighborProb(baseSeq)
    
    thisOrderParam=getOrderParameter(thisG1,thisPbond,0.5)
    
    orderParams[i]=-thisOrderParam
    pBonds[i]=thisPbond


# In[8]:


#calculate Tc according to fit
tcSamples=np.zeros(numSeqs)
for i in range(len(tcSamples)):
    
    tcSamples[i]=slope*orderParams[i]+intercept


# In[9]:


plt.figure(figsize=(8,6))

plt.hist(tcSamples,bins=np.arange(1.05,1.32,0.01),density=True)
plt.xlabel('Critical temperature $T_\mathrm{c}$',fontsize=24)
plt.ylabel('$P(T_\mathrm{c})$',fontsize=24)
plt.tick_params(axis='both',labelsize=18)


#add points for block sequences
markerset=['.','_','^','s','*','o']
markersizes=[16,12,10,8,14,10]
blockY=45
blockPoints=np.array([dataPairs_block[:,1],np.full(6,blockY)]).T
for i in range(len(blockSizeList)):
    if blockSizeList[i]==12:
        plt.plot(blockPoints[i,0],blockPoints[i,1],marker=markerset[i],markersize=markersizes[i],markerfacecolor='none',color='k')
    
    elif blockSizeList[i]==2:
        plt.plot(blockPoints[i,0],blockPoints[i,1],marker=markerset[i],mew=3,markersize=markersizes[i],markerfacecolor='none',color='k')

    else:
        plt.plot(blockPoints[i,0],blockPoints[i,1],marker=markerset[i],markersize=markersizes[i],color='k')
        
    plt.axvline(blockPoints[i,0],ymax=0.9,color='k',linestyle='--')

# plt.savefig('fig3_orderParamDist.svg')


# In[10]:


plt.figure(figsize=(8,6))

plt.hist(pBonds,bins=np.arange(0.49,1,0.025),density=True)
plt.xlabel('Correlation metric $\\langle P_\mathrm{corr} \\rangle$',fontsize=22)
plt.ylabel('Probability of $\\langle P_\mathrm{corr} \\rangle$',fontsize=22)
plt.tick_params(axis='both',labelsize=18)


# #add points for block sequences
markerset=['.','_','^','s','*','o']
markersizes=[16,12,10,8,14,10]
blockY=26
blockPoints=np.array([dataPairs_block[:,-1],np.full(6,blockY)]).T
for i in range(len(blockSizeList)):
    if blockSizeList[i]==12:
        plt.plot(blockPoints[i,0],blockPoints[i,1],marker=markerset[i],markersize=markersizes[i],markerfacecolor='none',color='k')
    
    elif blockSizeList[i]==2:
        plt.plot(blockPoints[i,0],blockPoints[i,1],marker=markerset[i],mew=3,markersize=markersizes[i],markerfacecolor='none',color='k')

    else:
        plt.plot(blockPoints[i,0],blockPoints[i,1],marker=markerset[i],markersize=markersizes[i],color='k')
        
    plt.axvline(blockPoints[i,0],ymax=0.9,color='k',linestyle='--')

# plt.savefig('s_corr_pBondDist.svg')


# In[ ]:




