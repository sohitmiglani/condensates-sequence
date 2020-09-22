#This script loads a sequence, calculates the condensation parameter chi, and use the linear chi-Tc fit to find the critical temperatur Tc.
#Tc is returned in units of epsilon/k_B, where epsilon is the specific binding energy.

#Fits were performed for Monte Carlo simulations on an FCC lattice with V=30^3 sites and nonspecific-interaction parameter J=0.05*epsilon
#The Tc values will be different for different systems, but the sequence-dependent hierarchy should be preserved.

#Sequences files have the following format:
#Line 1: header
#Line 2: length
#Line 3: the sequence, where A (B) motifs are represented as 1 (-1), separated by spaces.
#See the example file polySpecs_L24_l6.txt

#The only argument for the script is the name of the sequence file.


import numpy as np
import operator
import scipy.stats
import scipy.optimize
import sys


seqFile=sys.argv[1]


#calculate Chi

#params: g(1), correlation metric Pcorr, a/L
def getChi(g1,Pcorr,plusRatio):
    
    minusRatio=1-plusRatio
    normTerm=4*Pcorr
    
    renormalization=normTerm**(1/2) 
    entropy=g1/renormalization
    
    nPlus=length*plusRatio
    nMinus=length*minusRatio
    
    overlapTerm=1/(((minusRatio)**nPlus)*((plusRatio)**nMinus))
    
    return -np.log(entropy*overlapTerm)


#use chi and fit to calculate Tc
#argument: numpy array specifying sequence
#returns: chi, Tc
def getCriticalTemperature(seqArray):
    
    thisG1=getG1(seqArray,loopWeight)
    thisPcorr=getPcorr(seqArray)
    
    pluses=seqArray[seqArray==1]
    thisPlusRatio=len(pluses)/len(seqArray)
    
    thisChi=getChi(thisG1,thisPcorr,thisPlusRatio)
    
    thisTc=slope*thisChi+intercept
    
    
    return thisChi,thisTc


#functions called by getChi
#*************************************************************************************************

#calculate g(1) using SAW and loop entropy


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


#args: sequence to fit to using g(1)
#returns: a1
def fitLoopWeight(fitSeqBlockSize):
    
    filename='dos_L24_block'+str(fitSeqBlockSize)+'.out' 
    fittingData=np.loadtxt(dosPath_fitting+filename)
    
    
    datag0=fittingData[0]
    datag1=fittingData[1]/datag0
    datag2=fittingData[2]/datag0
    

    fitSeqArray=np.loadtxt(path_seqs+'polySpecs_L24_b'+str(fitSeqBlockSize)+'.txt',dtype=int,skiprows=2)

    
    pairs=getPairs(fitSeqArray)
    
    sol = scipy.optimize.root_scalar(lambda x:g1Norm(pairs,x)-datag1,x0=0.1,x1=100)

        
    return sol.root




#calculate correlation metric Pcorr

#assume a bond has formed
#choose 1 neighbor from each protein
#what is probability that these are compatible?
#argument: seq as array

def getPcorr(seq):
    
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


#parameters for scaling relations. Mu is lattice dependent.
#******************************************************************************************

#first, fit to known g(s) to determine loop weight

#exponents and scaling parameters
gamma=1.157
nu=0.588
muDict={'fcc':10.037}

lattice='fcc'
mu=muDict[lattice]



#loopweight fit to g(1) from domain size=12 polymer
loopWeight=0.7298937942762107



#chi-Tc fit using block sequences
slope=0.10894316705597923 
intercept=2.976691130509688


###################################################################################################

#load sequence and perform calculation
testSeq=np.loadtxt(seqFile,dtype=int,skiprows=2)
length=int(len(testSeq))
chi,tc=getCriticalTemperature(testSeq)
print(chi,tc)






