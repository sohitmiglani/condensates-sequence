#open a series of snapshot files, where the first line specifies number of each protein species, then lists coords in order
import numpy as np
import pandas as pd
from mayavi.mlab import *
from mayavi import mlab
import os
import sys

seq=sys.argv[1]
phi=sys.argv[2]
record=sys.argv[3]

#read data
##############################################################################

path='/snapshots/'


#multicanonical version
stepFile='polyOutput_multicanonical_'+seq+'_beta0_0.9267_beta1_0.9287_j0.05_record_'+record+'_v2.dat'



highRes=True
first=0


polyFilename='polySpecs_'+seq+'.txt'

#load polymer specification file and get sequence
sequence=np.loadtxt(path+polyFilename,dtype=int,skiprows=2)

#now read positions


#a list of monomer coordinates
timeSteps=[]


polyOutputFileName=path+stepFile

polyOutputFile=open(polyOutputFileName)
polyOutputLines=polyOutputFile.readlines()
proteinNumber=np.fromstring(polyOutputLines[1],dtype=int,sep=' ')


fullArray=pd.read_csv(polyOutputFileName,sep='\s+',header=None,skiprows=2,dtype=float)
fullArray=fullArray.values.astype(float)


#coordinate transformation to fcc lattice:
for i in range(fullArray.shape[0]):
	fullArray[i][0]=fullArray[i][0]-0.5*fullArray[i][1]+0.5*fullArray[i][2]
	fullArray[i][1]=(np.sqrt(3)/2)*fullArray[i][1]-(np.sqrt(3)/6)*fullArray[i][2]
	fullArray[i][2]=np.sqrt(2.0/3.0)*fullArray[i][2]

timeSteps.append(fullArray)


if highRes==True:
	drawLines=True
	iList=[first]



for i in iList:
	#change the color for overlapping points
	_, idx,inverse,counts = np.unique(timeSteps[i], return_index=True,return_inverse=True,return_counts=True,axis=0)

	#create monomer list for arbitrary protein numbers
	thisSequence=np.tile(sequence,proteinNumber)



	for k in range(counts.size):
		if counts[k]!=1:
			overlapping=np.where(inverse==k)
			thisSequence[overlapping]=0

	pointCoords=np.transpose(timeSteps[i])
	pointsObject=mlab.points3d(pointCoords[0],pointCoords[1],pointCoords[2],thisSequence,scale_mode='none',vmax=1,vmin=-1,scale_factor=0.5)

	#modify colors
	# lut = pointsObject.module_manager.scalar_lut_manager.lut.table.to_array()
	# lut[0]=[47,82,143,255]
	# lut[127]=[112,173,71,255]
	# lut[255]=[197,90,17,255]
	# pointsObject.module_manager.scalar_lut_manager.lut.table=lut



	# #now print lines connecting the points
	if(drawLines):



		splitArray=[]
		splitCounter=0
		firstSpecies=True



		thisStep=np.split(timeSteps[i],proteinNumber)

		for j in range(len(thisStep)):
			
			#remove contiguous overlapping units, so the lines don't disappear
			thisProteinLength=thisStep[j].shape
			thisProteinLength=thisProteinLength[0]


			lineNoAdjacentOverlapPoly=np.copy(thisStep[j])
			lineColors=np.arange(thisProteinLength)

			mask=[]
			for k in range(thisProteinLength):
				if np.array_equal(lineNoAdjacentOverlapPoly[k],lineNoAdjacentOverlapPoly[k-1]):
					mask.append(False)
				else:
					mask.append(True)

			lineNoAdjacentOverlapPoly=lineNoAdjacentOverlapPoly[mask]
			lineColors=lineColors[mask]


					#deal with crossings
			#print(lineNoAdjacentOverlapPoly)

			jumps=[]
			for k in range(1,lineNoAdjacentOverlapPoly.shape[0],1):
				if np.linalg.norm(lineNoAdjacentOverlapPoly[k]-lineNoAdjacentOverlapPoly[k-1])>3:
					jumps.append(k)

			if len(jumps)>0:
				lines=np.split(lineNoAdjacentOverlapPoly,jumps)
				for k in range(len(jumps)+1):                                                    #nb +1 if something breaks
					lineCoords=np.transpose(lines[k])
					lineObject=mlab.plot3d(lineCoords[0],lineCoords[1],lineCoords[2],color=(0,0,0))


			else:
				lines=lineNoAdjacentOverlapPoly
				lineCoords=np.transpose(lines)
				lineObject=mlab.plot3d(lineCoords[0],lineCoords[1],lineCoords[2],color=(0,0,0))



	if i==first:
		mlab.show(stop=True)
		v = mlab.view()
	



	mlab.view(*v)
	mlab.savefig("testSnap"+str(i)+".png")
	mlab.close()