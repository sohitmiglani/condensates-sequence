///choose a random head site and try to construct a self=avoiding walk of length L-1. Along the way, calculate Rosenbluth weights.
//return the a vector of site pointers, as well as the move success and the Rosenbluth weight=Prod_i free sites_i
#include "randomNumberGenerator.hpp"
using namespace std;

tuple<bool,vector<site*>,float> getInsertionConfig(vector<vector<vector<site>>> &space, auto polyInfo) 
{

	vector<site*> newConfig; //the vector of site pointers we're proposing
	bool trialSuccess=true;
	float weight=1;


	int length=get<0>(polyInfo);
	vector<int> sequence=get<1>(polyInfo);

	site *lastPoint;

	//choose random head site****************************************************************
	int linearDimension=space[0].size();
	uniform_int_distribution<int> startDist(0, (linearDimension-1));

	int headX=startDist(rng);
	int headY=startDist(rng);
	int headZ=startDist(rng);

	//we allow placement on empty site, or site with opposite motif

	vector<int> occupancyVector;
	if(space[headX][headY][headZ].occupancy.size()==1)
	{
	auto occIterator=space[headX][headY][headZ].occupancy.begin();
	occupancyVector=*occIterator;
	}
	if(space[headX][headY][headZ].occupancy.size()==0||space[headX][headY][headZ].occupancy.size()==1&&occupancyVector[1]!=sequence[0])
	{

		lastPoint=&space[headX][headY][headZ];        //pointer to previous units
		newConfig.push_back(lastPoint); 

	}
	else //failed to place head? give up, so we don't have to calculate free fracton every step for weight
	{
		trialSuccess=false;
	}

	//now try the body***********************************
	if(trialSuccess==true)
	{
		cout<<"successful head insertion"<<endl;
		int i=1;
		while(trialSuccess==true && i<length)
		{
			vector<site*> allNeighbors=lastPoint->neighbors;
			vector<site*> validNeighbors;

			//check nearest neighbors for valid sites
			for(int j=0;j<allNeighbors.size();j++)
			{
				vector<int> thisOccupancyVector;
				cout<<"this neighbor occupancy: "<<allNeighbors[j]->occupancy.size()<<endl;
				if(allNeighbors[j]->occupancy.size()==1)
				{
				auto thisOccupancyIterator=allNeighbors[j]->occupancy.begin();
				thisOccupancyVector=*thisOccupancyIterator;
				}
				if(allNeighbors[j]->occupancy.size()==0||allNeighbors[j]->occupancy.size()==1&&thisOccupancyVector[1]!=sequence[i])
				{

					cout<<"found valid neighbor"<<endl;
					validNeighbors.push_back(allNeighbors[j]);
					if(allNeighbors[j]->occupancy.size()==1)
					{
					cout<<"occupancy: "<<thisOccupancyVector[1]<<", seq: "<<sequence[i]<<endl;
					}
				}


			}
			//if there are compatible neighbors, randomly choose one and compute weights
			if(validNeighbors.size()>0)
			{
				int neighborIndex;
				if(validNeighbors.size()==1)
				{
					neighborIndex=0;
				}
				else
				{
					uniform_int_distribution<int> neighborDist(0, (validNeighbors.size()-1));
					neighborIndex=neighborDist(rng);
				}

				newConfig.push_back(validNeighbors[neighborIndex]);
				lastPoint=validNeighbors[neighborIndex];
				weight=weight*validNeighbors.size();


				i++;
			}
			else //no valid neighbors? end the trial as a failure
			{
				trialSuccess=false;
			}
		}
	}


	return make_tuple(trialSuccess,newConfig,weight);



}