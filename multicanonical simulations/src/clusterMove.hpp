//this function performs moves on a whole cluster
//move type 5: translation
//move type 6: rotation
#include "randomNumberGenerator.hpp"

using namespace std;

set<int> getCluster(int,unordered_map<int,polymer>,auto);

tuple<bool,vector<int>,int,vector<tuple<vector<int>,site*, site*>>> clusterMove(int moveType,int startProtein,unordered_map<int,polymer> &polyMap,vector<vector<vector<site>>> &space,double exchangeJ)
{	
	bool successfulMove=false;
	int neighborPairChange=0;
	vector<int> bondChanges{0,0,0}; //we don't allow cluster moves to change bonds, but this will keep structure more consistent with polyMoves
	vector<tuple<vector<int>,site*, site*>> moveRecord; //keep track of where you've moved, so it's reversible

		//if nonspecific interactions are not included, we don't need to calculate neighborpairs
	bool nonspecificStatus=true;
	if(abs(exchangeJ)<=0.000001)
	{
		nonspecificStatus=false;
	}


	set<int> clusterSet=getCluster(startProtein,polyMap,space);        //get all the proteins connected to the randomly chosen one (set of IDs)

	if(moveType==5)             //translation
	{
		//because each point has neighbors in the same order, choosing a consistent neighbor is the same as choosing a direction
		int numNeighbors=space[0][0][0].neighbors.size();
		uniform_int_distribution<int> neighborDist(0, (numNeighbors-1));
		int randNeighbor=neighborDist(rng);

		//check that new points do not have proteins outside cluster
		//first, get set of occupants in target sites
		set<int> targetSiteOccupants;
		auto iterator=clusterSet.begin();            //iterate over the set of polymers in this cluster
		while(iterator !=clusterSet.end())           
		{
			vector<site*> thisChain=polyMap[*iterator].chain;      //vector of pointers to sites
			for(int i=0;i<thisChain.size();++i)
			{
				site* thisNeighbor=thisChain[i]->neighbors[randNeighbor];
				set<vector<int>> thisOccupancy=thisNeighbor->occupancy;
				
				//iterate over occupancy set
				auto occupancyIterator=thisOccupancy.begin();
				while(occupancyIterator!=thisOccupancy.end())
				{
					vector<int> thisOccupant=*occupancyIterator;
					targetSiteOccupants.insert(thisOccupant[0]);
					++occupancyIterator;
				}
			}

			++iterator;
		}

		//check that target occupants don't have anything outside cluster
		vector<int> unionVector;
		set_union(clusterSet.begin(), clusterSet.end(),targetSiteOccupants.begin(), targetSiteOccupants.end(),back_inserter(unionVector));

		if(unionVector.size()==clusterSet.size())
		//finally move everything: iterate over cluster AND chain
		{


			int pairsI;
			if(nonspecificStatus==true)
			{
				//first, calculate initial neighbor number
				pairsI=0;
				auto startNeighborIterator=clusterSet.begin();            //iterate over the set of polymers in this cluster
				while(startNeighborIterator !=clusterSet.end())           
				{
					vector<site*> thisMoveChain=polyMap[*startNeighborIterator].chain;      //vector of pointers to sites
					for(int i=0;i<thisMoveChain.size();++i)
					{
						pairsI+=thisMoveChain[i]->getPairNumber();
					}
					++startNeighborIterator;
				}
			}




			auto moveIterator=clusterSet.begin();            //iterate over the set of polymers in this cluster
			while(moveIterator !=clusterSet.end())           
			{
				vector<site*> thisMoveChain=polyMap[*moveIterator].chain;      //vector of pointers to sites
				for(int i=0;i<thisMoveChain.size();++i)
				{
					//get correct occupancy vector
					set<vector<int>> movingSiteOccupancy=thisMoveChain[i]->occupancy;
					int movingSeq;
					auto movingSiteOccupancyIterator=movingSiteOccupancy.begin();
					while(movingSiteOccupancyIterator!=movingSiteOccupancy.end())
					{
						vector<int> thisMovingSiteOccupancy=*movingSiteOccupancyIterator;
						if(thisMovingSiteOccupancy[0]==*moveIterator&&thisMovingSiteOccupancy[2]==i)
						{
							movingSeq=thisMovingSiteOccupancy[1];
						}
						++movingSiteOccupancyIterator;
					}

					vector<int> occupancyVector={*moveIterator,movingSeq,i};
					site *thisTargetSite=thisMoveChain[i]->neighbors[randNeighbor];
					polyMap[*moveIterator].moveUnit(occupancyVector,thisMoveChain[i],thisTargetSite,moveRecord);
				}

				++moveIterator;
			}


			if(nonspecificStatus==true)
			{	
				//now, calculate final neighbor number
				int pairsF=0;
				auto endNeighborIterator=clusterSet.begin();            //iterate over the set of polymers in this cluster
				while(endNeighborIterator !=clusterSet.end())           
				{
					vector<site*> thisMoveChain=polyMap[*endNeighborIterator].chain;      //vector of pointers to sites
					for(int i=0;i<thisMoveChain.size();++i)
					{
						pairsF+=thisMoveChain[i]->getPairNumber();
					}
					++endNeighborIterator;
				}


				neighborPairChange=pairsF-pairsI;
			}
			
			successfulMove=true;
		}

	}

	return make_tuple(successfulMove,bondChanges,neighborPairChange,moveRecord);
}