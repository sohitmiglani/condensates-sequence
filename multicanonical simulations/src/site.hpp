#include "randomNumberGenerator.hpp"

using namespace std;

class site{
public:
	vector<int> position;
	vector<site*> neighbors;
	set<vector<int>> occupancy; 
	//empty if unoccupied, add (polymer ID,motif,chain unit) pair for each occupancy (up to 2)
	//polymer ID: [0, n], motif: +1 or -1, chain unit: [0:L-1]

	void addNeighbor(site* point) //pass a pointer, which you add to list of pointers: never access site itself!
	{
		neighbors.push_back(point);
	}

	site()
	{
		position={-1,-1,-1};
	}

	void printSite()
	{
		cout<<"position:"<<endl;
		for (int i = 0; i < 3; ++i)
		{
			cout<<position[i]<<" ";
		}
		cout<<endl;
		cout<<"occupancy:"<<endl;

		auto iterator=occupancy.begin();
		while(iterator !=occupancy.end())
		{
			vector<int> myOccupancy=*iterator;
			cout<<myOccupancy[0]<<" "<<myOccupancy[1]<<" "<<myOccupancy[2]<<", ";
			++iterator;
		}
		cout<<endl;
	}

	bool contiguousBondQ()         //returns true if site has bond between two subunits that are contiguous on the same polymer 
	{
		if(occupancy.size()<2)
		{
			return false;
		}
		else
		{
			vector<vector<int>> occupancyVector;
			auto iterator=occupancy.begin();

			while(iterator !=occupancy.end())       //make a vector of the occupancy vectors, max 2
			{
			vector<int> myOccupancyVector=*iterator;
			occupancyVector.push_back(myOccupancyVector);
			++iterator;
			}
			if(occupancyVector[0][0]==occupancyVector[1][0]&&abs(occupancyVector[0][2]-occupancyVector[1][2])==1) 
				//contiguous if same polymer and |unit1-unit2|=1
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}

	site* randomNeighbor(int identityIncluded)      //return a random neighbor of pivot, if identity included=1, include the pivot
	{
		vector<site*> possibleMoves=this->neighbors;   //vector of pointers to adjacent sites
		if(identityIncluded==1)
		{
			possibleMoves.push_back(this);
		}
		int numNeighbors=possibleMoves.size();
		bool stillSearching=true;
		
		uniform_int_distribution<int> randomMoveDist(0,(numNeighbors-1));
		int randomMove=randomMoveDist(rng);
		// if(possibleMoves[randomMove]->occupancy.size()<2)
		// {
		// 	stillSearching=false;
		// 	return possibleMoves[randomMove];
		// }
		
		// if(stillSearching) //if you try to find a new site and can't, return nullptr
		// {
		// 	return nullptr;
		// }
		return possibleMoves[randomMove];

	}

//if the neighboring monomer is from the same protein, you only add 1/2 to avoid double-counting
	float getPairNumber(int queryID)
	{
		float pairNumber=0;
		vector<site*> neighborsToCheck=this->neighbors;   //vector of pointers to adjacent sites
		for(int i=0;i<neighborsToCheck.size();++i)
		{

			auto occIterator=neighborsToCheck[i]->occupancy.begin();
			while(occIterator!=neighborsToCheck[i]->occupancy.end())
			{
				vector<int> thisOccupancy=*occIterator;
				int neighborID=thisOccupancy[0];
				if(neighborID==queryID)
				{
					pairNumber+=0.5;
				}
				else
				{
					pairNumber+=1;
				}

				++occIterator;
			}

		}


		return pairNumber;


	}

		//if the neighbors couldn't possibly be moving, just use this
	int getPairNumber()
	{
		int pairNumber=0;
		vector<site*> neighborsToCheck=this->neighbors;   //vector of pointers to adjacent sites
		for(int i=0;i<neighborsToCheck.size();++i)
		{
			pairNumber+=neighborsToCheck[i]->occupancy.size();
		}

		return pairNumber;
	}

};