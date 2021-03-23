///calculate number of neighbors in a configuration

using namespace std;

int getNeighborNumber(vector<polymer> polyVector)
{
	int neighborNumber=0;



	int polyLength=polyVector[0].length;
	for(int i=0;i<polyVector.size();++i)           //iterate over all occupied sites
	{
		for(int j=0;j<polyLength;++j)
		{
			for(int k=0;k<12;++k)
			{
				site* thisNeighbor=polyVector[i].chain[j]->neighbors[k];  //pointer to neighbor
				neighborNumber+=thisNeighbor->occupancy.size(); //add as many neighbors as there are occupants at that site
			}
		}
	}


	return int(neighborNumber/2); //we want the number of neighbor pairs
}