//what is energy of initial state? then we already calculate DeltaE, so it's easy to track subsequent steps

using namespace std;

float getEnergy(vector<vector<vector<site>>> space, vector<int> size,float exchangeJ) 
{

	int bonds=0;
	float neighborPairs=0;
	vector<int> siteMotifs;

	for(int i=0;i<size[0];++i)
	{
		for(int j=0; j<size[1];++j)
		{
			for(int k=0; k<size[2]; ++k)
			{


				site thisSite=space[i][j][k];
				int thisOccupancy=thisSite.occupancy.size();

				//count bonds
				if(thisOccupancy==2)
				{
					bonds+=1;
				}

				//count neighbors
				int numNeighbors=0;
				for(int q=0;q<thisSite.neighbors.size();++q)
				{
					site* thisNeighborPointer=thisSite.neighbors[q];
					int thisNeighborOccupancy=thisNeighborPointer->occupancy.size();
					
					vector<int> occupancy=thisSite.occupancy;
					
					auto iterator=occupancy.begin();
						while(iterator !=occupancy.end())
						{
						    vector<int> myOccupancy=*iterator;
						    cout<<myOccupancy[0]<<" "<<myOccupancy[1]<<" "<<myOccupancy[2]<<", ";
						    ++iterator;
						}
					
				}

				//have to double neighbor pairs if the site is doubly occupied
				neighborPairs+=numNeighbors*thisOccupancy;


			}
		}
	}

	//remove double counting of neighbors
	neighborPairs=neighborPairs/2;
	float nonspecificEnergy=(-1)*neighborPairs*exchangeJ;

	float energy=(-1)*float(bonds)+nonspecificEnergy;


	return energy;

}
