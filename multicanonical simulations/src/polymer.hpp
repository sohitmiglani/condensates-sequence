#include "randomNumberGenerator.hpp"
using namespace std;


//Latest code: main initializes polymer with default constructor, then uses void method "initialize" to actually place it

vector<int> getBondInfo(site*);
vector<int> vectorAddition(vector<int>, vector<int>, int);  

class polymer{
public:
	vector<site*> chain;
	int polymerID;
	int length;
	int dim;
	int numNeighbors;
	vector<int> sequence;

	//default constructor
	polymer()
	{

	}



	polymer(vector<vector<int>> givenChain, auto polyInfo, auto &space, int id)       //initialize a given polymer, for testing
	{
		polymerID=id;
		length=get<0>(polyInfo);
		sequence=get<1>(polyInfo);

		site *lastPoint;
		dim=givenChain[1].size();

		for(int i=0;i<length;++i)
		{
			int thisX=givenChain[i][0];
			int thisY=givenChain[i][1];
			int thisZ=givenChain[i][2];

			vector<int> myOccupancy={polymerID,sequence[i],i};
			space[thisX][thisY][thisZ].occupancy.insert(myOccupancy);        //now occuppied
			
			lastPoint=&space[thisX][thisY][thisZ];        //pointer to previous units
			chain.push_back(lastPoint);  

		}
	}

//initialize a polymer that is identical to the given polymer, but in a new space, for MC testing
	polymer(polymer originalPoly, auto &newSpace)
	{
		//copy all the simple class properties
		polymerID=originalPoly.polymerID;
		length=originalPoly.length;
		dim=originalPoly.dim;
		numNeighbors=originalPoly.numNeighbors;
		sequence=originalPoly.sequence;

		//set the chain and site occupancies in the *new* space, based on the old space
		for(int i=0;i<length;++i)
		{
			vector<int> myOccupancy={polymerID,sequence[i],i};
			vector <int> unitPosition=originalPoly.chain[i]->position;   //position in old space, to access correct point in new space

			newSpace[unitPosition[0]][unitPosition[1]][unitPosition[2]].occupancy.insert(myOccupancy); //insert occupancy in new space, not a pointer

			site* newSpaceSite=&newSpace[unitPosition[0]][unitPosition[1]][unitPosition[2]];
			chain.push_back(newSpaceSite);                         //add new space site pointer to chain
			
		}


	}



//initialize polymer with a given set of coordinates
	//space not necessary because we have access to site pointers directly
polymer(auto polyInfo, vector<vector<int>> coords, int id,auto &space)
{
	polymerID=id;
	length=get<0>(polyInfo);
	sequence=get<1>(polyInfo);

	dim=coords[0].size();

	for(int i=0;i<length;i++)
	{


		vector<int> thisCoord=coords[i];

		vector<int> myOccupancy={polymerID,sequence[i],i};
		space[thisCoord[0]][thisCoord[1]][thisCoord[2]].occupancy.insert(myOccupancy);
		chain.push_back(&space[thisCoord[0]][thisCoord[1]][thisCoord[2]]);

		cout<<"reversed chain:"<<chain[i]->position[0]<<" "<<chain[i]->position[1]<<" "<<chain[i]->position[2]<<endl;
	}

}

//places subunits at the neighbors of a given site, testing connectivity
	polymer(vector<int> center,auto latticeInfo,auto &space) 
	{
		numNeighbors=get<2>(latticeInfo);
		length=numNeighbors+1;
		dim=3;
		site *lastPoint;

		vector<int> centerOccupancy={-1,-1,0};
		space[center[0]][center[1]][center[2]].occupancy.insert(centerOccupancy);  //a site

		lastPoint=&space[center[0]][center[1]][center[2]];                //pointer to site
		chain.push_back(lastPoint);  

		vector<site*> possibleNeighbors=lastPoint->neighbors;         //this is a vector of site pointers
		site *thisNeighbor;								          //points to specific neighbor

		for(int i=0;i<numNeighbors;++i)
		{
			thisNeighbor=possibleNeighbors[i];

			vector<int> myOccupancy={i,1,i};
			thisNeighbor->occupancy.insert(myOccupancy);
			
			chain.push_back(thisNeighbor);
		}
	}


	void printChain()
	{
		for(int i=0;i<length;++i)
		{	
			site *thisSite=chain[i];
			for (int j = 0; j < dim; ++j)
			{
				int here=thisSite->position[j];
				cout<<here<<" ";
			}	
			cout<<endl;
		}
	}

	void writePolyChain(string filename)
	{
		
		ofstream polyOutput(filename,fstream::app);
		if (!polyOutput)
		{
	       	cout << "Unable to open polymer output file";
	   		exit(1); // terminate with error
	   	};
	 	

	    for(int i=0;i<length;++i)
	    {	
			site *thisSite=chain[i];
			for (int j = 0; j < dim; ++j)
			{
				int here=thisSite->position[j];
				polyOutput << here <<" ";
			}	
			polyOutput << endl;
	    }
	   
	}

	tuple<bool,vector<int>,int,vector<tuple<vector<int>,site*, site*>>> movePolymer(int moveType, double exchangeJ);
	//implemented in polyMoves.hpp. type 0 is end-move
	//arguments: move success, bond changes, neighbor changes, move record


//moves and makes record
	void moveUnit(vector<int> movingUnit, site *oldSite, site *newSite, vector<tuple<vector<int>,site*, site*>> &passedMoveRecord)
	{

		if(oldSite!=newSite)
		{
		oldSite->occupancy.erase(movingUnit);
		newSite->occupancy.insert(movingUnit);
		int unitNumber=movingUnit[2];
		chain[unitNumber]=newSite;

		//update the moving record. now we just need to pass this to each moveUnit
		passedMoveRecord.push_back(make_tuple(movingUnit,oldSite,newSite));
		}


	}

//makes no record, just moves
		void moveUnit(vector<int> movingUnit, site *oldSite, site *newSite)
	{

		if(oldSite!=newSite)
		{
		oldSite->occupancy.erase(movingUnit);
		newSite->occupancy.insert(movingUnit);
		int unitNumber=movingUnit[2];
		chain[unitNumber]=newSite;
		}


	}

//figure out the change in energy of inserting this polymer into the system. Put another way, how many bonds and nearest neighbor interactions does this polymer participate in?
	tuple<vector<int>,int> insertionEnergy(double exchangeJ)
	{	

		vector<int> bonds{0,0,0}; //delta +-, then delta ++/--, then delta triples
		float neighborPairs=0;
		int otherNeighbors=0;
		double selfNeighbors=0;

			//if nonspecific interactions are not included, we don't need to calculate neighborpairs
		bool nonspecificStatus=true;
		if(abs(exchangeJ)<=0.000001)
		{
			nonspecificStatus=false;
		}

		//add segments to set to get unique sites
		set<site*> uniqueSites;
		for(int i=0;i<length;++i)
		{
			uniqueSites.insert(chain[i]);
		}

		auto siteIterator=uniqueSites.begin();
		while(siteIterator!=uniqueSites.end())
		{
			site *thisSite=*siteIterator; //pointer to site, after getting value of set iterator
			bonds=vectorAddition(bonds,getBondInfo(thisSite),1);

			if(nonspecificStatus==true) //only need neighbor pairs if J=/=0
			{
				int numberSelfMonomers=1; //are there two monomers from this protein at this site?
				if(thisSite->occupancy.size()==2)
				{
					auto myOccupancyIt=thisSite->occupancy.begin();
					vector<int> selfOccupantIDs;
					while(myOccupancyIt!=thisSite->occupancy.end())
					{
						vector<int> myOccupancy=*myOccupancyIt;
						selfOccupantIDs.push_back(myOccupancy[0]);
						++myOccupancyIt;
					}
					if(selfOccupantIDs[0]==selfOccupantIDs[1])
					{
						numberSelfMonomers=2;
					}
				}

				float deltaNeighbors=thisSite->getPairNumber(polymerID);
				neighborPairs+=numberSelfMonomers*deltaNeighbors;

	

			}

			siteIterator++;
		}


		return make_tuple(bonds,int(neighborPairs));
	}

//remove occupancy vectors from sites
	void deletePoly()
	{	
		//we only go to chain.size() because we may have to delete an aborted polymer from a failed insertion
		for(int i=0;i<chain.size();++i)
		{
			site *thisSite=chain[i];
			vector<int> thisOccupant{polymerID,sequence[i],i};
			thisSite->occupancy.erase(thisOccupant);
		}
	}


//remove occupancy vectors while calculating weight to construct that polymer
	double deletePoly(bool getWeight)
	{
		double weight=1;

		for(int i=(length-1);i>-1;--i)
		{

			site *thisSite=chain[i];
			vector<int> thisOccupant{polymerID,sequence[i],i};

			//get weight if this monomer were to grow a new monomer
			if(i<(length-1))
			{
				//we are trying to place the next sequence, not the current one
				int nextMotif=sequence[i+1];

				vector<site*> allNeighbors=thisSite->neighbors;
				allNeighbors.push_back(thisSite);                    //include the current site so we can deal with contiguous bonds
				vector<site*> validNeighbors;

				for(int j=0;j<allNeighbors.size();j++)
				{
					vector<int> thisOccupancyVector;
					if(allNeighbors[j]->occupancy.size()==1)
					{
						auto thisOccupancyIterator=allNeighbors[j]->occupancy.begin();
						thisOccupancyVector=*thisOccupancyIterator;
					}
					if(allNeighbors[j]->occupancy.size()==0||allNeighbors[j]->occupancy.size()==1&&thisOccupancyVector[1]!=nextMotif)
					{

						validNeighbors.push_back(allNeighbors[j]);
					}

				}

				weight=weight*validNeighbors.size();

			}


			//delete this monomer
			thisSite->occupancy.erase(thisOccupant);
		}

	return weight;
	}


	//initialize polymer in open configuration to avoid loops and geometrical frustration. necessary to do this after construction to
	//avoid confusing make_pair with template matching when we store polymers in unordered_map
	void initialize(auto latticeInfo, auto polyInfo, auto &space, int ID, string openID)
	{
		polymerID=ID;
		length=get<0>(polyInfo);
		sequence=get<1>(polyInfo);

		auto size=get<1>(latticeInfo);
		numNeighbors=get<2>(latticeInfo);
		dim=size.size();

		site *lastPoint;

		//choose and place random start point (requires modification if space is not cubic)
		bool startPlaced=false;
		uniform_int_distribution<int> startDist(0, (size[0]-1));

		while(!startPlaced)
		{
			int startX=startDist(rng);
			int startY=startDist(rng);
			int startZ=startDist(rng);
			//we allow placement on empty site, or site with opposite motif

			vector<int> occupancyVector;
			if(space[startX][startY][startZ].occupancy.size()==1)
			{
			auto occIterator=space[startX][startY][startZ].occupancy.begin();
			occupancyVector=*occIterator;
			}
			if(space[startX][startY][startZ].occupancy.size()==0||space[startX][startY][startZ].occupancy.size()==1&&occupancyVector[1]!=sequence[0])
			{
				vector<int> myOccupancy={polymerID,sequence[0],0};
				space[startX][startY][startZ].occupancy.insert(myOccupancy);        //now occupied
				lastPoint=&space[startX][startY][startZ];        //pointer to previous units
				chain.push_back(lastPoint); 
				startPlaced=true;
			}
		}

		//choose random direction
		uniform_int_distribution<int> directionDist(0,(numNeighbors-1));
		int direction=directionDist(rng);


		for(int i=1;i<length;++i)
		{

			vector<int> thisOccupancy={polymerID,sequence[i],i};
			site *nextSite;
			bool foundNextSite=false;
			while(!foundNextSite) //make sure the next point is unoccupied or occupied with only 1 compatible motif
			{
				nextSite=lastPoint->neighbors[direction];        //pointer to next site

				vector<int> occupancyVector;
				if(nextSite->occupancy.size()==1)
				{
				auto occIterator=nextSite->occupancy.begin();  //iterator pointing to first occupancy vector 
				occupancyVector=*occIterator;
				}
				if(nextSite->occupancy.size()==0||nextSite->occupancy.size()==1&&occupancyVector[1]!=sequence[i])
				{
					foundNextSite=true;
				}
				else //try a new direction if occupied
				{
					direction=directionDist(rng);
				}

			}

			nextSite->occupancy.insert(thisOccupancy); //now occupied
			lastPoint=nextSite;
			chain.push_back(lastPoint); 
			
		}
	}

	//a method to initialize with specific coordinates after inserting default polymer into map
	void setConfig(auto polyInfo, vector<site*> coords, int id)
	{
		polymerID=id;
		length=get<0>(polyInfo);
		sequence=get<1>(polyInfo);

		dim=coords[0]->position.size();

		for(int i=0;i<length;i++)
		{

			vector<int> myOccupancy={polymerID,sequence[i],i};
			coords[i]->occupancy.insert(myOccupancy);
			chain.push_back(coords[i]);
		}
	}


//initialize polymer with valid walk and return the Rosenbluth weight, along with a success boolean
	//"valid" means no illegal overlaps (triples or ++/--)
	tuple<bool,double> growBiased(auto polyInfo,auto &space, int ID)
	{
		bool trialSuccess=true;
		double weight=1;


		length=get<0>(polyInfo);
		sequence=get<1>(polyInfo);
		polymerID=ID;

		dim=space[0][0][0].position.size();


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
			chain.push_back(lastPoint); 
			vector<int> myOccupancy={polymerID,sequence[0],0};
			lastPoint->occupancy.insert(myOccupancy);

		}
		else //failed to place head? give up, so we don't have to calculate free fracton every step for weight
		{
			trialSuccess=false;
		}

		//now try the body***********************************
		if(trialSuccess==true)
		{
			int i=1;
			while(trialSuccess==true && i<length)
			{
				vector<site*> allNeighbors=lastPoint->neighbors;
				allNeighbors.push_back(lastPoint);
				vector<site*> validNeighbors;

				//check nearest neighbors for valid sites
				for(int j=0;j<allNeighbors.size();j++)
				{
					vector<int> thisOccupancyVector;
					if(allNeighbors[j]->occupancy.size()==1)
					{
					auto thisOccupancyIterator=allNeighbors[j]->occupancy.begin();
					thisOccupancyVector=*thisOccupancyIterator;
					}
					if(allNeighbors[j]->occupancy.size()==0||allNeighbors[j]->occupancy.size()==1&&thisOccupancyVector[1]!=sequence[i])
					{
						validNeighbors.push_back(allNeighbors[j]);
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

					chain.push_back(validNeighbors[neighborIndex]);
					lastPoint=validNeighbors[neighborIndex];
					vector<int> myOccupancy={polymerID,sequence[i],i};
					lastPoint->occupancy.insert(myOccupancy);



					weight=weight*validNeighbors.size();


					i++;
				}
				else //no valid neighbors? end the trial as a failure
				{
					trialSuccess=false;
				}
			}
		}


		return make_tuple(trialSuccess,weight);
		}
};