#include "randomNumberGenerator.hpp"

//no nonspecific interactions!
using namespace std;
//this code does not check for defects, assumes energy penalties will take care of everything
//type 0: end move
//type 1: reptation
//type 2: contraction
//type 3: expansion
//type 4: corner moves
//we're keeping track of neighbor pair changes, bond changes

vector<int> getBondInfo(site*);
vector<int> vectorAddition(vector<int>, vector<int>, int);  

tuple<bool,vector<int>,int,vector<tuple<vector<int>,site*, site*>>> polymer::movePolymer(int moveType,double exchangeJ){
	bool successfulMove=false;
	vector<int> bondChanges{0,0,0}; //delta +-, then delta ++/--, then delta triples
	int neighborPairChange=0;
	vector<tuple<vector<int>,site*, site*>> moveRecord; //keep track of where you've moved, so it's reversible

	//if nonspecific interactions are not included, we don't need to calculate neighborpairs
	bool nonspecificStatus=true;
	if(abs(exchangeJ)<=0.000001)
	{
		nonspecificStatus=false;
	}

	switch(moveType)
	{
		case 0:          //end move
		{
			uniform_int_distribution<int> endDist(0, 1);
			int end=endDist(rng);//randomly choose end, picks 0 or 1
			site* movePoint;
			site* movePointPivot;
			if (end==0)
			{
				movePoint=chain[0];
				movePointPivot=chain[1];
			}
			else
			{
				end=length-1;
				movePoint=chain[end];
				movePointPivot=chain[end-1];
			}

			vector<int> myOccupancy={polymerID,sequence[end],end}; //this is the moving polymer unit, which we'll use to update the site occupancy

			site* randomMove=movePointPivot->randomNeighbor(0);    //pick a random neighbor with occupancy 0 or 1
			if(randomMove!=nullptr&&!(movePoint->contiguousBondQ()))
			{
				int pairsI;
				if(nonspecificStatus==true)
				{
					pairsI=movePoint->getPairNumber();
				}
				vector<int> startBonds=vectorAddition(getBondInfo(movePoint),getBondInfo(randomMove),1);

				moveUnit(myOccupancy,movePoint,randomMove,moveRecord);         //move the unit
				successfulMove=true;

				if(nonspecificStatus==true)
				{
				//find changes in neighbors
				int pairsF=randomMove->getPairNumber();
				neighborPairChange=pairsF-pairsI; //keep track of pair change for non-specific interaction
				}

				//bond changes
				vector<int> endBonds=vectorAddition(getBondInfo(movePoint),getBondInfo(randomMove),1);
				
				bondChanges=vectorAddition(endBonds,startBonds,-1);



			}
		}
			break;

		case 1:         //reptation
		{
			uniform_int_distribution<int> endDist(0, 1);
			int end=endDist(rng);//randomly choose end, picks 0 or 1
			
			if(end!=0)
			{end=length-1;}
			
			site* randomMove=chain[end]->randomNeighbor(0);          //get random neighbor of leading unit
			if(randomMove!=nullptr)
			{
				//first, check bonds at current chain and move site. And pairs at current site
				vector<int> startBonds={0,0,0};
				float pairsI=0;
				site* originallyOccupied; //we need a pointer to the site that may no longer be occupied at the end, for bond checking
				if(end==0)
				{
					originallyOccupied=chain[(length-1)];
				}
				else
				{
					originallyOccupied=chain[0];
				}

				set<site*> uniqueSites;  //set of pointers to sites, without repeats from contiguous bonds
				for(int i=0;i<length;++i)
				{
					uniqueSites.insert(chain[i]);
				}

				auto uniqueIterator=uniqueSites.begin();
				while(uniqueIterator!=uniqueSites.end())
				{
					site* thisSite=*uniqueIterator;
					startBonds=vectorAddition(startBonds,getBondInfo(thisSite),1);


					if(nonspecificStatus==true)
					{
						float pairChange=thisSite->getPairNumber(polymerID);
						if(thisSite->occupancy.size()==2)                               //if the moving protein has 2 bonds here, count both
						{
							vector<vector<int>> occupancyVector(thisSite->occupancy.begin(), thisSite->occupancy.end());  //initialize vector of occupancy vectors
							if(occupancyVector[0][0]==polymerID&&occupancyVector[1][0]==polymerID)
							{
								pairChange=pairChange*2;                        //if both the monomers could have moved, then check both their nonspecific interactions
							}
						}
						
						//pairChange=pairChange*(thisSite->occupancy.size()); //we allow each monomer to contribute
						pairsI+=pairChange;
					}

					++uniqueIterator;
				}

				if(uniqueSites.find(randomMove)==uniqueSites.end()) //if the place we're moving to is not our protein, add it too
				{
					startBonds=vectorAddition(startBonds,getBondInfo(randomMove),1);
				}

				//check for contiguous bonds before moving anything
				vector<bool> contigBondRecord;
				bool thisRecord;
				for(int i=0;i<length;++i)
				{	thisRecord=chain[i]->contiguousBondQ();
					contigBondRecord.push_back(thisRecord);}
				
				 //move body of the chain
				if(end==0)  //[0] end is moving, so we start at [length-1]
				{
					for(int i=length-1;i>=2;--i) 
					{
						if(!contigBondRecord[i])
						{
						vector <int> myOccupancy={polymerID,sequence[i],i};

						moveUnit(myOccupancy,chain[i],chain[i-1],moveRecord);
						}
						else
						{
						vector <int> myOccupancy1={polymerID,sequence[i],i};
						vector <int> myOccupancy2={polymerID,sequence[i-1],i-1};
						moveUnit(myOccupancy1,chain[i],chain[i-2],moveRecord);
						moveUnit(myOccupancy2,chain[i-1],chain[i-2],moveRecord);
						--i;

						}
					}
				}
				else               //[length-1] end is moving, so we start at [0]
				{
					for(int i=0;i<length-2;++i)
					{
						if(!contigBondRecord[i])
						{
						vector <int> myOccupancy={polymerID,sequence[i],i};
						moveUnit(myOccupancy,chain[i],chain[i+1],moveRecord);
						}
						else
						{
						vector <int> myOccupancy1={polymerID,sequence[i],i};
						vector <int> myOccupancy2={polymerID,sequence[i+1],i+1};
						moveUnit(myOccupancy1,chain[i],chain[i+2],moveRecord);
						moveUnit(myOccupancy2,chain[i+1],chain[i+2],moveRecord);
						++i;
						}
					 }
				}
				//move head and its neighbor. note the special cases for if the head is contiguous
				if(!contigBondRecord[end])
				{
					if(end==0)
					{
						vector <int> myOccupancy={polymerID,sequence[1],1};
						moveUnit(myOccupancy,chain[1],chain[end],moveRecord);
					}
					else
					{
						int thisIndex=length-2;
						vector <int> myOccupancy={polymerID,sequence[thisIndex],thisIndex};
						moveUnit(myOccupancy,chain[thisIndex],chain[end],moveRecord);
					}
					vector<int> leadingOccupancy={polymerID,sequence[end],end}; //the moving end
					moveUnit(leadingOccupancy,chain[end],randomMove,moveRecord);
					
				}
				else
				{
					vector<int> leadingOccupancy={polymerID,sequence[end],end}; //the moving end
					moveUnit(leadingOccupancy,chain[end],randomMove,moveRecord);
					if(end==0)
					{
						vector <int> myOccupancy={polymerID,sequence[1],1};
						moveUnit(myOccupancy,chain[1],randomMove,moveRecord);
					}
					else
					{
						int thisIndex=length-2;
						vector <int> myOccupancy={polymerID,sequence[thisIndex],thisIndex};
						moveUnit(myOccupancy,chain[thisIndex],randomMove,moveRecord);
					}
				}

				//now figure out end bonds, pairs and update
				vector<int> endBonds={0,0,0};
				float pairsF=0;
	
				set<site*> uniqueEndSites;  //set of pointers to sites, without repeats from contiguous bonds
				for(int i=0;i<length;++i)
				{
					uniqueEndSites.insert(chain[i]);
				}

				auto uniqueEndIterator=uniqueEndSites.begin();
				while(uniqueEndIterator!=uniqueEndSites.end())
				{
					site* thisSite=*uniqueEndIterator;
					endBonds=vectorAddition(endBonds,getBondInfo(thisSite),1);

					if(nonspecificStatus==true)
					{
						float pairChange=thisSite->getPairNumber(polymerID);
						if(thisSite->occupancy.size()==2)                               //if the moving protein has 2 bonds here, count both
						{
							vector<vector<int>> occupancyVector(thisSite->occupancy.begin(), thisSite->occupancy.end());  //initialize vector of occupancy vectors
							if(occupancyVector[0][0]==polymerID&&occupancyVector[1][0]==polymerID)
							{
								pairChange=pairChange*2;                        //if both the monomers could have moved, then check both their nonspecific interactions
							}
						}
						//pairChange=pairChange*(thisSite->occupancy.size()); //we allow each monomer to contribute
						pairsF+=pairChange;
					}

					++uniqueEndIterator;
				}

				if(uniqueEndSites.find(originallyOccupied)==uniqueEndSites.end()) //if we've abandoned a site, add it too
				{
					endBonds=vectorAddition(endBonds,getBondInfo(originallyOccupied),1);
				}

				neighborPairChange=int(pairsF-pairsI);
				bondChanges=vectorAddition(endBonds,startBonds,-1);



				successfulMove=true;

			}
		}
		break;

		case 2:   //contraction move: allow two inner units to form contiguous bond
		{
			//choose direction, site
			uniform_int_distribution<int> directionDist(0, 1);
			int neighborDirection=directionDist(rng);//randomly choose end, picks 0 or 1
			if(neighborDirection==0) 
				{neighborDirection=-1;}                       //direction -1 is towards [0], direction 1 is towards [l-1]

			site* interiorSite=nullptr;
			int interiorSiteIndex;
			

			if(chain.size()==2) //no interior sites, really
			{
			uniform_int_distribution<int> interiorSiteDist(0,1);
			interiorSiteIndex= interiorSiteDist(rng); //randomly choose interior site

			}
			else
			{
			uniform_int_distribution<int> interiorSiteDist(1, (length-2));
			interiorSiteIndex= interiorSiteDist(rng); //randomly choose interior site

			}

			int neighborIndex=interiorSiteIndex+neighborDirection;


			if(neighborIndex>=0&&neighborIndex<=(length-1)&&!chain[interiorSiteIndex]->contiguousBondQ()&&!chain[neighborIndex]->contiguousBondQ()) //first two cases are for L=2
			{
				interiorSite=chain[interiorSiteIndex];
			}
				
			

			//check for contiguous bonds before moving anything
			vector<bool> contigBondRecord;
			bool thisRecord;
			for(int i=0;i<length;++i)
			{	thisRecord=chain[i]->contiguousBondQ();
				contigBondRecord.push_back(thisRecord);}


			if(interiorSite!=nullptr)     //proceed if you've chosen an eligible site
			{


				//first, check bonds at current chain and pairs at current site
				vector<int> startBonds{0,0,0};
				float pairsI=0;
				set<site*> uniqueStartSites;  //set of pointers to sites, without repeats from contiguous bonds
				for(int i=0;i<length;++i)
				{
					uniqueStartSites.insert(chain[i]);
				}
				auto uniqueStartIterator=uniqueStartSites.begin();
				while(uniqueStartIterator!=uniqueStartSites.end())
				{
					site* thisSite=*uniqueStartIterator;
					startBonds=vectorAddition(startBonds,getBondInfo(thisSite),1);

					if(nonspecificStatus==true)
					{
						float pairChange=thisSite->getPairNumber(polymerID);

						if(thisSite->occupancy.size()==2)                               //if the moving protein has 2 bonds here, count both
						{
							vector<vector<int>> occupancyVector(thisSite->occupancy.begin(), thisSite->occupancy.end());  //initialize vector of occupancy vectors
							if(occupancyVector[0][0]==polymerID&&occupancyVector[1][0]==polymerID)
							{
								pairChange=pairChange*2;                        //if both the monomers could have moved, then check both their nonspecific interactions
							}
						}
						//pairChange=pairChange*(thisSite->occupancy.size()); //we allow each monomer to contribute
						pairsI+=pairChange;
					}

					++uniqueStartIterator;
				}
				site* originallyOccupied; //in case we abandon a site completely
				if(neighborDirection==-1)
				{
					originallyOccupied=chain[0];
				}
				else
				{
					originallyOccupied=chain[(length-1)];
				}


				if(neighborDirection==-1) //move units from 0 to 0+1, up to interior site
				{
					for(int i=0;i<interiorSiteIndex;++i)
					{
						if(!contigBondRecord[i])
						{
							vector<int>myOccupancy={polymerID,sequence[i],i};
							moveUnit(myOccupancy,chain[i],chain[i+1],moveRecord);
						}
						else
						{
							vector <int> myOccupancy1={polymerID,sequence[i],i};
							vector <int> myOccupancy2={polymerID,sequence[i+1],i+1};
							moveUnit(myOccupancy1,chain[i],chain[i+2],moveRecord);
							moveUnit(myOccupancy2,chain[i+1],chain[i+2],moveRecord);
							++i;
						}
					}
				}
				else      //neighbor is going towards [l-1], so we move stuff down
				{
					for(int i=length-1;i>interiorSiteIndex;--i)
					{
						if(!contigBondRecord[i])
						{
							vector<int>myOccupancy={polymerID,sequence[i],i};
							moveUnit(myOccupancy,chain[i],chain[i-1],moveRecord);
						}
						else
						{
							vector <int> myOccupancy1={polymerID,sequence[i],i};
							vector <int> myOccupancy2={polymerID,sequence[i-1],i-1};
							moveUnit(myOccupancy1,chain[i],chain[i-2],moveRecord);
							moveUnit(myOccupancy2,chain[i-1],chain[i-2],moveRecord);
							--i;
						}
					}
				}

				//now check bonds and pairs at end
				vector<int> endBonds={0,0,0};
				float pairsF=0;
				set<site*> uniqueEndSites;  //set of pointers to sites, without repeats from contiguous bonds
				for(int i=0;i<length;++i)
				{
					uniqueEndSites.insert(chain[i]);
				}

				auto uniqueEndIterator=uniqueEndSites.begin();
				while(uniqueEndIterator!=uniqueEndSites.end())
				{
					site* thisSite=*uniqueEndIterator;
					endBonds=vectorAddition(endBonds,getBondInfo(thisSite),1);

					if(nonspecificStatus==true)
					{
						float pairChange=thisSite->getPairNumber(polymerID);

						if(thisSite->occupancy.size()==2)                               //if the moving protein has 2 bonds here, count both
						{
							vector<vector<int>> occupancyVector(thisSite->occupancy.begin(), thisSite->occupancy.end());  //initialize vector of occupancy vectors
							if(occupancyVector[0][0]==polymerID&&occupancyVector[1][0]==polymerID)
							{
								pairChange=pairChange*2;                        //if both the monomers could have moved, then check both their nonspecific interactions
							}
						}

						//pairChange=pairChange*(thisSite->occupancy.size()); //we allow each monomer to contribute
						pairsF+=pairChange;
					}

					++uniqueEndIterator;
				}

				if(uniqueEndSites.find(originallyOccupied)==uniqueEndSites.end()) //if we've abandoned a site, add it too
				{
					endBonds=vectorAddition(endBonds,getBondInfo(originallyOccupied),1);
				}

				neighborPairChange=int(pairsF-pairsI);
				bondChanges=vectorAddition(endBonds,startBonds,-1);
				successfulMove=true;
			}

		}
		break;

		case 3:      //expansion move: allow contiguous bond to be broken
		{
			//choose direction, site
			uniform_int_distribution<int> directionDist(0, 1);
			int neighborDirection=directionDist(rng);//randomly choose end, picks 0 or 1
			if(neighborDirection==0) 
				{neighborDirection=-1;}                       //direction -1 is expanding towards [0], direction 1 is towards [l-1]

			site* interiorSite=nullptr;
			int interiorSiteIndex;

			site* newHeadSite=nullptr;

						
			if(chain.size()==2) //no interior sites, really
			{
			uniform_int_distribution<int> interiorSiteDist(0,1);
			interiorSiteIndex= interiorSiteDist(rng); //randomly choose interior site

			}
			else
			{
			uniform_int_distribution<int> interiorSiteDist(1, (length-2));
			interiorSiteIndex= interiorSiteDist(rng); //randomly choose interior site

			}



			//choose site for expanding head
			if(neighborDirection==-1)
			{
				newHeadSite=chain[0]->randomNeighbor(0);    
			}
			else if(neighborDirection==1)
			{
				int localEnd=length-1;
				newHeadSite=chain[localEnd]->randomNeighbor(0);
			}

			//now check condition to forbid asymmetric head expansion
			// bool headConditionSatisfied=true;
			// bool headBond1=interiorSiteIndex==1&&chain[0]->contiguousBondQ();
			// bool headBond2=interiorSiteIndex==(length-2)&&chain[(length-1)]->contiguousBondQ();
			// bool direction1Condition=headBond1&&neighborDirection==1;
			// bool direction2Condition=headBond2&&neighborDirection==-1;
			// if(direction1Condition||direction2Condition)
			// 	{headConditionSatisfied=false;}

			//we want to check condition that we only expand if the neighbor direction and expansion direction match
			//find the chain index contiguous to the neighbor you're moving
			bool fullExpansionCondition=false;
			if(chain[interiorSiteIndex]->contiguousBondQ())
			{
				vector<vector<int>> occupancyVector;
				auto iterator=chain[interiorSiteIndex]->occupancy.begin();
				while(iterator !=chain[interiorSiteIndex]->occupancy.end())       //make a vector of the occupancy vectors, max 2
				{
					vector<int> myOccupancyVector=*iterator;
					occupancyVector.push_back(myOccupancyVector);
					++iterator;
				}
				vector<int>contigSiteIndices={occupancyVector[0][2],occupancyVector[1][2]};
				sort(contigSiteIndices.begin(), contigSiteIndices.end());

				bool directionalExpansionCondition1=interiorSiteIndex==contigSiteIndices[0]&&neighborDirection==1;
				bool directionalExpansionCondition2=interiorSiteIndex==contigSiteIndices[1]&&neighborDirection==-1;
				fullExpansionCondition=directionalExpansionCondition1||directionalExpansionCondition2;
			}

			if(newHeadSite!=nullptr&&fullExpansionCondition)   //allowed to do this if you find a bond and the head can move
			{
				interiorSite=chain[interiorSiteIndex];
			}
				

			if(interiorSite!=nullptr)     //proceed if you've chosen an eligible site and direction
			{

				//first, check bonds at current chain and pairs at current site
				vector<int> startBonds={0,0,0}; //start with new site we'll occupy
				float pairsI=0;
				set<site*> uniqueStartSites;  //set of pointers to sites, without repeats from contiguous bonds
				for(int i=0;i<length;++i)
				{
					uniqueStartSites.insert(chain[i]);
				}

				auto uniqueStartIterator=uniqueStartSites.begin();
				while(uniqueStartIterator!=uniqueStartSites.end())
				{
					site* thisSite=*uniqueStartIterator;
					startBonds=vectorAddition(startBonds,getBondInfo(thisSite),1);

					if(nonspecificStatus==true)
					{
						float pairChange=thisSite->getPairNumber(polymerID);

						if(thisSite->occupancy.size()==2)                               //if the moving protein has 2 bonds here, count both
						{
							vector<vector<int>> occupancyVector(thisSite->occupancy.begin(), thisSite->occupancy.end());  //initialize vector of occupancy vectors
							if(occupancyVector[0][0]==polymerID&&occupancyVector[1][0]==polymerID)
							{
								pairChange=pairChange*2;                        //if both the monomers could have moved, then check both their nonspecific interactions
							}
						}

						pairsI+=pairChange;
					}

					++uniqueStartIterator;
				}

				if(uniqueStartSites.find(newHeadSite)==uniqueStartSites.end()) //if we've abandoned a site, add it too
				{
					startBonds=vectorAddition(startBonds,getBondInfo(newHeadSite),1);
				}

				//check for contiguous bonds before moving anything
				vector<bool> contigBondRecord;
				bool thisRecord;
				for(int i=0;i<length;++i)
				{	thisRecord=chain[i]->contiguousBondQ();
					contigBondRecord.push_back(thisRecord);}


				bool specialHeadExpansion=false;            //if you decide to expand from the head, you don't need to move everything else

				if(neighborDirection==-1)                             //pick the correct index and move interior site to avoid normal contiguous move
				{
					vector<int> myOccupancy={polymerID,sequence[interiorSiteIndex-1],(interiorSiteIndex-1)};
					if(interiorSiteIndex!=1)
					{moveUnit(myOccupancy,chain[interiorSiteIndex-1],chain[(interiorSiteIndex-2)],moveRecord);}
					else      //nb special case for if you're trying to expand a contiguous bond at the head
						{moveUnit(myOccupancy,chain[interiorSiteIndex-1],newHeadSite,moveRecord);
						specialHeadExpansion=true;}
				}
				else if(neighborDirection==1)
				{
					vector<int> myOccupancy={polymerID,sequence[interiorSiteIndex+1],(interiorSiteIndex+1)};
					if(interiorSiteIndex!=(length-2))
					{moveUnit(myOccupancy,chain[interiorSiteIndex+1],chain[(interiorSiteIndex+2)],moveRecord);}
					else                           //special case of expanding contiguous bond at head
					{moveUnit(myOccupancy,chain[interiorSiteIndex+1],newHeadSite,moveRecord);
					specialHeadExpansion=true;}
				}

				if(!specialHeadExpansion) //only do all this if you didn't finish by just moving the head
				{
					if(neighborDirection==-1)  //[0] end is moving, so we start at the neighbor of interior site and move towards 0
					{
						for(int i=(interiorSiteIndex-2);i>=2;--i) 
						{
							if(!contigBondRecord[i])
							{
							vector <int> myOccupancy={polymerID,sequence[i],i};
							moveUnit(myOccupancy,chain[i],chain[i-1],moveRecord);
							}
							else
							{
							vector <int> myOccupancy1={polymerID,sequence[i],i};
							vector <int> myOccupancy2={polymerID,sequence[i-1],i-1};
							moveUnit(myOccupancy1,chain[i],chain[i-2],moveRecord);
							moveUnit(myOccupancy2,chain[i-1],chain[i-2],moveRecord);
							--i;
							}
						}
					}
					else               //[length-1] end is moving, so we start at neighbor of interior site and move to [l-1]
					{
						for(int i=(interiorSiteIndex+2);i<length-2;++i)
						{
							if(!contigBondRecord[i])
							{
							vector <int> myOccupancy={polymerID,sequence[i],i};
							moveUnit(myOccupancy,chain[i],chain[i+1],moveRecord);
							}
							else
							{
							vector <int> myOccupancy1={polymerID,sequence[i],i};
							vector <int> myOccupancy2={polymerID,sequence[i+1],i+1};
							moveUnit(myOccupancy1,chain[i],chain[i+2],moveRecord);
							moveUnit(myOccupancy2,chain[i+1],chain[i+2],moveRecord);
							++i;
							}
						 }
					}
					//move head and its neighbor. note the special cases for if the head is contiguous
					int end;
					if(neighborDirection==-1)
						{end=0;}
					else if(neighborDirection==1)
						{end=length-1;}

					if(!contigBondRecord[end])
					{
						if(end==0)        //move next-to-end
						{
							vector <int> myOccupancy={polymerID,sequence[1],1};
							moveUnit(myOccupancy,chain[1],chain[end],moveRecord);
						}
						else
						{
							int thisIndex=length-2;
							vector <int> myOccupancy={polymerID,sequence[thisIndex],thisIndex};
							moveUnit(myOccupancy,chain[thisIndex],chain[end],moveRecord);
						}
						vector<int> leadingOccupancy={polymerID,sequence[end],end}; //the moving end
						moveUnit(leadingOccupancy,chain[end],newHeadSite,moveRecord);
						
					}
					else      //move contiguous bond head and penultimate together
					{
						vector<int> leadingOccupancy={polymerID,sequence[end],end}; //the moving end
						moveUnit(leadingOccupancy,chain[end],newHeadSite,moveRecord);
						if(end==0)
						{
							vector <int> myOccupancy={polymerID,sequence[1],1};
							moveUnit(myOccupancy,chain[1],newHeadSite,moveRecord);
						}
						else
						{
							int thisIndex=length-2;
							vector <int> myOccupancy={polymerID,sequence[thisIndex],thisIndex};
							moveUnit(myOccupancy,chain[thisIndex],newHeadSite,moveRecord);

						}
					}

				}

				//now check bonds and pairs at end
				vector<int> endBonds{0,0,0}; 
				float pairsF=0;
				set<site*> uniqueEndSites;  //set of pointers to sites, without repeats from contiguous bonds
				for(int i=0;i<length;++i)
				{
					uniqueEndSites.insert(chain[i]);
				}
				auto uniqueEndIterator=uniqueEndSites.begin();
				while(uniqueEndIterator!=uniqueEndSites.end())
				{
					site* thisSite=*uniqueEndIterator;
					endBonds=vectorAddition(endBonds,getBondInfo(thisSite),1);

					if(nonspecificStatus==true)
					{
						float pairChange=thisSite->getPairNumber(polymerID);

						if(thisSite->occupancy.size()==2)                               //if the moving protein has 2 bonds here, count both
						{
							vector<vector<int>> occupancyVector(thisSite->occupancy.begin(), thisSite->occupancy.end());  //initialize vector of occupancy vectors
							if(occupancyVector[0][0]==polymerID&&occupancyVector[1][0]==polymerID)
							{
								pairChange=pairChange*2;                        //if both the monomers could have moved, then check both their nonspecific interactions
							}
						}

						pairsF+=pairChange;
					}

					++uniqueEndIterator;
				}

				neighborPairChange=int(pairsF-pairsI);
				bondChanges=vectorAddition(endBonds,startBonds,-1);

				successfulMove=true;



			}

		}
		break;

		case 4:    //corner move
		{
			if(chain.size()>2)
			{
				uniform_int_distribution<int> interiorSiteDist(1, (length-2));
				int randSiteIndex=interiorSiteDist(rng); //randomly choose interior site

				if(chain[randSiteIndex]->contiguousBondQ())        //move whole bond at once, if the protein is longer than 2 monomers
				{
						//find out which subunits are overlapping
						vector<vector<int>> occupancyVector;
						auto iterator=chain[randSiteIndex]->occupancy.begin();
						while(iterator !=chain[randSiteIndex]->occupancy.end())       //make a vector of the occupancy vectors, max 2
						{
						vector<int> myOccupancyVector=*iterator;
						occupancyVector.push_back(myOccupancyVector);
						++iterator;
						}

						vector<int> contiguousUnitIndices={occupancyVector[0][2],occupancyVector[1][2]};
						sort(contiguousUnitIndices.begin(),contiguousUnitIndices.end());

						if(contiguousUnitIndices[0]!=0&&contiguousUnitIndices[1]!=(length-1))       //forbid move if it's the two on the end
						{
							vector<site*> leftNeighbors=chain[(contiguousUnitIndices[0]-1)]->neighbors;
							vector<site*> rightNeighbors=chain[(contiguousUnitIndices[1]+1)]->neighbors;
							
							vector<site*> intersection;
							//find intersection of left and right neighbors
							sort(leftNeighbors.begin(), leftNeighbors.end());
			    			sort(rightNeighbors.begin(), rightNeighbors.end());
			    			set_intersection(leftNeighbors.begin(),leftNeighbors.end(),rightNeighbors.begin(),rightNeighbors.end(),back_inserter(intersection));

			    			//pick a random eligible neighbor, check its validity, and move if it's available
			    			
			    			uniform_int_distribution<int> randMoveDist(0,(intersection.size()-1));
			    			int randomMove=randMoveDist(rng);
							//move both subunits
							if(intersection[randomMove]->occupancy.size()==0&&intersection[randomMove]!=chain[randSiteIndex])       //check to see if destination site is free 
			    			{
			    				int pairsI;
			    				if(nonspecificStatus==true)
			    				{
			    					pairsI=2*chain[randSiteIndex]->getPairNumber();
			    				}

			    				site* movingSite=chain[randSiteIndex];          //pointer to the relevant site. don't reference chain[site] because it might
			    				//change before the second move
			    				moveUnit(occupancyVector[0],movingSite,intersection[randomMove],moveRecord);
			    				moveUnit(occupancyVector[1],movingSite,intersection[randomMove],moveRecord);

			    				if(nonspecificStatus==true)
			    				{
				    				int pairsF=2*(intersection[randomMove]->getPairNumber());
				    				neighborPairChange=pairsF-pairsI;
			    				}
			    				//this version cannot change bonds
			    				bondChanges={0,0,0};

			    				successfulMove=true;

			    			}
		    			}
				}
				else          //otherwise, just find intersection of neighbors and move
				{

					vector<site*> leftNeighbors=chain[randSiteIndex-1]->neighbors;
					vector<site*> rightNeighbors=chain[randSiteIndex+1]->neighbors;
					vector<site*> intersection;
					//find intersection of left and right neighbors
					sort(leftNeighbors.begin(), leftNeighbors.end());
	    			sort(rightNeighbors.begin(), rightNeighbors.end());
	    			set_intersection(leftNeighbors.begin(),leftNeighbors.end(),rightNeighbors.begin(),rightNeighbors.end(),back_inserter(intersection));

	
	    			//pick a random eligible neighbor, check its validity, and move if it's available
	    			uniform_int_distribution<int> randMoveDist(0,(intersection.size()-1));
			    	int randomMove=randMoveDist(rng);
	    			if(intersection[randomMove]->occupancy.size()<2&&chain[randSiteIndex]!=intersection[randomMove])       //check to see if moving site is free 
	    			{

	    				vector<int> startBonds=vectorAddition(getBondInfo(chain[randSiteIndex]),getBondInfo(intersection[randomMove]),1);

	    				site* originalSitePointer=chain[randSiteIndex]; //will need to check this after movement, so see if you broke a bond

	    				
	    				int pairsI;
	    				if(nonspecificStatus==true)
	    				{
	    					pairsI=chain[randSiteIndex]->getPairNumber();
	    				}

	  					vector<int> thisOccupant={polymerID,sequence[randSiteIndex],randSiteIndex};
	    				moveUnit(thisOccupant,chain[randSiteIndex],intersection[randomMove],moveRecord);

	    				if(nonspecificStatus==true)
	    				{
		    				int pairsF=intersection[randomMove]->getPairNumber();
		    				neighborPairChange=pairsF-pairsI;
	    				}


	    				vector<int> endBonds=vectorAddition(getBondInfo(originalSitePointer),getBondInfo(intersection[randomMove]),1);
	    				bondChanges=vectorAddition(endBonds,startBonds,-1);

	    				successfulMove=true;
	    			} 	
		    				

	    		}
			}
		}
		break;
	}
	return make_tuple(successfulMove,bondChanges,neighborPairChange,moveRecord);
}