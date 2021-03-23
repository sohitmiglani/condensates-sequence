//this function accepts a vector of polymer objects and the lattice, by reference
//it  copies them, proposes a move, and evaluates the energy before and after
//if the new energy is lower or equal, it accepts the move
//if the energy is higher, it accepts the move with probability e^-E/kt
#include "randomNumberGenerator.hpp"

using namespace std;




void monteCarlo(unordered_map<int,polymer> &polyMap,auto &realSpace,double beta,string polymerFile,string writeRecordFile,long long int &moveCounter,int writeInterval, int &tempCounter, double exchangeJ,double mu, float volume, auto polyInfo, int &polymerIDindex, vector<int> &keyVector, float numNeighbors, float length, double &energy,map<int,double> preweighting,bool writeFullConfig)
{

	vector<double> bondEnergies={beta,100,100};  
	//energy changes for making a compatible bond (epsilon=1, so beta*epsilon=beta), making an incompatible bond, and triple occupancy
	double betaJ=beta*exchangeJ; //lower the energy by exchangeJ kT for each pair of neighbors
	double betaMu=beta*mu;             //note that mu can be positive or negative

	int deltaNumProteins=0; //change if there's a successful insertion/deletion
	double prefactorWeighting=0;  //this will be used to accommodate the weighting factors from V/N, and Rosenbluth configurational bias 

	int neighborPairChanges;
	vector<int> bondChanges;
	vector<tuple<vector<int>,site*, site*>> moveRecord;
	int randomMoveType;
	bool tryInsertion;

	vector<site*> deletedCoords;  //keep track of where we deleted 
	int deletionID;                     //record ID of deleted polymer if we need to restore it


	bool successfulMove =false;
	while(successfulMove==false)
	{	
		deltaNumProteins=0;
		prefactorWeighting=0;

		int randKeyVectorIndex;//random polymer
		int randomKey; //the ID of the random polymer

		//number 0 to 4 for single polymer move, 5 for cluster translation, 6 for exchange with reservoir
		if(polyMap.size()>0)
		{
			discrete_distribution<int> moveDist {108,108,9,108,108,1,324};

			//nvt
			// discrete_distribution<int> moveDist {108,108,9,108,108,1,0};

			//cubic nvt
			// discrete_distribution<int> moveDist {108,108,18,108,108,1,0};



			randomMoveType=moveDist(rng);  
			uniform_int_distribution<int> polyDist(0, (polyMap.size()-1));
			randKeyVectorIndex=polyDist(rng);
			randomKey=keyVector[randKeyVectorIndex];

		}
		else
		{
			randomMoveType=6;
		}





		if(randomMoveType<=4)       //individual polymers
		{
			//nb that we pass exchangeJ so we can adaptively choose whether or not to both with neighbor pairs
			tuple<bool,vector<int>,int,vector<tuple<vector<int>,site*, site*>>> moveData;
			moveData=polyMap[randomKey].movePolymer(randomMoveType,exchangeJ);
			successfulMove=get<0>(moveData);
			bondChanges=get<1>(moveData);
			neighborPairChanges=get<2>(moveData);
			moveRecord=get<3>(moveData);


		}
		else if(randomMoveType==5)   //move cluster to which randomPolymer belongs
		{
			tuple<bool,vector<int>,int,vector<tuple<vector<int>,site*, site*>>> moveData;
			moveData=clusterMove(randomMoveType,randomKey,polyMap,realSpace,exchangeJ);
			successfulMove=get<0>(moveData);
			bondChanges=get<1>(moveData);
			neighborPairChanges=get<2>(moveData);
			moveRecord=get<3>(moveData);

		}

		else  //exchange with reservoir
		{
			//randomly choose between insertion and deletion
			bernoulli_distribution particleExchangeDist(0.5);
			tryInsertion=particleExchangeDist(rng);
			if(tryInsertion)        //insertion move
			{	
				deltaNumProteins=1;

				//first try to generate a valid insertion configuration with associated Rosenbluth weight
				tuple<bool,double> moveData;
				polyMap.insert(make_pair(polymerIDindex,polymer()));

				moveData=polyMap[polymerIDindex].growBiased(polyInfo,realSpace,polymerIDindex);



				successfulMove=get<0>(moveData);
				//prod_i=2^L k_i where k_i is number of valid sites 
				double rosenbluth=get<1>(moveData);


				//if the move succeeds,  calculate energy change
				if(successfulMove==true)
				{

					//get energy change
					tuple<vector<int>,int> insertionData;
					insertionData=polyMap[polymerIDindex].insertionEnergy(exchangeJ);
					bondChanges=get<0>(insertionData);
					neighborPairChanges=get<1>(insertionData);


					//add key to keyVector  
					keyVector.push_back(polymerIDindex);
					//increment ID, so next protein will have a new ID
					polymerIDindex++;



				}
				else //remove everything
				{	
					//delete occupancies
					polyMap[polymerIDindex].deletePoly();
					//erease from polyMap
					polyMap.erase(polymerIDindex);


				}

				//see deletion weight


				
				//already inserted it, so size of polyVector is N+1
				//numNeighbors+1 because we allow staying on same site

				prefactorWeighting=log(volume/double(polyMap.size()))+log(rosenbluth)-(length-1)*log(numNeighbors+1);
			}
			else   //deletion move
			{

				if(polyMap.size()>0) //cannot delete without any polymers...
				{
					deltaNumProteins=-1;

					//this always can work
					successfulMove=true;
					
					//already chosen a random polymer above. save its coords and ID, in case move is rejected
					for(int k=0;k<polyMap[randomKey].chain.size();++k)  //deep copy the pointers
					{	
						deletedCoords.push_back(polyMap[randomKey].chain[k]);
					}
					deletionID=randomKey;


					//get Delta E from removing it
					tuple<vector<int>,int> deletionEnergy;
					deletionEnergy=polyMap[randomKey].insertionEnergy(exchangeJ);
					bondChanges=get<0>(deletionEnergy);
					neighborPairChanges=get<1>(deletionEnergy);
					//insertion energy gives you the number of bonds/pairs, so we need to multiply by -1
					neighborPairChanges=(-1)*neighborPairChanges;
					for(int k=0;k<bondChanges.size();++k)
					{
						bondChanges[k]=(-1)*bondChanges[k];
					}
					//get weight and delete from sites. We have overloaded deletePoly, so passing a bool returns the weight
					double rosenbluth=polyMap[randomKey].deletePoly(true);


					//delete from map and keyVector
					polyMap.erase(randomKey);
					keyVector.erase(keyVector.begin()+randKeyVectorIndex);



					//already deleted it, so we need to add to get N+1
					//numNeighbors+1 because we allow staying on same site
					prefactorWeighting=log(double(polyMap.size()+1)/volume)-log(rosenbluth)+(length-1)*log(numNeighbors+1);

				}
			}

		}

		if(successfulMove==false)                             //this is necessary to not undersample failed moves
		{	

			if(moveCounter%writeInterval==0&&moveCounter!=0)
			{

				//write bond strength and number of proteins
				ofstream writeOutput(writeRecordFile,fstream::app);
				// ofstream polyOutput(polymerFile,fstream::app);
				int thisNp=polyMap.size();
				writeOutput<<moveCounter<<" "<<beta<<" "<<thisNp<<" "<<energy<<endl;

				if(writeFullConfig==true)
				{
					for(int j=0;j<keyVector.size();++j)
					{
						int thisKey=keyVector[j];
						polyMap[thisKey].writePolyChain(polymerFile);
					}
				}




			}
			moveCounter++;
			tempCounter++;
		}
	}



	//include Delta U from bonds and nonspecific interactions, plus beta*mu*Delta N
	double deltaEBond=double(-beta*bondChanges[0]+100*bondChanges[1]+100*bondChanges[2]);
	double deltaEPairs=double(neighborPairChanges)*betaJ*(-1);

	double particleWeight=double(deltaNumProteins)*betaMu;

	//take preweighting into account
	double preweightingTerm=1;
	if(deltaNumProteins!=0)
	{
		int trialNp=polyMap.size();
		int oldNp=trialNp-deltaNumProteins;

		preweightingTerm=preweighting[oldNp]/preweighting[trialNp];
		

	}



	double argument=particleWeight-deltaEBond-deltaEPairs+prefactorWeighting+log(preweightingTerm);


	if(argument>=0) //if it lowers or maintains the energy, accept and update E
	{
		double deltaETotal=(-1)*(float(bondChanges[0])+exchangeJ*double(neighborPairChanges));
		energy+=deltaETotal;

	}

	if(argument<0)          //if it raises energy, accept with probability e^-(deltaE)
	{
		double moveProb=exp(argument);

		// random_device device;
  //   	mt19937 generator(device());
    	bernoulli_distribution d(moveProb);

    	bool accepted=d(rng);


    	if(accepted) //update energy
    	{
			double deltaETotal=(-1)*(float(bondChanges[0])+exchangeJ*double(neighborPairChanges));
			energy+=deltaETotal;


    	}


    	if(!accepted) //if you reject the move, you have to change it back
    	{

    		if(randomMoveType==6) //exchange particles with reservoir
    		{
    			if(tryInsertion) //reverse insertion by deleting polymer
    			{
    				//decrease ID so we're back at the one we just inserted
    				polymerIDindex--;
    				polyMap[polymerIDindex].deletePoly(); //removes from sites
    				polyMap.erase(polymerIDindex);

    				//remove from keyVector    				
    				keyVector.pop_back();



    			}
    			else //reverse deletion by re-inserting polymer
    			{

					polyMap.insert(make_pair(deletionID,polymer()));
					polyMap[deletionID].setConfig(polyInfo,deletedCoords,deletionID);
					//add key to keyVector  
					keyVector.push_back(deletionID);    


				}
    		}
    		else
    		{
    			for(int i=0;i<moveRecord.size();++i)
    			{
    			auto thisRecord=moveRecord[i];
    			vector<int> currentOccupancy=get<0>(thisRecord);
    			site* original =get<1>(thisRecord);
    			site* updated =get<2>(thisRecord);
    			polyMap[currentOccupancy[0]].moveUnit(currentOccupancy,updated,original);
    			}
    		}
    		
     		
    	}
	}

}
