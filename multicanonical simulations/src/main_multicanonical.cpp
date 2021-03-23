#include "head.hpp"
#include "randomNumberGenerator.hpp"

using namespace std;

 
mt19937 rng;

int main(int argc, char* argv[]) {

	// random_device rd;
	// rng.seed(rd());

	bool writeFullConfig=false;

	string filename_lattice;
	string filename_poly;
	string test_id;
	string tScheduleFile;
	int numProteins;
	long long int steps;              
	int writeInterval;
	double exchangeJ;         //non-specific interaction in units of binding energy epsilon (effectively the t scale)
	double mu;               //chemical potential in units of binding energy epsilon
	string preweightingConfig;
	int seedMod;

	if (argc < 4) { 
        cerr << "Enter test ID, lattice file, polymer file, T schedule, and # of proteins." <<endl;
        cerr<<"Try lattice_fcc.txt, polySpecs.txt."<<endl;
        return 1;
    }




	test_id=argv[1];
	filename_lattice=argv[2];
	filename_poly=argv[3];
	tScheduleFile=argv[4];
	numProteins=atoi(argv[5]);
	steps=atoll(argv[6]);
	writeInterval=atoi(argv[7]);
	exchangeJ=atof(argv[8]);            //units of energy where epsilon=1
	preweightingConfig=argv[9];


	if (argc<10)
    {
    	seedMod=0;
    }
    else
	{seedMod=atoi(argv[10]);}


	//use seedMod to have a different random seed for simultaneous simulations
	srand(time(0)+seedMod);
	rng.seed(rand());

//lattice info includes dim, size, #neighbors, edge vectors
	tuple<int,vector<int>,int,vector<vector<int>>> latticeInfo;
	latticeInfo=readLattice(filename_lattice);
	vector<int> size=get<1>(latticeInfo);
	float volume=1;
	for(int i=0;i<size.size();i++)
	{
		volume=volume*size[i];
	}
	float numNeighbors=get<2>(latticeInfo);


//poly info has length and sequence
	tuple<int,vector<int>> polyInfo;
	polyInfo=readPoly(filename_poly);
	float length=get<0>(polyInfo);

//tSchedule has MC steps at each temperture
	vector<tuple<double,int>> tSchedule;
	tSchedule=readTSchedule(tScheduleFile);

//read preweighting vector
	tuple<double,map<int,double>> preweightingInfo;
	map<int,double> preweighting;
	preweightingInfo=readPreweightingFile(preweightingConfig);

	mu=get<0>(preweightingInfo);
	preweighting=get<1>(preweightingInfo);



//write parameter file
	//int numProteins=10;
	paramRecord(test_id,latticeInfo,polyInfo,tSchedule,numProteins,steps,writeInterval,exchangeJ,mu);

//create list of sites with neighbors
	vector<vector<vector<site>>> realSpaceSystem;
	realSpaceSystem=makeSpace(latticeInfo);

//initialize vector of 10 randomly-placed,self-avoiding polymers, (size,polyInfo,space,ID)
//then write it to file
string polymerFile="polyOutput_"+test_id+".dat";
string writeRecordFile="writeRecord_"+test_id+".dat";



unordered_map<int, polymer> polymerMap;
vector<int> keyVector;
// //standard random initialization**********************************************************************
for(int i=0;i<numProteins;++i)
{
	polymerMap.insert(make_pair(i,polymer()));
	polymerMap[i].initialize(latticeInfo,polyInfo,realSpaceSystem,i,"open");

	//keep track of existing keys for simple access to random elements via ID
	keyVector.push_back(i);
}
//***************************************************************************************************

//ID for next protein inserted
int polymerIDindex=numProteins;


//initialize temperature schedule
double beta;
tuple<double,int> thisTStep=tSchedule[0];
beta=get<0>(thisTStep);
int tempInterval=get<1>(thisTStep);
int tempCounter=0;
int currentTempInterval=0;
int numTempIntervals=tSchedule.size();

//save initial configuration
ofstream writeOutput(writeRecordFile,fstream::app);
// ofstream polyOutput(polymerFile,fstream::app);
int thisNp=polymerMap.size();
double energy=getEnergy(realSpaceSystem,size,exchangeJ);
writeOutput<<0<<" "<<beta<<" "<<thisNp<<" "<<energy<<endl;
// for(int j=0;j<keyVector.size();++j)
// {
// int thisKey=keyVector[j];
// polymerMap[thisKey].writePolyChain(polymerFile);
// }


for(long long int i=0;i<steps;++i)
{
	if(tempCounter>=tempInterval&&currentTempInterval<(numTempIntervals-1)) //if you've done as many steps as specified MC steps for this T, move to next T
	{
		currentTempInterval++;
		thisTStep=tSchedule[currentTempInterval];
		beta=get<0>(thisTStep);
		tempInterval=get<1>(thisTStep);
		tempCounter=0;
	}

	monteCarlo(polymerMap,realSpaceSystem,beta,polymerFile,writeRecordFile,i,writeInterval,tempCounter,exchangeJ,mu,volume,polyInfo,polymerIDindex,keyVector,numNeighbors,length,energy,preweighting,writeFullConfig);
	tempCounter++;

	if(i%writeInterval==0&&i!=0)
	{
		//write bond strength and number of proteins
		ofstream writeOutput(writeRecordFile,fstream::app);
		// ofstream polyOutput(polymerFile,fstream::app);
		int thisNp=polymerMap.size();
		writeOutput<<i<<" "<<beta<<" "<<thisNp<<" "<<energy<<endl;

		//write full polymer configuration
		if(writeFullConfig==true)
		{
			for(int j=0;j<keyVector.size();++j)
			{
			int thisKey=keyVector[j];
			polymerMap[thisKey].writePolyChain(polymerFile);
			}
		}

	}
}



	return 0;
}
