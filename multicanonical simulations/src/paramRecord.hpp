//this function creates a parameter output file with a test ID to match the model output

#include "head.hpp"

using namespace std;

void paramRecord(string test_id, auto latticeInfo, auto polyInfo, auto tSchedule,int numProteins,int steps, int writeInterval, double exchangeJ, double mu)
{

	cout.precision(15);
	//file handling
		string filename="params_"+test_id+".dat";
		ofstream params(filename,fstream::app);
		 if (!params)
		 {
	        cout << "Unable to open file";
	        exit(1); // terminate with error
	    };

	    string latticeFormat="lattice: dim,size,neighbors,edges";
	 	string polyFormat="polymer: length, sequence, number of proteins";

	    params<<latticeFormat<<endl;

//reload all info
	    int dim=get<0>(latticeInfo);
	    vector<int> size=get<1>(latticeInfo);
	    int neighbors=get<2>(latticeInfo);
	    vector<vector<int>> edges=get<3>(latticeInfo);

	    int length=get<0>(polyInfo);
	    vector<int> sequence=get<1>(polyInfo);
//write lattice info
	    params<<dim<<endl;
	    for (int i = 0; i < dim; ++i)
	    	params<<size[i]<<" ";
	    params<<endl;
	    params<<neighbors<<endl;

	    for (int i = 0; i < neighbors; ++i)
	    {
	    	for (int j = 0; j < dim; ++j)
	    	{
	    		params<<edges[i][j]<<" ";
	    	}
	    	params<<endl;
	    }
//write polymer info
	    params<<polyFormat<<endl;
	    params<<length<<endl;
	    for(int i=0;i<length;++i)
	    	params<<sequence[i]<<" ";
	    params<<endl<<numProteins<<endl;

//write simulation details
	    params<<"simulation params: steps, write interval, J, mu:"<<endl;
	    params<<steps<<endl;
	    params<<writeInterval<<endl;
	    params<<exchangeJ<<endl;
	    params<<fixed<<setprecision(15)<<mu<<endl;

//write annealing schedule
	    params<<"t schedule:"<<endl;
	    double thisT;
	    int thisInterval;
	    tuple<double,int> currentTStep;
	    for(int i=0;i<tSchedule.size();++i)
	    {	
	    	currentTStep=tSchedule[i];
	    	thisT=get<0>(currentTStep);
	    	thisInterval=get<1>(currentTStep);
	    	params<<thisT<<" "<<thisInterval<<endl;
	    }
	    params.close();

}