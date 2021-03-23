using namespace std;

tuple<double,map<int,double>> readPreweightingFile(string filename)
{
	double mu;

	ifstream preweightingSource(filename);

	 if (!preweightingSource) {
	        cout << "Unable to open preweighting schedule";
	        exit(1); // terminate with error
	    };

	string dummy;
	string intermediate;

	//skip these things which are just for the record
	getline(preweightingSource,dummy);
	getline(preweightingSource,dummy);
	getline(preweightingSource,dummy);
	getline(preweightingSource,dummy);	
	getline(preweightingSource,dummy);
	getline(preweightingSource,dummy);

	//mu
	getline(preweightingSource,dummy);
	getline(preweightingSource,intermediate);
	mu=stod(intermediate);


	getline(preweightingSource,dummy);



	map<int,double> preweighting;
	int Np;
	double weight;

	while(preweightingSource>>Np>>weight)
	{

		preweighting[Np]=weight;

	}



	preweightingSource.close();

	return make_tuple(mu,preweighting);
}


