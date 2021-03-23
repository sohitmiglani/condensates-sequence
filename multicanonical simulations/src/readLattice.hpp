using namespace std;

tuple<int,vector<int>,int,vector<vector<int>>> readLattice(string filename)
{
	ifstream latticeSource(filename);

	 if (!latticeSource) {
	        cout << "Unable to open lattice file";
	        exit(1); // terminate with error
	    };
	string dummy;
	getline(latticeSource,dummy);
	int dim;
	int neighbors;

	latticeSource >> dim;

	vector<int> size;
	int intermediate;
	for (int i=0;i<dim;++i)		
	{
		latticeSource>>intermediate;
		size.push_back(intermediate);
	}


	latticeSource >> neighbors;
	vector <vector<int>> edges;	
	vector <int> eachEdge(3);

	while(latticeSource>>eachEdge[0]>>eachEdge[1]>>eachEdge[2])
	{
		edges.push_back(eachEdge);
	};
	latticeSource.close();

	return make_tuple(dim,size,neighbors,edges);
}


