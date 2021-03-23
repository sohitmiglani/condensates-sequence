using namespace std;

tuple<int,vector<int>> readPoly(string filename)
{
	ifstream polySource(filename);

	 if (!polySource) {
	        cout << "Unable to open poly file";
	        exit(1); // terminate with error
	    };
	string dummy;
	getline(polySource,dummy);
	int length;

	polySource >> length;

	vector <int> sequence(length);

	for(int i=0;i<length;++i)
	{
		polySource>>sequence[i];
	};
	polySource.close();

	return make_tuple(length,sequence);
}


