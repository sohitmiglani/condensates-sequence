using namespace std;

vector<tuple<double,int>> readTSchedule(string filename)
{
	ifstream scheduleSource(filename);

	 if (!scheduleSource) {
	        cout << "Unable to open t schedule";
	        exit(1); // terminate with error
	    };
	string dummy;
	getline(scheduleSource,dummy);

	vector<tuple<double,int>> schedule;
	double intermediateE;
	int intermediateStep;

	while(scheduleSource>>intermediateE>>intermediateStep)
	{
		tuple<double,int> intermediatePair{intermediateE,intermediateStep};
		schedule.push_back(intermediatePair);
	}


	scheduleSource.close();

	return schedule;
}


