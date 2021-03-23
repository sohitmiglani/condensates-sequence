//this function takes a specified polymer and the space it lives in, and returns the set of polymers in a connected cluster
//we use a Breadth-first traversal of the graph
//traverse by ID, not vector index

using namespace std;

set<int> getCluster(int thisID,unordered_map<int,polymer> polyMap, auto space)
{
	set<int> proteinCluster;
	deque<int> searchQueue;            //FIFO search queue
	searchQueue.push_back(thisID);

	while(searchQueue.size()>0)
	{	
		int searchingIndex=searchQueue[0];
		int thisLength=polyMap[searchingIndex].length;
		int thisID=polyMap[searchingIndex].polymerID;
		proteinCluster.insert(thisID);  //add the protein we're searching to the cluster
		set<int> idToSearch;    //the set of IDs of proteins connected to this search protein

		for(int i=0;i<thisLength;++i)
		{
			site *currentSite=polyMap[searchingIndex].chain[i];   //pointer to the chain site

			auto iterator=currentSite->occupancy.begin();            //iterate over the occupancy vectors of this site
			while(iterator !=currentSite->occupancy.end())
			{
				vector<int> myOccupancy=*iterator;
				if(myOccupancy[0]!=thisID)              //don't bother adding the current search protein, which is already in cluster
				{
					idToSearch.insert(myOccupancy[0]);
				}
				++iterator;
			}
		}

		//now search cluster to see if we've found new proteins
		auto searchIDIterator=idToSearch.begin();
		while(searchIDIterator!=idToSearch.end())
		{
			int currentSearchID=*searchIDIterator;
			if(proteinCluster.find(currentSearchID) == proteinCluster.end())        //add polyID to search queue if not in cluster set
			{
				searchQueue.push_back(currentSearchID);
			}
			++searchIDIterator;
		}

		searchQueue.pop_front();        //remove the current protein from the search queue
	}

	return proteinCluster;
}