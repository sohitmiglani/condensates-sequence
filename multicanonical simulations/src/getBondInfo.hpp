//how many triples, +-, or ++/-- do you have at the given site?
//call before doing the move!

using namespace std;

vector<int> getBondInfo(site *thisSite) 
{

	vector<int> siteMotifs;


	auto occupancyIterator=thisSite->occupancy.begin();
	while(occupancyIterator!=thisSite->occupancy.end())
	{
		vector<int> thisOccupancy=*occupancyIterator;
		siteMotifs.push_back(thisOccupancy[1]);
		++occupancyIterator;
	}


	//now check combos
	//1: triple occupancy
	int triples=0;
	if(siteMotifs.size()>2)
	{
		triples=1;
	}


	//2: single bonds
	int sameBond=0;
	int diffBond=0;

	if(siteMotifs.size()==2)
	{
		if(siteMotifs[0]==siteMotifs[1])
		{
			sameBond=1;
		}
		else
		{
			diffBond=1;
		}
	}


	vector<int> theseBonds{diffBond,sameBond,triples};
	return theseBonds;		

}