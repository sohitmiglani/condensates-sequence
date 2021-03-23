using namespace std;

// void clearSpace(vector<vector<vector<site>>> &space)
// {
// 	for(int i=0;i<space.size();++i)                      //iterate over all space
// 	{
// 		for(int j=0;j<space[0].size();++j)
// 		{
// 			for(int k=0;k<space[0][0].size();++k)
// 			{
// 				space[i][j][k].occupancy.clear();
// 			}
// 		}
// 	}
// }

void clearSpace(vector<vector<vector<site>>> &space,vector<polymer> existingPolys)
{
	int polyLength=existingPolys[0].length;
	for(int i=0;i<existingPolys.size();++i)
	{
		for(int j=0;j<polyLength;++j)
		{
			existingPolys[i].chain[j]->occupancy.clear();
		}
	}
	
}