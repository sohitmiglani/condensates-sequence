using namespace std;

vector<vector<vector<site>>> makeSpace(tuple<int,vector<int>,int,vector<vector<int>>> latticeInfo)
{
	//convert lattice data into simple params
	int dim=get<0>(latticeInfo);
	vector<int> size(dim);
	size=get<1>(latticeInfo);
	int neighbors=get<2>(latticeInfo);
	vector<vector<int>> edges;
	edges=get<3>(latticeInfo);

//initialize an mxnxl vector of default sites, in order of z, y, x
vector<site> c(size[2],site());
vector<vector<site>> b(size[1],c);
vector<vector<vector<site>>> points(size[0],b);

for(int i=0;i<size[0];++i)
{
	for(int j=0;j<size[1];++j)
	{
		for(int k=0;k<size[2];++k)
		{
			points[i][j][k].position={i,j,k};
			for(int l=0;l<neighbors;++l)
			{
				auto thisEdge=edges[l];
				//cout<<(i+thisEdge[0]+size[0])%size[0]<<" "<<(j+thisEdge[1]+size[1])%size[1]<<" "<<(k+thisEdge[2]+size[2])%size[2]<<endl;
				points[i][j][k].addNeighbor(&points[(i+thisEdge[0]+size[0])%size[0]][(j+thisEdge[1]+size[1])%size[1]][(k+thisEdge[2]+size[2])%size[2]]);
			}
		}
	}
}

	return points;
}