///element-wise addition of vectors, useful for calculating total bond changes
//use sign =+1 for addition, sign=-1 for subtraction

using namespace std;

vector<int> vectorAddition(vector<int> v1,vector<int> v2, int sign) 
{

	vector<int> v3;
	for(int i=0;i<v1.size();++i)
	{	
		int thisTerm=v1[i]+sign*v2[i];
		v3.push_back(thisTerm);
	}

	return v3;

}