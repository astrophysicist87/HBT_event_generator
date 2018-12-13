#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

ostream & operator << (ostream &out, const vector<int> &v)
{
	int k = 0;
	out << "{" << v[k];

	for (k = 1; k < v.size(); ++k)
		out << ", " << v[k];
	out << "}";

    return out;
}

int main(int argc, char *argv[])
{
	const int N = 116;

	vector<int> long_vec;

	for (int k = 0; k < N; ++k)
		long_vec.push_back(k);

	cout << "long_vec = " << long_vec << endl;

	// do partition
	const int n = 17;	// main partition length
	const int l = (int)ceil(double(N)/double(n));	// length of list of partitions
	
	vector<int> lowerLimits, upperLimits;
	for (int il = 0; il < l; ++il)
	{
		lowerLimits.push_back(il*n);
		upperLimits.push_back((il+1)*n-1);
	}
	// make sure last element is correct
	upperLimits[l-1] = long_vec[N-1];

	cout << "lowerLimits = " << lowerLimits << endl;

	cout << "upperLimits = " << upperLimits << endl
			<< "------------------------------------------------------"
			<< endl;

	// now, print corresponding sub-vectors
	for (int il = 0; il < l; ++il)
	{
		vector<int> tmp;
		for (int ik = lowerLimits[il]; ik <= upperLimits[il]; ++ik)
			tmp.push_back(long_vec[ik]);
		cout << "Sub-vector #" << il << " = " << tmp << endl;
	}

	return (0);
}

//End of file
