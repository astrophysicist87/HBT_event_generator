#include <iostream>
#include <string>
#include <chrono>
#include <random>

using namespace std;

int main()
{
	const int nrolls=1;  // number of experiments
	const int nstars=100;    // maximum number of stars to distribute

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();

	default_random_engine generator (seed);
	normal_distribution<double> distribution(5.0,2.0);

	int p[10]={};

	for (int i = 0; i < nrolls; ++i)
	{
		double number = distribution(generator);
		if ( ( number >= 0.0 ) and ( number < 10.0 ) )
			++p[int(number)];
	}

	cout << "normal_distribution (5.0,2.0):" << endl;

	int norm = 0;
	for (int i = 0; i < 10; ++i)
	{
		cout << p[i] << "   " << i << "-" << (i+1) << ": ";
		cout << string(p[i]*nstars/nrolls,'*') << endl;
		norm += p[i];
	}

	cout << norm << endl;

	return 0;
}
