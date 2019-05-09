#ifndef ESTIMATE_ERROR_H
#define ESTIMATE_ERROR_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

using namespace std;

double estimate_ratio_error(
			double num, double den,
			double num2, double den2, double numden,
			double nSamples,
			bool verbose, std::ostream &out
			)
{
	double result = 0.0;

	const double prefactor = nSamples / ( nSamples - 1.0 );

	// numerator, denominator, and covariances
	double EA = num / nSamples;
	double EB = den / nSamples;
	double EA2 = num2 / nSamples;
	double EB2 = den2 / nSamples;
	double EAB = numden / nSamples;

	// set variances and covariance
	double sigA2 = prefactor * (EA2 - EA*EA);
	double sigB2 = prefactor * (EB2 - EB*EB);
	double sigAB = prefactor * (EAB - EA*EB);

	// want standard error, not variance itself
	sigA2 /= nSamples;
	sigB2 /= nSamples;
	sigAB /= nSamples;

	// set relative widths
	double cA = sigA2 / ( EA*EA+1.e-100 );
	double cB = sigB2 / ( EB*EB+1.e-100 );
	double cAB = sigAB / ( EA*EB+1.e-100 );

	double disc = cA + cB - 2.0*cAB;
	if (verbose and disc < -1.e-6)
		out << "Warning in estimate_error(): "
			<< "disc < 0!" << endl;

	result =
		( disc < 0.0 or den < 1.e-25 )
		? 1.e+6
		: abs( num / den ) * sqrt( cA + cB - 2.0*cAB );
	/*try
	{
		result = abs( num / den ) * sqrt( cA + cB - 2.0*cAB );
	}
	catch (...)
	{
		out << "Exception in estimate_error(): "
			<< num << "   " << den << "   "
			<< cA << "   " << cB << "   " << - 2.0*cAB
			<< "   " << cA + cB - 2.0*cAB << endl;
		result = 1.0e+6;
	}*/

	if (verbose)
	{
		out << setprecision(8);
		out << "\t\t EA = " << EA << endl;
		out << "\t\t EB = " << EB << endl;
		out << "\t\t EA2 = " << EA2 << endl;
		out << "\t\t EB2 = " << EB2 << endl;
		out << "\t\t EAB = " << EAB << endl;

		out << "\t\t sigA2(std.err) = " << sigA2 << endl;
		out << "\t\t sigB2(std.err) = " << sigB2 << endl;
		out << "\t\t sigAB(std.err) = " << sigAB << endl;

		out << "\t\t cA = " << cA << endl;
		out << "\t\t cB = " << cB << endl;
		out << "\t\t cAB = " << cAB << endl;

		out << "\t\t disc = " << disc << endl;
	}

	return (result);
}

#endif
