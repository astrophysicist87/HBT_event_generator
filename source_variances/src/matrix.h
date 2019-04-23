#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

inline double det2D(double [] a)
{
	// a = {a11, a12, a21, a22}
	//	 = {0,   1,   2,   3  }
	return (
			a[0]*a[3] - a[1]*a[2]
			);
}

inline double det3D(double [] a)
{
	// a = {a11, a12, a13, a21, a22, a23, a31, a32, a33}
	//	 = {0,   1,   2,   3,   4,   5,   6,   7,   8  }
	return (
			a[0]*a[4]*a[8] + a[1]*a[5]*a[6] + a[2]*a[3]*a[7]
			-a[0]*a[5]*a[7] - a[2]*a[4]*a[6] - a[1]*a[3]*a[8]
			);
}

#endif
