//============================================================================
// Name        : zft.cpp
// Author      : loywong
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cmath>
#include <iostream>
#include "../LoyMath/Complex.h"
using namespace std;
using namespace LoyMath;

#define TF(x) ((1 - zpow(-1000, (x)))/(1 - zpow(-1, (x))))

int testComplex()
{
    char str[64];
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	Cplx<double> a(1.0, 2.0);
	Cplx<double> b(3, 4);
	Cplx<double> c = a ^ b;
    c = a + b;
    cout << c.ToString(str, 't') << endl;
    c = a - b;
    cout << c.ToString(str, 't') << endl;
    c = a * b;
    cout << c.ToString(str, 't') << endl;
    c = a / b;
    cout << c.ToString(str, 't') << endl;
    c = a.Conj();
    cout << c.ToString(str, 't') << endl;
    a.ConjMe();
    cout << a.ToString(str, 't') << endl;
    c = a.Conj();
    cout << c.ToString(str, 't') << endl;
    c *= a;
    cout << c.ToString(str, 't') << endl;
    c += a;
    cout << c.ToString(str, 't') << endl;
//	int i = sizeof(Cplx<double>);
	cout << c.ToString(str, 't') << endl;
	Cplx<double> resp = TF(0.00598398601);
    double gain = 20 * std::log10(resp.Abs());
	cout << gain << endl;
	return 0;
}
