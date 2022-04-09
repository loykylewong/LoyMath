//============================================================================
// Name        : mattest.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;
#include "../LoyMath/Matrix.h"
#include <cmath>
#include <chrono>

using namespace LoyMath;

int testMatrix() {
//	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	char s[256];
//	float x = 0.f, y, z;

	Mat<2,2> a(4, 1.f, 2.f, 3.f, 4.f);
	Mat<2,1> b(2, 2.f, 3.f);
	Mat<2,1> c;
	Mat<2,4> d;
	cout << a(0,0) << '\t' << a(0,1) << endl;
	cout << a(1,0) << '\t' << a(1,1) << endl;
	cout << a.PrintStr(s) << endl;
	cout << s << endl;
	cout << a.ToStr(s) << endl;
//	a(0,0) = 5.5f;
	cout << a.ToStr(s) << endl;

	cout << a.ToStr(s) << endl;
	c = a * b;
	cout << c.ToStr(s) << endl;
	d = a.PadRight(Mat<2,2>(1.f));
	cout << d.ToStr(s) << endl;
	Mat<2,2> e;
	d.SubMat(0,2,e);
	cout << e.ToStr(s) << endl;

	Mat<3,3> f(9,
			0.f,1.f,2.f,
			3.f,4.f,5.f,
			6.f,7.f,8.f);
	f.SwapRow(0,2);
	cout << f.ToStr(s) << endl;

	Mat<4,2> g = d.Trans();
	cout << g.ToStr(s) << endl;
	f.TrMe();
	cout << f.ToStr(s) << endl;
	f.TrMe();
	f.SwapRow(0,2);
	f(0,0) = 1.f;
	Mat<3,3> h = f.Inv();
	cout << h.ToStr(s, "%-8f") << endl;
	cout << (f*h).ToStr(s) << endl;
	cout << Mat<2,2>(4,1.f,2.f,2.f,2.f).Inv().ToStr(s) << endl;

	Mat<6,6,float> i(36,
			 1.f,  1.f,   1.f,    1.f,    1.f,    1.f,
			 1.f,  2.f,   3.f,    4.f,    5.f,    6.f,
			 1.f,  4.f,   9.f,   16.f,   25.f,   36.f,
			 1.f,  8.f,  27.f,   64.f,  125.f,  216.f,
			 1.f, 16.f,  81.f,  256.f,  625.f, 1296.f,
			 1.f, 32.f, 243.f, 1024.f, 3125.f, 7776.f);
	cout << i.Inv().ToStr() << endl;

	Mat<6,1> x(6, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f);
	Mat<6,1> y(6, 1.f, 3.f, 8.f, 17.f, 24.f, 33.f);
	Mat<6,3> X;
	for(int r = 0; r < 6; r++)
		for(int c = 0; c < 3; c++)
			X(r,c) = powf(x(r,0), c);
	Mat<3,1> coef = (X.Trans() * X).Inv() * X.Trans() * y;
	cout << coef.ToStr(s) << endl;
	cout << X.Trans().ToStr(s) << endl;
	cout << (X.Trans() * X).ToStr(s) << endl;
	cout << (X.Trans() * X).Inv().ToStr(s) << endl;
    auto tp0 = chrono::steady_clock::now();
    auto XIP = (X.Trans() * X).Inv() * X.Trans();
    auto tp1 = chrono::steady_clock::now();
	cout << ((X.Trans() * X).Inv() * X.Trans()).ToStr(s) << endl;
	cout << ((X.Trans() * X).Inv() * X.Trans() * y).ToStr(s) << endl;
    auto tp2 = chrono::steady_clock::now();
    auto XIP2 = X.PInv(0.001f);
    auto tp3 = chrono::steady_clock::now();
    double d0 = chrono::duration<double, ::micro>(tp1 - tp0).count();
    double d1 = chrono::duration<double, ::micro>(tp3 - tp2).count();
    cout << d0 << endl << d1 << endl;
    
    cout << X.PInv().ToStr(s) << endl;
    
    Mat<3,3> m33(9, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);
    Mat<3,3> m33i, m33pi;
    m33i = m33.Inv();
    m33pi = m33.PInv();
    
    Mat<3,2> m32(6, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0);
    Mat<2,3> m32pi;
    m32pi = m32.PInv();
    
    Mat<2,3> m23(6, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0);
    Mat<3,2> m23pi;
    m23pi = m23.PInv();

    const int m = 10, n = 5;
    Mat<m,m,double> mm, mmi, mmpi;
    Mat<m,n,double> mn;
    Mat<n,m,double> nm;
	for(int r = 0; r < m; r++)
		for(int c = 0; c < m; c++)
		{
            mm(r,c) = sin((double)(r + c));
		}
    for(int r = 0; r < m; r++)
        for(int c = 0; c < n; c++)
        {
            mm(r,c) = sin((double)(r + c));
        }
    mmi = mm.Inv();
    mmpi = mm.PInv();
    nm = mn.PInv();
    
    Mat<2,2,double> md22(4, 1.0, 2.0, 3.0, 4.0);
    cout << md22.ToStr() << endl << md22.Det() << endl;
    Mat<3,3,double> md33(9, 1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0,
                            7.0, 8.0, 0.0);
    cout << md33.ToStr() << endl << md33.Det() << endl;
    
    Mat<4,4,float> md44(16, 1.0f, 2.0f,  3.0f,  4.0f,
                            1.0f, 3.0f,  5.0f,  7.0f,
                            1.0f, 4.0f,  9.0f, 16.0f,
                            1.0f, 8.0f, 27.0f, 64.0f);
    cout << md44.ToStr() << endl << md44.Det() << endl;
    cout << md44.CycAdjoint(0, 0).ToStr() << md44.CycAdjoint(0, 0).Det() << endl;
    cout << md44.CycAdjoint(1, 1).ToStr() << md44.CycAdjoint(1, 1).Det() << endl;
    cout << md44.CycAdjoint(2, 2).ToStr() << md44.CycAdjoint(2, 2).Det() << endl;
    cout << md44.CycAdjoint(3, 3).ToStr() << md44.CycAdjoint(3, 3).Det() << endl;
    cout << md44.CycAdjoint(0, 1).ToStr() << md44.CycAdjoint(0, 1).Det() << endl;
    cout << md44.CycAdjoint(0, 2).ToStr() << md44.CycAdjoint(0, 2).Det() << endl;
    cout << md44.CycAdjoint(0, 3).ToStr() << md44.CycAdjoint(0, 3).Det() << endl;
    cout << md44.CycAdjoint(1, 0).ToStr() << md44.CycAdjoint(1, 0).Det() << endl;
    cout << md44.CycAdjoint(2, 0).ToStr() << md44.CycAdjoint(2, 0).Det() << endl;
    cout << md44.CycAdjoint(3, 0).ToStr() << md44.CycAdjoint(3, 0).Det() << endl;
    cout << md44.CycAdjoint(1, 2).ToStr() << md44.CycAdjoint(1, 2).Det() << endl;

    Vec<3> v0(3, 3.f, 4.f, 0.f);
    Vec<3> v1(3, 4.f, 3.f, 0.f);
    Vec<3> v2 = Cross(v0, v1);
    float f3 = Dot(v0, v1);
    // 4d-vector cross
    Vec<4> v4(4, 1.f, 2.f, 3.f, 4.f);
    Vec<4> v5(4, 1.f, 3.f, 5.f, 7.f);
    Vec<4> v6(4, 1.f, 4.f, 9.f, 16.f);
    Mat<4,4> mc;
    mc.CopyFrom(v4, 0, 0, 0, 1);
    mc.CopyFrom(v5, 0, 0, 0, 2);
    mc.CopyFrom(v6, 0, 0, 0, 3);
    cout << mc.Det() << endl;
    int sign = 1;
    cout << mc.Adjoint(0, 0, sign).ToStr() << mc.Adjoint(0, 0, sign).Det() * sign << endl;
    cout << mc.Adjoint(1, 0, sign).ToStr() << mc.Adjoint(1, 0, sign).Det() * sign << endl;
    cout << mc.Adjoint(2, 0, sign).ToStr() << mc.Adjoint(2, 0, sign).Det() * sign << endl;
    cout << mc.Adjoint(3, 0, sign).ToStr() << mc.Adjoint(3, 0, sign).Det() * sign << endl;

    return 0;
}
