//
//  testLinApprox.cpp
//  mymath
//
//  Created by Loy Kyle Wong on 2019/8/8.
//  Copyright © 2019 Loy Kyle Wong. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <cmath>

void linEquSolve(float *mat, int n)
{
    int k, j, i;
    // to up triangular mat
    for(k = 0; k < n - 1; k++)
    {
        // find major row
        int major = k;
        for(j = k + 1; j < n; j++)
        {
            if(mat[j * (n + 1) + k] > mat[major * (n + 1) + k])
                major = j;
        }
        // swap major row
        if(major != k)
        {
            float t;
            for(i = 0; i < n + 1; i++)
            {
                t = mat[k * (n + 1) + i];
                mat[k * (n + 1) + i] = mat[major * (n + 1) + i];
                mat[major * (n + 1) + i] = t;
            }
        }
        // eliminating column k, form row k + 1 to n - 1
        for(j = k + 1; j < n; j++)
        {
            float c = mat[j * (n + 1) + k] / mat[k * (n + 1) + k];
            for(i = k; i < n + 1; i++)
            {
                mat[j * (n + 1) + i] -= mat[k * (n + 1) + i] * c;
            }
        }
    }
    // to 1
    for(k = 0; k < n; k++)
    {
        float c = mat[k * (n + 1) + k];
        for(i = k; i < n + 1; i++)
        {
            mat[k * (n + 1) + i] /= c;
        }
    }
    //
    for(k = n - 1; k >= 1; k--)
    {
        for(j = k - 1; j >= 0; j--)
        {
            float c = mat[j * (n + 1) + k];
            for(i = k; i < n + 1; i++)
            {
                mat[j * (n + 1) + i] -= mat[k * (n + 1) + i] * c;
            }
        }
    }
}

void Poly2ndApprox(float *x, float *y, int n, float *coef)
{
    int i, j, k;
    float equ[3][4] = {0.f};
    for(j = 0; j < 3; j++)
    {
        for(i = 0; i < 3; i++)
        {
            for(k = 0; k < n; k++)
                equ[j][i] += std::pow(x[k], (float)(i + j));
        }
        for(k = 0; k < n; k++)
            equ[j][i] += std::pow(x[k], (float)(j)) * y[k];
    }
    linEquSolve(&equ[0][0], 3);
    coef[0] = equ[0][3];
    coef[1] = equ[1][3];
    coef[2] = equ[2][3];
}

void LinApprox(float *x, float *y, int n, float *k, float *b)
{
    int i;
    float sumx = 0.0f, sumy = 0.0f, sumx2 = 0.0f, sumxy = 0.0f;
    float equ[2][3];
    for(i = 0; i < n; i++)
    {
        sumx += x[i];
        sumy += y[i];
        sumx2 += x[i] * x[i];
        sumxy += x[i] * y[i];
    }
    equ[0][0] = n;
    equ[0][1] = sumx;
    equ[0][2] = sumy;
    equ[1][0] = sumx;
    equ[1][1] = sumx2;
    equ[1][2] = sumxy;
    linEquSolve(&equ[0][0], 2);
    *b = equ[0][2];
    *k = equ[1][2];
}

void testLinApprox()
{
    float x[] = {1.f, 2.f, 3.f, 4.f, 5.f};
    float y1[] = {-.11f, -.21f, -.31f, -.41f, -.51f};   // y1 = -0.1*x-0.01
    float y2[] = {2.01f, 3.02f, 5.9f, 11.1f, 18.0f};      // y2 = x^2-2*x+3
    float k = 0.f, b = 0.f;
    float c[3] = {0.f};
    LinApprox(x, y1, 5, &k, &b);
    using namespace std;
    cout << "lin fit: k = " << k << ", b = " << b << ", ref is: k = -0.1, b = -0.01." << endl;
    Poly2ndApprox(x, y2, 5, c);
    printf("2nd poly fit: c0 = %f, c1 = %f, c2 = %f, ref is c0 = 3, c1 = -2, c2 = 1.", c[0], c[1], c[2]);
    cout << endl;
    
}
