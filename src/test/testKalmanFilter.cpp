//
//  testKalmanFilter.cpp
//  mymath
//
//  Created by Loy Kyle Wong on 10/01/2017.
//  Copyright © 2017 Loy Kyle Wong. All rights reserved.
//

//#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>

#include "../LoyMath/KalmanFilter.h"

using namespace std;
using namespace LoyMath;

const int KDLEN = 3501;
float *acc, *omg, *per, *ver, *pel, *vel;
Mat<4,1,float> *xout;

void testKalmanFilter()
{
    float x, dx, psi ,dpsi;
//    ifs >> acc >> omg >> per >> ver >> pel >> vel >> x >> dx >> psi >> dpsi;
//    ifs >> acc >> omg >> per >> ver >> pel >> vel >> x >> dx >> psi >> dpsi;
    acc = new float[KDLEN];
    omg = new float[KDLEN];
    per = new float[KDLEN];
    ver = new float[KDLEN];
    pel = new float[KDLEN];
    vel = new float[KDLEN];
    xout = new Mat<4, 1, float>[KDLEN];

    float Ts = 0.001f,
    W = 0.034f,
    ad = 0.35f,
    an = 0.004f,
    gd = 0.175f,
    gn = 8.7e-5f,
    res = 28.f, r = 0.0063f,
    pi = M_PI;
    float sa = an*an * 0.5f / Ts + ad*ad,
    sg = gn*gn * 0.5f / Ts + gd*gd,
    qx = 2.f * pi * r / res,
    qv = 2.f * pi * r / res / Ts;

    Mat<4,4,float> A(16,
                     1.f,  Ts, 0.f, 0.f,
                     0.f, 1.f, 0.f, 0.f,
                     0.f, 0.f, 1.f, 0.f,
                     0.f, 0.f, 0.f, 0.f);
    Mat<4,2,float> B(8,
                     .5f * Ts * Ts, 0.f,
                                Ts, 0.f,
                               0.f,  Ts,
                               0.f, 1.f);
    Mat<4,4,float> H(16,
                     1.f, 0.f,  W * .5f,      0.f,
                     0.f, 1.f,      0.f,  W * .5f,
                     1.f, 0.f, -W * .5f,      0.f,
                     0.f, 1.f,      0.f, -W * .5f);

    Mat<4,4,float> Q(16,
                     .25f*sa*Ts*Ts*Ts*Ts, .5f*sa*Ts*Ts*Ts, 0.f, 0.f,
                     .5f *sa*Ts*Ts*Ts   ,     sa*Ts*Ts   , 0.f, 0.f,
                     0.f, 0.f, sg*Ts*Ts, sg*Ts,
                     0.f, 0.f, sg*Ts   , sg
                     );

    Mat<4,4,float> R(16,
                     1.f/12.f*qx*qx, 0.f/12.f*qx*qv, 0.f, 0.f,
                     0.f/12.f*qx*qv, 1.f/12.f*qv*qv, 0.f, 0.f,
                     0.f, 0.f, 1.f/12.f*qx*qx, 0.f/12.f*qx*qv,
                     0.f, 0.f, 0.f/12.f*qx*qv, 1.f/12.f*qv*qv
                     );

    Mat<4,1,float> X0;
    Mat<4,4,float> P0;

    KalmanFilter<4, 2, 4, float> kf(A, B, H, Q, R, X0, P0);

    ifstream ifs("../test/kalman_data/kd");
    for(int i = 0; i < KDLEN; i++)
    {
        ifs >> acc[i] >> omg[i] >> per[i] >> ver[i] >> pel[i] >> vel[i] >> x >> dx >> psi >> dpsi;
    }
    ifs.close();

    auto tp0 = std::chrono::steady_clock::now();
    for(int i = 0; i < KDLEN; i++)
    {
        kf.Predict(Mat<2, 1, float>(2, acc[i], omg[i]));
        xout[i] = kf.Correct(Mat<4, 1, float>(4, per[i], ver[i], pel[i], vel[i]));
    }
    auto tp1 = std::chrono::steady_clock::now();
    float t = std::chrono::duration<float, std::micro>(tp1 - tp0).count();
    std::cout << "time cost: " << t << "us" << endl;

    ofstream ofs("../test/kalman_data/kdc");
    for(int i = 0; i < KDLEN; i++)
    {
        ofs << xout[i](0,0) << ' ' << xout[i](1,0) << ' ' << xout[i](2,0) << ' ' << xout[i](3,0) << endl;
    }
    ofs.close();

    delete[] acc ;
    delete[] omg ;
    delete[] per ;
    delete[] ver ;
    delete[] pel ;
    delete[] vel ;
    delete[] xout;
}
