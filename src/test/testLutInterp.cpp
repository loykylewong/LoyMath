/*
 * testLutInterp.cpp
 *
 *  Created on: Dec 21, 2022
 *      Author: loywong
 */

#include <cstdint>
#include <iostream>
#include <chrono>
#include <cmath>

#include "../LoyMath/LutInterp.h"

using namespace std;
using namespace LoyMath;

void lut_interp_u16_test()
{
    const size_t w = 640, h = 480, f = 10;
    const size_t test_size = w * h * f;

    const uint8_t lut_aw = 8;
    const size_t lut_size = (size_t)1 << lut_aw;
    LutInterp<uint16_t, int16_t, lut_aw> li;
    for(size_t i = 0; i < lut_size + 1; i++)
    {
        // x: [0, 65536)->[0, 2pi); y: [-32767, 32767]->[-1, 1]
        li.Lut(i) = (int16_t)std::round(32767.0f * std::sin(2.0f * (float)M_PI * i / lut_size));
    }
    float af = 0.123f * (float)M_PI;
    uint16_t ai = std::round(af * 32768.f / (float)M_PI);
    cout << li(ai) << endl;
    cout << (float)li(ai) / 32767.f << endl;
    // -------- prepare storage --------
    float    *fx  = new float   [test_size];
    float    *fy  = new float   [test_size];
    uint16_t *ix  = new uint16_t[test_size];
    int16_t  *iy  = new int16_t [test_size];
    float    *iyf = new float   [test_size];
    // -------- data set --------
    for(size_t i = 0; i < test_size; i++)
    {
        ix[i] = i & 0xffff;
        fx[i] = (float)(int16_t)ix[i] / 32768.0f * (float)M_PI;
    }
    // -------- perf test for std::sin(float) --------
    auto tp0 = chrono::system_clock::now();
    for(size_t i = 0; i < test_size; i++)
    {
        fy[i] = std::sin(fx[i]);
    }
    auto tp1 = chrono::system_clock::now();
    double t0 = chrono::duration<double, std::milli>(tp1 - tp0).count();
    // -------- perf test for lutinterp --------
    tp0 = chrono::system_clock::now();
    for(size_t i = 0; i < test_size; i++)
    {
        iy[i] = li(ix[i]);
    }
    tp1 = chrono::system_clock::now();
    double t1 = chrono::duration<double, std::milli>(tp1 - tp0).count();
    // -------- precision check --------
    float err_max = 0.0f, err_rms = 0.0f;
    int16_t err_max_x = 0;
    for(size_t i = 0; i < test_size; i++)
    {
        iyf[i] = iy[i] / 32767.0f;
        float err = iyf[i] - fy[i];
        if(std::abs(err) > err_max)
        {
            err_max = std::abs(err);
            err_max_x = ix[i];
        }
        err_rms += err * err;
    }
    err_rms = std::sqrt(err_rms / (float)test_size);
    printf("Time of %lu*%lu*%lu std::sin(float)  : %7.3lfms.\r\n", w, h, f, t0);
    printf("Time of %lu*%lu*%lu Lut<u16,s16,%2u> : %7.3lfms.\r\n", w, h, f, lut_aw, t1);
    printf("RMSE  : %8.6f.\r\n", err_rms);
    printf("PeakE : %8.6f Lut(%d)\r\n", err_max, err_max_x);

    delete[] ix;
    delete[] iy;
    delete[] fx;
    delete[] fy;
    delete[] iyf;
}

void testLutInterp()
{
    lut_interp_u16_test();
}


