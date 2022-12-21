#include <cstdint>
#include <iostream>
#include <random>
#include <chrono>

#include "../LoyMath/Cordic.h"

using namespace std;
using namespace LoyMath;

void cordic_u16_perf_test()
{
    const size_t w = 640, h = 480, f = 10;
    const size_t test_size = w * h * f;
    const int16_t max = (int16_t)((1UL << 14) * 0.707);
    // -------- prepare storage --------
    float   *fi  = new float  [test_size];
    float   *fq  = new float  [test_size];
    float   *fa  = new float  [test_size];
    int16_t *ii  = new int16_t[test_size];
    int16_t *iq  = new int16_t[test_size];
    int16_t *ia  = new int16_t[test_size];
    float   *iaf = new float  [test_size];
    // -------- random data set --------
    std::default_random_engine rnd_eng;
    std::uniform_int_distribution<int16_t> dist(1000, max);
    for(size_t i = 0; i < test_size; i++)
    {
        ii[i] = dist(rnd_eng);
        iq[i] = dist(rnd_eng);
        fi[i] = (float)ii[i];
        fq[i] = (float)iq[i];
    }
    // -------- test for std::atan2 --------
    auto tp0 = chrono::system_clock::now();
    for(size_t i = 0; i < test_size; i++)
    {
        fa[i] = std::atan2(fq[i], fi[i]);
    }
    auto tp1 = chrono::system_clock::now();
    double t0 = chrono::duration<double, std::milli>(tp1 - tp0).count();
    // -------- test for int16_t cordic --------
    tp0 = chrono::system_clock::now();
    for(size_t i = 0; i < test_size; i++)
    {
        ia[i] = Cordic<int16_t>::Atan2(iq[i], ii[i]);
    }
    tp1 = chrono::system_clock::now();
    double t1 = chrono::duration<double, std::milli>(tp1 - tp0).count();
    // -------- precision check --------
    float err_max = 0.0f, err_rms = 0.0f;
    int16_t err_max_i = 0, err_max_q = 0;
    for(size_t i = 0; i < test_size; i++)
    {
        iaf[i] = ia[i] / 32768.0f * (float)M_PI;
        float err = iaf[i] - fa[i];
        if(std::abs(err) > err_max)
        {
            err_max = std::abs(err);
            err_max_i = ii[i];
            err_max_q = iq[i];
        }
        err_rms += err * err;
    }
    err_rms = std::sqrt(err_rms / (float)test_size);
    printf("Time of %lu*%lu*%lu std::atan2(float)  : %7.3lfms.\r\n", w, h, f, t0);
    printf("Time of %lu*%lu*%lu Cordic<u16>::Atan2 : %7.3lfms.\r\n", w, h, f, t1);
    printf("Cordic::Atan2 RMSE      : %6.4fdeg.\r\n", err_rms * 180.f / (float)M_PI);
    printf("Cordic::Atan2 PeakE     : %6.4fdeg @atan2(%u, %u)\r\n", err_max * 180.f / (float)M_PI, err_max_i, err_max_q);

    delete[] fi;
    delete[] fq;
    delete[] fa;
    delete[] ii;
    delete[] iq;
    delete[] ia;
    delete[] iaf;
}

void testCordic()
{
    cordic_u16_perf_test();
}
