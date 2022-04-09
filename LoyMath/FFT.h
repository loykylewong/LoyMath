/*
 * Fft.h
 *
 *  Created on: Oct 27, 2019
 *      Author: loywong
 */

#ifndef LOYMATH_FFT_H_
#define LOYMATH_FFT_H_

#include "Complex.h"
#include <cmath>

namespace LoyMath {

template <typename T = float, int M = 8>
class FFT
{
//    using namespace std;
private:
	const T pi = 3.1415926535897932384626;
	int len;
	Cplx<T> w[1 << (M - 1)];
	int bitRev(int data, int n)
	{
		const unsigned int lmask[] = {
				0x55555555, 0x33333333, 0x0f0f0f0f, 0x00ff00ff, 0x0000ffff};
		const unsigned int umask[] = {
				0xaaaaaaaa, 0xcccccccc, 0xf0f0f0f0, 0xff00ff00, 0xffff0000};
		int g, gi;
		for(g = 1, gi = 0; g < n; g <<= 1, gi++)
		{
			data = ((data & umask[gi]) >> g) | ((data & lmask[gi]) << g);
		}
		return data >> (g - n);
	}
    inline void swap(Cplx<T> &a, Cplx<T> &b)
    {
        Cplx<T> t = a;
        a = b;
        b = t;
    }
public:

	FFT(void)
	{
		this->len = 1 << M;
		for(int i = 0; i < (len >> 1); i++)
		{
			w[i].real = std::cos(2 * pi * (T)i / (T)len);
			w[i].imag = std::sin(2 * pi * (T)i / (T)len);
		}
	}
	// y can be the same as x
    void Trans(Cplx<T> (&x)[1<<M], Cplx<T> (&y)[1<<M], int a = -1, int b = 1)
	{
		Cplx<T> u, v;
		T gain = a == -1 ? (T)1.0 / (T)len :
				 a ==  0 ? (T)1.0 / std::sqrt((T)len) : (T)1.0;
		for(int i = 0; i < len; i++)
		{
			y[i] = x[i] * gain;
		}
		for(int step = 0, gLen = (len >> 1); step < M; step++, gLen >>= 1)
		{
			for(int grp = 0; grp < len; grp += (gLen << 1))
			{
				for(int i = grp, j = grp + gLen, k = 0;
						i < grp + gLen;
						i++, j++, k += (1 << step))
				{
					u = y[i] + y[j];
					if(b < 0)
						v = (y[i] - y[j]) * w[k].Conj();
					else
						v = (y[i] - y[j]) * w[k];
					y[i] = u;
					y[j] = v;
				}
			}
		}
		for(int i = 0; i < len; i++)
		{
			int j = bitRev(i, M);
			if(i < j)
				swap(y[i], y[j]);
		}
	}

};

}

#endif /* LOYMATH_FFT_H_ */
