/*
 * testFFT.cpp
 *
 *  Created on: Oct 27, 2019
 *      Author: Loywong Kyle Wong
 */


#include "../LoyMath/FFT.h"

using namespace LoyMath;

void testFFT(void)
{
	FFT<float, 8> fft;
	Cplx<float> x[256], y[256], z[256];
	for(int i = 0; i < 256; i++)
	{
		x[i].real = i < 128 ? 1.0f : -1.0f;
		x[i].imag = 0.f;
	}
	fft.Trans(x, y);	// watch result in watch...
	fft.Trans(y, z, 1, -1);
	fft.Trans(x, y, 0, 1);
	fft.Trans(x, x);
}

