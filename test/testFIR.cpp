//
//  testFIR.cpp
//  mymath
//
//  Created by Loy Kyle Wong on 11/02/2017.
//  Copyright © 2017 Loy Kyle Wong. All rights reserved.
//

#include <stdio.h>
#include <cmath>
#include "../LoyMath/FIR.h"
using namespace LoyMath;

const int dlen = 4096;
float wave[dlen];
float wout[dlen];

void testFIR()
{
    FIRSymm<15> sfir(8,
                     0.001734049016867023063070973876165226102f,
                     0.008416143153479254540583198718195490073f,
                     0.010871617236390278965485478579466871452f,
                     -0.016010360359646268885169817508540290873f,
                     -0.056770464972300163641882875253941165283f,
                     -0.012048177265742469888021304313951986842f,
                     0.180164024639846431785628055877168662846f,
                     0.38360027088864467881279551875195465982f
                   );
    
    float omg;
    float ang = 0.f;
    // initialize wave
    for(int i = 0; i < dlen; i++)
    {
        // freq from 0 to M_PI
        omg = (float)M_PI * (float)i / (float)dlen;
        ang += omg;
        wave[i] = sinf(ang);
    }
    for(int i = 0; i < dlen; i++)
    {
        wout[i] = sfir.Filter(wave[i]); // debug & watch
    }

}

