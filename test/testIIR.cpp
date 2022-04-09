//
//  testIIR.cpp
//  mymath
//
//  Created by Loy Kyle Wong on 13/02/2017.
//  Copyright © 2017 Loy Kyle Wong. All rights reserved.
//

#include <cmath>
#include "../LoyMath/IIR.h"
using namespace LoyMath;

const double pi = 3.1415926535897932384626;

void testIIR()
{
    const int dlen = 4096;
    float *wave = new float[dlen];
    float *wout = new float[dlen];
    IIRSos<9, float> iir(23,
                  -1.942480944379443741354407393373548984528,
                  1.0,
                  -1.971465158132082340358692817972041666508,
                  0.989254183046181401195440230367239564657,
                  0.098852005821808253127613852484500966966,

                  -1.928601969568733309756680682767182588577,
                  1.0,
                  -1.952775102227502834750794136198237538338,
                  0.967322422577252871711550596955930814147,
                  0.203749602921530881260991918679792433977,

                  -1.879323229886710855396358965663239359856,
                  1.0,
                  -1.935962702964093562485459187882952392101,
                  0.945158787255304000396449737309012562037,
                  0.07620426269759740522946600549403228797,

                  -1.622582557211260434471000735356938093901,
                  1.0,
                  -1.922991305196560363199864696071017533541,
                  0.926922042208547458663758789043640717864,
                  0.010414826042333949138174453707961220061,

                  1.0,
                  -0.958960676763685238022105750133050605655,
                  0.064198533136987256941807800103561021388
                  );
//    IIR1st<> iir(
//                 1,
//                 -0.766838340633112180988462114328285679221,
//                 0.116580829683443965016920174093684181571);

    float omg;
    float ang = 0.f;
    // initialize wave
    for(int i = 0; i < dlen; i++)
    {
        // freq from 0 to M_PI
        omg = (float)pi * (float)i / (float)dlen;
        ang += omg;
        wave[i] = std::sin(ang);
    }
    for(int i = 0; i < dlen; i++)
    {
        wout[i] = iir.Filter(wave[i]); // debug & watch
    }

    delete[] wave;
    delete[] wout;
}
