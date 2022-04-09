//
//  main.cpp
//  testQuaternion
//
//  Created by Loy Kyle Wong on 21/12/2016.
//  Copyright © 2016 Loy Kyle Wong. All rights reserved.
//

#include <iostream>
#include "../LoyMath/Quaternion.h"
#include <cmath>
using namespace LoyMath;

int testQuaternion()
{
    float r, p, y;
    float rl, pl, yl;
    float ux, uy, uz, zeta;
    float uxl, uyl, uzl, zetal;
    
    Quat<> q;
    Quat<> ql;
    
    // q:  45deg about y'
    // ql: 45deg about y
    for(int i = 0; i < 10; i++)
    {
        q *= Quat<>(0.f, M_PI_4 / 10, 0.f);
        ql.LMult(Quat<>(0.f, M_PI_4 / 10, 0.f));
        q.ToTBAngleZYX(r, p, y);
        ql.ToTBAngleZYX(rl, pl, yl);
    }
    q.ToAxisAndAngle(ux, uy, uz, zeta);
    ql.ToAxisAndAngle(uxl, uyl, uzl, zetal);

    // q:  90deg about z'
    // ql: 90deg about z
    for(int i = 0; i < 10; i++)
    {
        q *= Quat<>(0.f, 0.f, M_PI_2 / 10);
        ql.LMult(Quat<>(0.f, 0.f, M_PI_2 / 10));
        q.ToTBAngleZYX(r, p, y);
        ql.ToTBAngleZYX(rl, pl, yl);
    }
    q.ToAxisAndAngle(ux, uy, uz, zeta);
    ql.ToAxisAndAngle(uxl, uyl, uzl, zetal);
    
    Vec<3> v(3, 1.0, 0.0, 0.0);
    // rotate (1,0,0)
    Vec<3> rotatedVec = (q * Quat<>(0, v) * q.Conj()).Vector();
    // (1,0,0) in rotated frame to coordinates in origin frame
    Vec<3> coorInNED = (q * Quat<>(0, v) * q.Conj()).Vector();
    // (1,0,0) in origin frame to coordinates in rotated frame
    Vec<3> coorInBody = (q.Conj() * Quat<>(0, v) * q).Vector();
    return 0;
}
