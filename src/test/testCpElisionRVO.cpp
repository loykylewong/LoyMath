//
//  testCpElisionRVO.cpp
//  mymath
//
//  Created by Loy Kyle Wong on 2019/11/7.
//  Copyright © 2019 Loy Kyle Wong. All rights reserved.
//

#include <iostream>
#include "../LoyMath/Complex.h"
using namespace LoyMath;

struct CCC
{
    CCC() = default;
    CCC(const CCC&)
    {
        std::cout << "A copy was made.\n";
    }
};

CCC FFF()
{
    return CCC();
}

void testCERVO()
{
    std::cout << "Hello World!\n";
    CCC obj = FFF();
    
    Cplx<short> bbb(100, 100);
    Cplx<short> aaa = bbb.Conj();
}
