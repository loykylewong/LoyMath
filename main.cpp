//
//  main.cpp
//  loymaths
//
//  Created by Loy Kyle Wong on 2021/11/12.
//

#include <iostream>

void testComplex();
void testMatrix();
void testLinApprox();
void testKalmanFilter();
void testQuaternion();
void testFFT();
void testFIR();
void testIIR();

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    testMatrix();
    testComplex();
    testLinApprox();
    testKalmanFilter();
    testQuaternion();
    testFFT();
    testFIR();
    testIIR();
    return 0;
}
