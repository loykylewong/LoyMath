//
//  KalmanFilter.h
//  mymath
//
//  Created by Loy Kyle Wong on 10/01/2017.
//  Copyright © 2017 Loy Kyle Wong. All rights reserved.
//

#ifndef KalmanFilter_h
#define KalmanFilter_h

#include "Matrix.h"

using namespace std;

namespace LoyMath {
    
    // N: state number; L: input number; M: measurement number
    template<int N, int L, int M, typename T = float>
    class KalmanFilter
    {
    private:
        Mat<N,N,T> A;
        Mat<N,L,T> B;
        Mat<M,N,T> H;
        Mat<N,1,T> xpri, xpos;
        Mat<N,N,T> Ppri, Ppos;
        Mat<N,M,T> K;
        Mat<N,N,T> AT;
        Mat<N,M,T> HT;
    public:
        Mat<N,N,T> Q;
        Mat<M,M,T> R;
        
        KalmanFilter(const Mat<N,N,T> &A, const Mat<N,L,T> &B, const Mat<M,N,T> &H,
                     const Mat<N,N,T> &Q, const Mat<M,M,T> &R,
                     const Mat<N,1,T> &X0, const Mat<N,N,T> &P0)
        {
            this->A = A; this->B = B, this->H = H;
            this->Q = Q; this->R = R;
            this->xpos = X0; this->Ppos = P0;
            AT = A.Trans();
            HT = H.Trans();
        }
        
        void Predict(const Mat<L,1,T> &u)
        {
            // x- = A * x + B * u
            xpri = A * xpos;
            xpri += B * u;
            // P- = A P A^T + Q
            Ppri = A * Ppos * AT;
            Ppri += Q;
        }
        
        Mat<N,1,T>& Correct(const Mat<M,1,T> &z)
        {
            // K = P- * H^T * (H * P- * H^T + R)^-1
            Mat<M,M,T> r;
            Mat<M,1,T> zk = z;
            Mat<N,N,T> i((T)1.0);
            r = H * Ppri * HT;
            r += R;
            K = Ppri * HT * r.PInv((T)1e-5);
            // x = x- + K * (z - H * x-)
            zk -= H * xpri;
            xpos = xpri;
            xpos += K * zk;
            // P = (I - K * H) * P-
            i -= K * H;
            Ppos = i * Ppri;
            return xpos;
        }
    };
    
}

#endif /* KalmanFilter_h */
