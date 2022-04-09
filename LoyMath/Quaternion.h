//
//  Quaternion.h
//  testQuaternion
//
//  Created by Loy Kyle Wong on 21/12/2016.
//  Copyright © 2016 Loy Kyle Wong. All rights reserved.
//

#ifndef Quaternion_h
#define Quaternion_h

#include <stdio.h>
#include <cmath>
#include "Matrix.h"

namespace LoyMath
{
    template<typename T=float>
    class Quat
    {
    private:
        T q, i, j, k;
    public:
        Quat()
        {
            q = 1.0; i = j = k = 0.0;
        }
//        Quat(T i, T j, T k)
//        {
//            q = (T)0.0; i = i; j = j; k = k;
//        }
        Quat(Quat<T> &q)
        {
            *this = q;
        }
        Quat(T q, T i, T j, T k)
        {
            this->q = q; this->i = i; this->j = j; this->k = k;
        }
        Quat(T q, Vec<3,T> &v)
        {
            this->q = q; i = v(0); j = v(1); k = v(2);
        }
        Quat(T ax, T ay, T az)
        {
            T an = std::sqrt(ax * ax + ay * ay + az * az);
            q = std::sin(an * (T)0.5);
            i = q * ax / an;
            j = q * ay / an;
            k = q * az / an;
            q = std::cos(an * (T)0.5);
        }
        Quat(Vec<3,T> &a)
        {
            T an = std::sqrt(a(0) * a(0) + a(1) * a(1) + a(2) * a(2));
            q = std::sin(an * (T)0.5);
            i = q * a(0) / an;
            j = q * a(1) / an;
            k = q * a(2) / an;
            q = std::cos(an * (T)0.5);
        }
        void Reset()
        {
            q = (T)1.0; i = j = k = (T)0.0;
        }
        
        Vec<3,T> Vector()
        {
            return Vec<3,T>(3, i, j, k);
        }
        
        Quat<T> Conj() const
        {
            Quat<T> rtn(q, -i, -j, -k);
            return rtn;
        }
        
        Quat<T>& operator+=(const Quat<T> &r)
        {
            q += r.q;
            i += r.i;
            j += r.j;
            k += r.k;
            return *this;
        }

        Quat<T>& operator-=(const Quat<T> &r)
        {
            q -= r.q;
            i -= r.i;
            j -= r.j;
            k -= r.k;
            return *this;
        }

        Quat<T>& RevMinus(const Quat<T> &r)
        {
            q = r.q - q;
            i = r.i - i;
            j = r.j - j;
            k = r.k - k;
            return *this;
        }
                
        Quat<T>& operator*=(const Quat<T> &r)
        {
            T lq = q, li = i, lj = j, lk = k;
            q = lq * r.q - li * r.i - lj * r.j - lk * r.k;
            i = lq * r.i + li * r.q + lj * r.k - lk * r.j;
            j = lq * r.j + lj * r.q + lk * r.i - li * r.k;
            k = lq * r.k + lk * r.q + li * r.j - lj * r.i;
            return *this;
        }
        
        Quat<T>& LMult(const Quat<T> &l)
        {
            T rq = q, ri = i, rj = j, rk = k;
            q = l.q * rq - l.i * ri - l.j * rj - l.k * rk;
            i = l.q * ri + l.i * rq + l.j * rk - l.k * rj;
            j = l.q * rj + l.j * rq + l.k * ri - l.i * rk;
            k = l.q * rk + l.k * rq + l.i * rj - l.j * ri;
            return *this;
        }
        
        Quat<T> operator*(const Quat<T> &r)
        {
            Quat<T> rtn;
            rtn.q = q * r.q - i * r.i - j * r.j - k * r.k;
            rtn.i = q * r.i + i * r.q + j * r.k - k * r.j;
            rtn.j = q * r.j + j * r.q + k * r.i - i * r.k;
            rtn.k = q * r.k + k * r.q + i * r.j - j * r.i;
            return rtn;
        }
        
        void ToAxisAndAngle(T &x, T &y, T &z, T &ang) const
        {
            T bs = (T)1.0 - q * q;
            T b = bs < (T)0.0 ? (T)0.0 : sqrt(bs);
            x = (i == (T)0.0 && b == (T)0.0) ? (T)0.0 : i / b;
            y = (j == (T)0.0 && b == (T)0.0) ? (T)0.0 : j / b;
            z = (k == (T)0.0 && b == (T)0.0) ? (T)0.0 : k / b;
            ang = (T)2.0 * atan2(b, q);
        }
        
        void ToTBAngleZYX(T &roll, T &pitch, T &yaw) const
        {
            T jsqr = j * j;
            
            roll = atan2((T)2.0 * (q * i + j * k), (T)1.0 - (T)2.0 * (i * i + jsqr));
            
            T t = (T)2.0 * (q * j - k * i);
            t = t > (T)1.0 ? (T)1.0 : t < (T)-1.0 ? (T)-1.0 : t;
            pitch = asin(t);
            
            yaw = atan2((T)2.0 * (q * k + i * j), (T)1.0 - (T)2.0 * (jsqr + k * k));
        }
//        void ToAtti(T &roll, T &pitch, T &yaw)
//        {
//            T qq = q * q, ii = i * i, jj = j * j, kk = k * k;
//            
//            roll = Atan2((T)2.0 * (q * i + j * k), qq - ii - jj + kk);
//            
//            T t = (T)2.0 * (q * j - k * i);
//            t = t > (T)1.0 ? (T)1.0 : t < (T)-1.0 ? (T)-1.0 : t;
//            pitch = Asin(t);
//            
//            yaw = Atan2((T)2.0 * (q * k + i * j), qq + ii - jj - kk);
//        }

    };
}


#endif /* Quaternion_h */
