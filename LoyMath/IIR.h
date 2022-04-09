//
//  IIR.h
//  mymath
//
//  Created by Loy Kyle Wong on 13/02/2017.
//  Copyright © 2017 Loy Kyle Wong. All rights reserved.
//

#ifndef IIR_h
#define IIR_h

#include <stdarg.h>

namespace LoyMath {
    template<typename T = float>
    class IIR1st
    {
    private:
        T num1;
        T den1;
        T gain;
        T z0;
    public:
        IIR1st()
        {
            num1 = (T)0;
            den1 = (T)0;
            gain = (T)0;
            z0 = (T)0;
        }
        IIR1st(T num1, T den1, T gain)
        {
            this->num1 = num1;
            this->den1 = den1;
            this->gain = gain;
            this->z0 = (T)0;
        }
        void UpdateCoef(T num1, T den1, T gain)
        {
            this->num1 = num1;
            this->den1 = den1;
            this->gain = gain;
        }
        T Filter(T data)
        {
            T rtn;
            rtn = z0 + data;
            z0 = (data * num1) - (rtn * den1);
            return rtn * gain;
        }
    };
    
    template<typename T = float>
    class IIR2nd
    {
    private:
        T num1, num2;
        T den1, den2;
        T gain;
        T z0, z1;
    public:
        IIR2nd()
        {
            num1 = (T)0;
            num2 = (T)0;
            den1 = (T)0;
            den2 = (T)0;
            gain = (T)0;
            z0 = (T)0;
            z1 = (T)0;
        }
        IIR2nd(T num1, T num2, T den1, T den2, T gain)
        {
            this->num1 = num1;
            this->num2 = num2;
            this->den1 = den1;
            this->den2 = den2;
            this->gain = gain;
            z0 = (T)0;
            z1 = (T)0;
        }
        void UpdateCoef(T num1, T num2, T den1, T den2, T gain)
        {
            this->num1 = num1;
            this->num2 = num2;
            this->den1 = den1;
            this->den2 = den2;
            this->gain = gain;
        }
        T Filter(T data)
        {
            T rtn;
            rtn = z0 + data;
            z0 = z1 + (data * num1) - (rtn * den1);
            z1 = (data * num2) - (rtn * den2);
            return rtn * gain;
        }
    };
    
    // Order must >= 3
    template<int Order, typename T = float>
    class IIRSos
    {
    private:
        static const int cNum = 3 * Order - (Order / 2);
        static const int sNum = (1 + Order) / 2;
        IIR2nd<T> stg[Order / 2];
        IIR1st<T> slast;
        T sdata[sNum];
    public:
        IIRSos()
        {
            for(int i = 0; i < Order; i++)
            {
                sdata[i] = (T)0;
            }
        }
        IIRSos(int num, ...)
        {
            va_list vl;
            va_start(vl, num);
            if(num != cNum)
                throw "Coefs number not met!";
            float c[cNum];
            for(int i = 0; i < cNum; i++)
            {
                c[i] = (T)va_arg(vl, double);
                if(--num == 0) break;
            }
            for(int o = 0, i = 0; o < Order; o+=2)
            {
                if(o != Order - 1)
                {
                    stg[o / 2].UpdateCoef(c[i], c[i+1], c[i+2], c[i+3], c[i+4]);
                    i += 5;
                }
                else
                {
                    slast.UpdateCoef(c[i], c[i+1], c[i+2]);
                }
            }
            va_end(vl);
            
            for(int i = 0; i < Order; i++)
            {
                sdata[i] = (T)0;
            }
        };
        IIRSos(T *coef)
        {
            for(int o = 0, i = 0; o < Order; o+=2)
            {
                if(o != Order - 1)
                {
                    stg[o / 2].UpdateCoef(coef[i], coef[i+1], coef[i+2], coef[i+3], coef[i+4]);
                    i += 5;
                }
                else
                {
                    slast.UpdateCoef(coef[i], coef[i+1], coef[i+2]);
                }
            }
        }
        T Filter(T data)
        {
            sdata[0] = data;
            int o;
            for(o = 0; o < Order - 2; o+=2)
            {
                sdata[o / 2 + 1] = stg[o / 2].Filter(sdata[o / 2]);
            }
            if(o == Order - 1) // slast, last 1st-order stg
            {
                return slast.Filter(sdata[o / 2]);
            }
            else // last 2nd-order stg
            {
                return stg[o / 2].Filter(sdata[o / 2]);
            }

        }
    };
}

#endif /* IIR_h */
