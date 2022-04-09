//
//  FIR.h
//  mymath
//
//  Created by Loy Kyle Wong on 11/02/2017.
//  Copyright © 2017 Loy Kyle Wong. All rights reserved.
//

#ifndef FIR_h
#define FIR_h

#include <stdarg.h>

namespace LoyMath {
    
    template <int Order, typename T = float>
    class FIRDirect
    {
    private:
        static const int cNum = Order + 1;
        T coef[cNum];
        T prod[cNum];
        T z[cNum];
    public:
        FIRDirect()
    		{
            for(int i = 0; i < cNum; i++)
            {
            		this->z[i] = (T)0;
            }
    		}
        FIRDirect(int num, ...)
        {
        		if(num != cNum)
        			throw "Coef num not met!\n";
            va_list vl;
            va_start(vl, num);
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = (T)va_arg(vl, double);
            }
            for(int i = 0; i < cNum; i++)
            {
            		this->z[i] = (T)0;
            }
            va_end(vl);
        };
        FIRDirect(float *coef)
        {
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = coef[i];
                this->z[i] = (T)0;
            }
        }
        void UpdateCoef(int num, ...)
        {
			if(num != cNum)
				throw "Coef num not met!\n";
            va_list vl;
            va_start(vl, num);
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = (T)va_arg(vl, double);
                if(--num == 0) break;
            }
            va_end(vl);
        };
        void UpdateCoef(float *coef)
        {
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = coef[i];
            }
        }
        T Filter(T data)
        {
            for(int i = 0; i < cNum; i++)
                prod[i] = data * coef[i];
            
            for(int i = Order; i >= 1; i--)
            {
                z[i] = z[i - 1] + prod[Order - i];
            }
            z[0] = prod[Order];
            return z[Order];
        }
    };
    
    template <int Order, typename T = float>
    class FIRSymm
    {
    private:
        static const int cNum = (Order + 2) / 2;
        T coef[cNum];
        T prod[cNum];
        T z[Order + 1];
        
        inline int min(int a, int b)
        {
            return a < b ? a : b;
        }
    public:
        FIRSymm()
        {
            for(int i = 0; i < Order + 1; i++)
            {
            		this->z[i] = (T)0;
            }
        }
        FIRSymm(int num, ...)
        {
        		if(num != cNum)
        			throw "Coef num not met!\n";
            va_list vl;
            va_start(vl, num);
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = (T)va_arg(vl, double);
            }
            for(int i = 0; i < Order + 1; i++)
            {
            		this->z[i] = (T)0;
            }
            va_end(vl);
        };
        FIRSymm(float *coef)
        {
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = coef[i];
            }
            for(int i = 0; i < Order + 1; i++)
            {
            		this->z[i] = (T)0;
            }
        }
        void UpdateCoef(int num, ...)
        {
            va_list vl;
            va_start(vl, num);
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = (T)va_arg(vl, double);
                if(--num == 0) break;
            }
            va_end(vl);
        };
        void UpdateCoef(float *coef)
        {
            for(int i = 0; i < cNum; i++)
            {
                this->coef[i] = coef[i];
            }
        }
        T Filter(T data)
        {
            for(int i = 0; i < cNum; i++)
                prod[i] = data * coef[i];
            
            for(int i = Order; i >= 1; i--)
            {
                z[i] = z[i - 1] + prod[min(Order - i, i)];
            }
            z[0] = prod[0];
            return z[Order];
        }
    };
}

#endif /* FIR_h */
