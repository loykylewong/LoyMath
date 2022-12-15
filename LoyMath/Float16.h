#ifndef __FLOAT16_H__
#define __FLOAT16_H__

#include <cinttypes>
#include <type_traits>
#include <limits>

namespace LoyMath
{

// T : base type, all calculation are done by base type.
template<typename T = float, typename = std::enable_if_t<
    std::is_same_v<T, float> /*|| std::is_same_v<T, double>*/ > >
class fp16
{
private:
    union _fp16
    {
        struct
        {
            uint16_t man : 10;
            uint16_t exp :  5;
            uint16_t sig :  1;
        };
        uint16_t bits;
    } value;
    union _fp32
    {
        struct
        {
            uint32_t man : 23;
            uint32_t exp :  8;
            uint32_t sig :  1;
        };
        uint32_t bits;
        float value;
        _fp32() : bits(0) {}
        _fp32(float v) : value(v) {}
        int16_t sexp() const { return (int16_t)exp - 127; }
    };
    union _fp64
    {
        struct
        {
            uint64_t man : 52;
            uint64_t exp : 11;
            uint64_t sig :  1;
        };
        uint64_t bits;
        double value;
        _fp64() : bits(0) {}
        _fp64(double v) : value(v) {}
        int16_t sexp() const { return (int16_t)exp - 1023; }
    };
    template<typename U = T, typename = std::enable_if_t<std::is_same_v<U, float>>>
    _fp16 from_base(T value) const
    {
        _fp16 f16;
        _fp32 f32(value);
        uint16_t man = 1 + (f32.man >> 12); // 11 bit for further rounding
        // ---- round, may cause: ----
        //   1. subnormal to_base normal
        //   2. normal to_base overflow
        int16_t exp = (int16_t)f32.exp - (127 - 15) + ((man == 0x800) ? 1 : 0);
        man = (man == 0x800) ? 0x0 : man >> 1;
        if(exp <= 0)  // zero or subnormal
        {
            f16.sig = f32.sig;
            f16.exp = 0;
            f16.man = (0x400 | man) >> (1 - exp);
        }
        else if(exp < 31) // normal
        {
            f16.sig = f32.sig;
            f16.exp = exp;
            f16.man = man;
        }
        else if(exp >= 31 && f32.exp < 255) // overflow
        {
            f16.sig = f32.sig;
            f16.exp = 31;
            f16.man = 0;
        }
        else    // inf, nan and indet
        {
            f16.sig = f32.sig;
            f16.exp = 31;
            f16.man = (f32.man >> 13) | (f32.man & 1);  // keep man for inf, nan, snan and ndef
        }
        return f16;
    }
    template<typename U = T, typename = std::enable_if_t<std::is_same_v<U, float>>>
    T to_base(const _fp16 &f16) const
    {
    	_fp32 f32;
        uint32_t man = (uint32_t)f16.man << 13;
        int16_t exp = (int16_t)f16.exp + (127 - 15);
        if(f16.exp == 0) // zero or subnormal
        {
        	if(f16.man == 0)	// zero
        	{
        		f32.sig = f16.sig;
        		f32.exp = 0;
        		f32.man = 0;
        	}
        	else 				// subnormal to normal
        	{
        		while(!(man & (1UL << 23)))
        		{
        			man <<= 1;
        			exp--;
        		}
        		f32.sig = f16.sig;
        		f32.exp = exp + 1;
        		f32.man = man & 0x7fffff;
        	}
        }
        else if(f16.exp < 31)	// normal
        {
        	f32.sig = f16.sig;
        	f32.exp = exp;
        	f32.man = man;
        }
        else					// inf, nan and indet
        {
        	f32.sig = f16.sig;
        	f32.exp = 255;
        	f32.man = man;
        }
        return f32.value;
    }
    fp16(const _fp16 &f16) : value(f16) {}
public:
    // ---- constructors and assign ----
    // default constructor
    fp16() : value(0) {}
    // construct from base type
    fp16(T value)
    {
        this->value = from_base(value);    
    }
    // assign from base type
    fp16 &operator=(T value)
    {
        this->value = from_base(value);
    }
    // assign from another
    fp16 &operator=(const fp16& r)
    {
        this->value = r.value;
    }
    // ---- converter----
    // convert to base type
    T to_base() const
    {
        return to_base(value);
    }
    T to_base2() const
    {
    	if(value.exp < 31)
    	{
    		int16_t e = value.exp == 0 ? -14 : (int16_t)value.exp - 15;
    		int32_t ex = 1 << (e < 0 ? -e : e);
    		T man = ((T)(value.exp != 0) + (T)value.man / (T)1024);
    		T exp = e < 0 ? (T)1 / (T)ex : (T)ex;
    		T val = value.sig ? -man * exp : man * exp;
    		return val;
    	}
    	else
    	{
    		if(value.man == 0)
    		{
    			return value.sig ? -std::numeric_limits<T>::infinity() :
    					            std::numeric_limits<T>::infinity() ;
    		}
    		else if(value.man < 512)
    		{
    			return value.sig ? -std::numeric_limits<T>::signaling_NaN() :
    					            std::numeric_limits<T>::signaling_NaN() ;
    		}
    		else if(value.man > 512)
    		{
    			return value.sig ? -std::numeric_limits<T>::quiet_NaN() :
    					            std::numeric_limits<T>::quiet_NaN() ;
    		}
    		else
    		{
    			union { float f; uint32_t bits; };
    			bits = (1UL << 31) | (255UL << 23) | (1UL << 22);
    			return value.sig ? f : std::numeric_limits<T>::quiet_NaN();
    		}
    	}
    }
    // ---- implicit converter ----
    // convert operator to base type
    operator T() const
    {
        return to_base(value);
    }
    // ---- status check ----
    bool isnormal()
    {
        return value.exp > 0 && value.exp < 31;
    }
    bool issub()
    {
        return value.exp == 0 && value.man != 0;
    }
    bool isinf()
    {
        return value.exp == 31 && value.man == 0;
    }
    bool isnan()
    {
        return value.exp == 31 && value.man != 0;
    }
    bool issnan()
    {
        return value.exp == 31 && value.man < 512;
    }
    bool isqnan()
    {
        return value.exp == 31 && value.man >= 512;
    }
    bool isindet()
    {
        return value.sig == 1 && value.exp == 31 && value.man == 512;
    }
    bool iszero()
    {
        return value.exp == 0 && value.man == 0;
    }
    bool ispos()
    {
        return value.sig == 0;
    }
    bool isneg()
    {
        return value.sig != 0;
    }
    int8_t sign()
    {
        return value.sig ? -1 : 1;
    }
    // ---- compares ----
    bool operator==(const fp16& r)
    {
        return to_base(value) == to_base(r);
    }
    bool operator!=(const fp16& r)
    {
        return to_base(value) != to_base(r);
    }
    bool operator<(const fp16& r)
    {
        return to_base(value) < to_base(r);
    }
    bool operator<=(const fp16& r)
    {
        return to_base(value) <= to_base(r);
    }
    bool operator>=(const fp16& r)
    {
        return to_base(value) >= to_base(r);
    }
    bool operator>(const fp16& r)
    {
        return to_base(value) > to_base(r);
    }
    // ---- assign arith oprs ----
    fp16 &operator+=(const fp16& r)
    {
        this->value = from_base(to_base(value) + to_base(r));
    }
    fp16 &operator-=(const fp16& r)
    {
        this->value = from_base(to_base(value) - to_base(r));
    }
    fp16 &operator*=(const fp16& r)
    {
        this->value = from_base(to_base(value) * to_base(r));
    }
    fp16 &operator/=(const fp16& r)
    {
        this->value = from_base(to_base(value) / to_base(r));
    }
    // ---- arith oprs ----
    fp16 operator+(const fp16& r) const
    {
        return fp16(from_base(to_base(value) + to_base(r)));
    }
    fp16 operator-(const fp16& r) const
    {
        return fp16(from_base(to_base(value) - to_base(r)));
    }
    fp16 operator*(const fp16& r) const
    {
        return fp16(from_base(to_base(value) * to_base(r)));
    }
    fp16 operator/(const fp16& r) const
    {
        return fp16(from_base(to_base(value) / to_base(r)));
    }
};

}

#endif
