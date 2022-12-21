/*
 * Complex.h
 *
 *  Created on: Aug 30, 2016
 *      Author: Loywong Kyle Wong
 */

#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#include <cmath>

namespace LoyMath {

template <typename T = float>
class Cplx
{
public:
	T real;
	T imag;
    Cplx() = default;
    Cplx(T a, T b = (T)0) : real(a), imag(b) {}
	Cplx(T a, T b, char format)
	{
		if(format == 'p')	// polar format
		{
			this->real = a * std::cos(b);
			this->imag = a * std::sin(b);
		}
		else
		{
			this->real = a;
			this->imag = b;
		}
	}

	T Abs() const
	{
		return sqrt(this->real * this->real + this->imag * this->imag);
	}

	T Arg() const
	{
		return std::atan2(this->imag, this->real);
	}

	char *ToString(char *str, char format) const
	{
		if(format == 'p')	// polar format
			sprintf(str, "%g@%g", this->Abs(), this->Arg());
		else
		{
		    if(this->imag >= 0)
		        sprintf(str, "%g+%gj", this->real, this->imag);
		    else
                sprintf(str, "%g-%gj", this->real, -this->imag);
		}

		return str;
	}

	friend Cplx<T> operator+(const Cplx<T> &l, const Cplx<T> &r)
	{
		return Cplx<T>(l.real + r.real, l.imag + r.imag);
	}
    Cplx<T> &operator+=(const Cplx<T> &r)
    {
        this->real += r.real;
        this->imag += r.imag;
        return *this;
    }

	friend Cplx<T> operator-(const Cplx<T>& l, const Cplx<T>& r)
	{
		return Cplx<T>(l.real - r.real, l.imag - r.imag);
	}
    Cplx<T> &operator-=(const Cplx<T> &r)
    {
        this->real -= r.real;
        this->imag -= r.imag;
        return *this;
    }

	friend Cplx<T> operator*(const Cplx<T>& l, const Cplx<T>& r)
	{
		Cplx<T> rtn;
		rtn.real = l.real * r.real - l.imag * r.imag;
		rtn.imag = l.real * r.imag + l.imag * r.real;
		return rtn;
	}
    Cplx<T> &operator*=(const Cplx<T> &r)
    {
        Cplx<T> l = *this;
        this->real = l.real * r.real - l.imag * r.imag;
        this->imag = l.real * r.imag + l.imag * r.real;
        return *this;
    }

	friend Cplx<T> operator/(const Cplx<T>& l, const Cplx<T>& r)
	{
		T sqabs;
		Cplx<T> rtn;
		sqabs = r.real * r.real + r.imag * r.imag;
		rtn.real = (l.real * r.real + l.imag * r.imag) / sqabs;
		rtn.imag = (l.imag * r.real - l.real * r.imag) / sqabs;
		return rtn;
	}
    Cplx<T> &operator/=(const Cplx<T> &r)
    {
        T lr = real, li = imag;
        T sqabs = r.real * r.real + r.imag * r.imag;
        real = (lr * r.real + li * r.imag) / sqabs;
        imag = (li * r.real - lr * r.imag) / sqabs;
        return *this;
    }

	friend Cplx<T> operator^(const Cplx<T>& b, const Cplx<T>& e)
	{
		T babs, barg, abs, arg;
		babs = b.Abs();
		barg = b.Arg();
		abs = std::pow(babs, e.real) * std::exp(-e.imag * barg);
		arg = e.real * barg + e.imag * std::log(babs);
		return Cplx<T>(abs, arg, 'p');
	}
    Cplx<T> &operator^=(const Cplx<T> &e)
    {
        T babs, barg, abs, arg;
        babs = this->Abs();
        barg = this->Arg();
        abs = std::pow(babs, e.real) * std::exp(-e.imag * barg);
        arg = e.real * barg + e.imag * std::log(babs);
        this->real = abs * std::cos(arg);
        this->imag = abs * std::sin(arg);
        return *this;
    }

	Cplx<T> Conj(void)
	{
        return Cplx<T>(real, -imag);
	}

    Cplx<T> &ConjMe(void)
    {
        imag = -imag;
        return *this;
    }
};

template <typename T> Cplx<T> exp(const Cplx<T>& a)
{
	T abs;
	Cplx<T> rtn;
	abs = Exp(a.real);
	rtn.real = abs * Cos(a.imag);
	rtn.imag = abs * Sin(a.imag);
	return rtn;
}
template <typename T> Cplx<T> zpow(int i, T omega)
{
	T arg = omega * i;
	return Cplx<T>(std::cos(arg), std::sin(arg));
}
}

#endif

