// this file is OBSOLETED.

/*
 * Mat.h
 *
 *  Created on: Oct 20, 2015
 *      Author: Loywong Kyle Wong
 */

#ifndef MAT_H_
#define MAT_H_

#include <stdarg.h>
#include <stdio.h>
#include "mathwrapper.h"

namespace LoyMath {
    
enum class MatError : int
{
    ArgOutOfRange       = 0x00010000,
    IndexOutOfRange     = 0x00020000,
    NotASquareMatrix    = 0x00030000,
    SingularMatrix      = 0x00040000
};

template <int R, int C, typename T = float>
class Mat
{
private:
//	int idx(int r, int c);
	inline void swap(T &a, T &b)
	{
		T t = a;
		a = b;
		b = t;
	}
    inline T cpSign(T a, T b)
    {
        return b >= 0.0 ? Abs(a) : -Abs(a);
    }
    inline T max(T a, T b)
    {
        return a >= b ? a : b;
    }
//    T pythag(T a, T b)
//    {
//        T at = Abs(a), bt = Abs(b), ct, result;
//        if(at > bt)
//        {
//            ct = bt / at;
//            result = at * Sqrt((T)1.0 + ct * ct);
//        }
//        else if(bt > 0.0)
//        {
//            ct = at / bt;
//            result = bt * Sqrt((T)1.0 + ct * ct);
//        }
//        else
//            result = 0.0;
//        return(result);
//    }
//	inline T abs(T a)
//	{
//		return a < 0 ? -a : a;
//	}
public:
    T elem[R][C] = {0};
	Mat(){}
	Mat(T diag)
	{
		for(int r = 0; r < R; r++)
			for(int c = 0; c < C; c++)
				if(r == c)
					elem[r][c] = diag;
				else
					elem[r][c] = 0;
	}
	Mat(const Mat<R,C,T> &m)
	{
		*this = m;
	}
	Mat(int num, ...)
	{
		va_list vl;
		va_start(vl, num);
		for(int r = 0; r < R; r++)
		{
			for(int c = 0; c < C; c++)
            {
				elem[r][c] = (T)va_arg(vl, double);
                if(--num == 0) {c = C; r = R; }
            }
		}
		va_end(vl);
	}
	virtual ~Mat(){}
	inline int Row(){return R;}
	inline int Col(){return C;}
    T &operator()(int row)
    {
        if(row >= R)
            throw MatError::ArgOutOfRange;
        else
            return elem[row][0];
    }
	T &operator()(int row, int col)
	{
		if(row >= R || col >= C)
            throw MatError::ArgOutOfRange;
		else
			return elem[row][col];
	}
	void SwapRow(int r1, int r2)
	{
		if(r1 < 0 || r1 >= R || r2 < 0 || r2 >= R)
            throw MatError::ArgOutOfRange;
		T tmp;
		for(int c = 0; c < C; c++)
		{
			tmp = elem[r1][c];
			elem[r1][c] = elem[r2][c];
			elem[r2][c] = tmp;
		}
	}
	// r1 += r2 * coef
	void RowTrans(int r1, int r2, T coef)
	{
		if(r1 < 0 || r1 >= R || r2 < 0 || r2 >= R)
            throw MatError::ArgOutOfRange;
		for(int c = 0; c < C; c++)
		{
			elem[r1][c] += elem[r2][c] * coef;
		}
	}
	template <int C2>
	Mat<R,C+C2,T> PadRight(const Mat<R,C2,T> &m)
	{
		Mat<R,C+C2,T> rtn;
		for(int r = 0; r < R; r++)
		{
			for(int c = 0; c < C; c++)
				rtn.elem[r][c] = elem[r][c];
			for(int c = 0; c < C2; c++)
				rtn.elem[r][C + c] = m.elem[r][c];
		}
		return rtn;
	}
	template <int R2, int C2>
	Mat<R2,C2,T> &SubMat(int startRow, int startCol, Mat<R2,C2,T> &m)
	{
		if(startRow + R2 > R || startCol + C2 > C)
            throw MatError::IndexOutOfRange;
		for(int r = 0; r < R2; r++)
			for(int c = 0; c < C2; c++)
			{
				m.elem[r][c] = elem[r + startRow][c + startCol];
			}
		return m;
	}
    Mat<R,C,T> &operator+=(T v)
    {
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                elem[r][c] += v;
            }
        return *this;
    }
	Mat<R,C,T> &operator+=(const Mat<R,C,T> &m)
	{
		for(int r = 0; r < R; r++)
			for(int c = 0; c < C; c++)
			{
				elem[r][c] += m.elem[r][c];
			}
		return *this;
	}
    Mat<R,C,T> &operator-=(T v)
    {
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                elem[r][c] -= v;
            }
        return *this;
    }
	Mat<R,C,T> &operator-=(const Mat<R,C,T> &m)
	{
		for(int r = 0; r < R; r++)
			for(int c = 0; c < C; c++)
			{
				elem[r][c] -= m.elem[r][c];
			}
		return *this;
	}
    Mat<R,C,T> &RevMinus(const Mat<R,C,T> &m)
    {
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                elem[r][c] = m.elem[r][c] - elem[r][c];
            }
        return *this;
    }
    Mat<R,C,T> operator+(const Mat<R,C,T> &m) const
    {
        Mat<R,C,T> rtn;
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                rtn.elem[r][c] = elem[r][c] + m.elem[r][c];
            }
        return rtn;
    }
    Mat<R,C,T> operator-(const Mat<R,C,T> &m) const
    {
        Mat<R,C,T> rtn;
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                rtn.elem[r][c] = elem[r][c] - m.elem[r][c];
            }
        return rtn;
    }
    Mat<R,C,T> &operator*=(T v)
    {
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                elem[r][c] *= v;
            }
        return *this;
    }
    // self = self * m
	Mat<R,C,T> &operator*=(const Mat<C,C,T> &m)
	{
		Mat<R,C,T> tmp(*this);
		for(int r = 0; r < R; r++)
			for(int c = 0; c < C; c++)
			{
				elem[r][c] = 0;
				for(int i = 0; i < C; i++)
				{
					elem[r][c] += tmp.elem[r][i] * m.elem[i][c];
				}
			}
		return *this;
	}
    // self = m * self
    Mat<R,C,T> &LMult(const Mat<R,R,T> &m)
    {
        Mat<R,C,T> tmp(*this);
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                elem[r][c] = 0;
                for(int i = 0; i < R; i++)
                {
                    elem[r][c] += m.elem[r][i] * tmp.elem[i][c];
                }
            }
        return *this;
    }
	template <int C2>
	Mat<R,C2,T> operator*(const Mat<C,C2,T> &m) const
	{
		Mat<R, C2, T> rtn;
		for(int r = 0; r < R; r++)
			for(int c = 0; c < C2; c++)
			{
				rtn.elem[r][c] = 0;
				for(int i = 0; i < C; i++)
				{
					rtn.elem[r][c] += elem[r][i] * m.elem[i][c];
				}
			}
		return rtn;
	}
    
	// return Transpose
	Mat<C,R,T> Trans() const
	{
		Mat<C,R,T> rtn;
		for(int r = 0; r < R; r++)
			for(int c = 0; c < C; c++)
			{
				rtn.elem[c][r] = elem[r][c];
			}
		return rtn;
	}
	// convert to transpose
	Mat<C,R,T> &TrMe()
	{
		if(R != C)
            throw MatError::NotASquareMatrix;
		for(int r = 1; r < R; r++)
			for(int c = 0; c < r; c++)
			{
				if(r != c)
				{
					swap(elem[r][c], elem[c][r]);
				}
			}
        return *this;
	}
	// inverse
    Mat<C,R,T> Inv()
    {
        if(C != R)
            throw MatError::NotASquareMatrix;
        else
        {
            Mat<C,R+R,T> dm;
            Mat<C,R,T> rtn;
            float coef;
            for(int r = 0; r < C; r++)
                for(int c = 0; c < R+R; c++)
                    if(c < R)
                        dm.elem[r][c] = elem[r][c];
                    else if(c - R == r)
                        dm.elem[r][c] = 1;
                    else
                        dm.elem[r][c] = 0;
            for(int r = 0; r < C - 1; r++)	// up triangle
            {
                int major = r;
                for(int i = r + 1; i < C; i++)	// find major row
                    if(Abs(dm.elem[i][r]) > Abs(dm.elem[major][r]))
                        major = i;
                if(major != r) // swap row
                    for(int c = 0; c < R+R; c++)
                        swap(dm.elem[r][c], dm.elem[major][c]);
                if(dm.elem[r][r] == (T)0)
                {
                    throw MatError::SingularMatrix;
                }
                for(int i = r + 1; i < C; i++)
                {
                    coef = -dm.elem[i][r] / dm.elem[r][r];
                    dm.elem[i][r] = 0;
                    for(int c = r + 1; c < R+R; c++)
                        dm.elem[i][c] += coef * dm.elem[r][c];
                }
            }
            for(int r = 0; r < C; r++)	// diag to 1
            {
                if(dm.elem[r][r] == (T)0)
                {
                    throw MatError::SingularMatrix;
                }
                coef = 1.f / dm.elem[r][r];
                dm.elem[r][r] = 1;
                for(int c = r + 1; c < R+R; c++)
                    dm.elem[r][c] *= coef;
            }
            for(int r = C - 1; r > 0; r--)	// to I
            {
                for(int i = r - 1; i >= 0; i--)
                {
                    float coef = -dm.elem[i][r];
                    for(int c = r; c < R+R; c++)
                        dm.elem[i][c] += dm.elem[r][c] * coef;
                }
            }
            dm.SubMat(0, R, rtn);
            return rtn;
        }
    }

//	Mat<R,C,T> Inv()
//	{
//		if(R != C)
//            throw MatError::NotASquareMatrix;
//		else
//		{
//			Mat<R,C+C,T> dm;
//            Mat<R,C,T> rtn;
//            float coef;
//			for(int r = 0; r < R; r++)
//				for(int c = 0; c < C+C; c++)
//					if(c < C)
//						dm.elem[r][c] = elem[r][c];
//					else if(c - C == r)
//						dm.elem[r][c] = 1;
//					else
//						dm.elem[r][c] = 0;
//			for(int r = 0; r < R - 1; r++)	// up triangle
//			{
//				int major = r;
//				for(int i = r + 1; i < R; i++)	// find major row
//					if(Abs(dm.elem[i][r]) > Abs(dm.elem[major][r]))
//						major = i;
//				if(major != r) // swap row
//					for(int c = 0; c < C+C; c++)
//						swap(dm.elem[r][c], dm.elem[major][c]);
//                if(dm.elem[r][r] == (T)0)
//                {
//                    throw MatError::SingularMatrix;
//                }
//				for(int i = r + 1; i < R; i++)
//				{
//					coef = -dm.elem[i][r] / dm.elem[r][r];
//                    dm.elem[i][r] = 0;
//                    for(int c = r + 1; c < C+C; c++)
//                        dm.elem[i][c] += coef * dm.elem[r][c];
//				}
//			}
//			for(int r = 0; r < R; r++)	// diag to 1
//			{
//                if(dm.elem[r][r] == (T)0)
//                {
//                    throw MatError::SingularMatrix;
//                }
//                coef = 1.f / dm.elem[r][r];
//                dm.elem[r][r] = 1;
//                for(int c = r + 1; c < C+C; c++)
//                    dm.elem[r][c] *= coef;
//			}
//			for(int r = R - 1; r > 0; r--)	// to I
//			{
//				for(int i = r - 1; i >= 0; i--)
//				{
//					float coef = -dm.elem[i][r];
//					for(int c = r; c < C+C; c++)
//						dm.elem[i][c] += dm.elem[r][c] * coef;
//				}
//			}
//			dm.SubMat(0, C, rtn);
//			return rtn;
//		}
//	}
    
    T NormM1() const
    {
        T rtn = 0.0;
        for(int r = 0; r < R; r++)
        {
            for(int c = 0; c < C; c++)
            {
                rtn += Abs(elem[r][c]);
            }
        }
        return rtn;
    }
    
    T NormMInf() const
    {
        T rtn = 0.0, a;
        for(int r = 0; r < R; r++)
        {
            for(int c = 0; c < C; c++)
            {
                a = Abs(elem[r][c]);
                if(a > rtn)
                    rtn = a;
            }
        }
        return rtn;
    }
    
    T Norm2() const
    {
        T rtn = 0.0;
        for(int r = 0; r < R; r++)
        {
            for(int c = 0; c < C; c++)
            {
                rtn += elem[r][c] * elem[r][c];
            }
        }
        return Sqrt(rtn);
    }
    
//    Mat<C,R,T> PInv(T eps = (T)0.0)
//    {
//        if(eps == (T)0.0)
//        {
//            if(sizeof(T) == sizeof(float))
//                eps = 1e-6;
//            else if(sizeof(T) == sizeof(double))
//                eps = 1e-14;
//        }
//        Mat<C,R,T> X(this->Trans());
//        Mat<R,R,T> I((T)1.0), E;
//        T nm = (*this * X).Norm2();
//        X *= (T)1.0 / nm;
//        nm = (T)1.0;
//        while(nm > eps)
//        {
////            nmd = nm * eps;
////            nmu = nm + nmd;
////            nmd = nm - nmd;
////            X *= II - (*this * X);
////            nm = X.NormMInf();
//            E = I - (*this * X);
//            nm = E.NormMInf();
//            E += I;
//            X *= E;
//        }
//        return X;
//    }

    Mat<C,R,T> PInv(T eps = (T)0.0)
    {
        try
        {
            return this->Inv();
        }
        catch (MatError err)
        {
            if(eps == (T)0.0)
            {
                if(sizeof(T) == sizeof(float))
                    eps = 1e-6;
                else if(sizeof(T) == sizeof(double))
                    eps = 1e-14;
            }
            Mat<R,R,T> II((T)2.0), temp;
            Mat<C,R,T> X(this->Trans());
            T nm = (*this * X).Norm2();
            X *= (T)1.0 / nm;
            do
            {
                temp = *this * X;
                X *= II - temp;
                nm = (*this - (temp * *this)).NormMInf();
            } while((nm > eps));
            return X;
        }
    }
    
	long PrintStr(char *str)
	{
		char *s = str;
		for(int r = 0; r < R; r++)
		{
			for(int c = 0; c < C; c++)
			{
				s += sprintf(s, "%g", this->elem[r][c]);
				if(c < C - 1) *s++ = '\t';
			}
			*s++ = '\n';
		}
		*s = '\0';
		return s - str;
	}
	char *ToStr(char *buf)
	{
		char *s = buf;
		for(int r = 0; r < R; r++)
		{
			for(int c = 0; c < C; c++)
			{
				s += sprintf(s, "%g", this->elem[r][c]);
				if(c < C - 1) *s++ = '\t';
			}
			*s++ = '\n';
		}
		*s = '\0';
		return buf;
	}
};

template<int R, typename T=float>
using Vec = Mat<R,1,T>;
//class Vec : Mat<R,1,T>
//{
//public:
//    T& operator()(int row)
//    {
//        if(row >= R)
//            throw MatError::ArgOutOfRange;
//        else
//            return this->elem[row][0];
//    }
//
//};
    
}

#endif /* MATF_H_ */
