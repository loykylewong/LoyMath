/*
 * Matrix.h
 *
 *  Created on: Aug 30, 2016
 *      Author: Loywong Kyle Wong
 */

#ifndef MAT_H_
#define MAT_H_

#include <stdarg.h>
//#include <stdio.h>
#include <string.h>
#include <sstream>
#include <cmath>

namespace LoyMath {
    
enum class MatError : int
{
    ArgOutOfRange       = 0x00010000,
    IndexOutOfRange     = 0x00020000,
    NotASquareMatrix    = 0x00030000,
    SingularMatrix      = 0x00040000,
    RowNumNotEqual3     = 0x00050000,
};

template <int R, int C, typename T = float>
class Mat
{
private:
//    int idx(int r, int c);
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
    template<typename T1>
    static constexpr T1 max(T1 a, T1 b)
    {
        return a >= b ? a : b;
    }
    template<typename T1>
    static constexpr T1 min(T1 a, T1 b)
    {
        return a <= b ? a : b;
    }
    // inverse 2 order
    inline void inv2(Mat<2,2,T> &m)
    {
        T den = -elem[0][1]*elem[1][0] + elem[0][0]*elem[1][1];
        if(den == (T)0.0)
            throw MatError::SingularMatrix;
        m.elem[0][0] = elem[1][1] / den;
        m.elem[0][1] = elem[0][1] / den;
        m.elem[1][0] = elem[1][0] / den;
        m.elem[1][1] = elem[0][0] / den;
    }
    // inverse 3 order
    inline void inv3(Mat<3,3,T> &m)
    {
        T den =
            -elem[0][2]*elem[1][1]*elem[2][0] + elem[0][1]*elem[1][2]*elem[2][0] +
             elem[0][2]*elem[1][0]*elem[2][1] - elem[0][0]*elem[1][2]*elem[2][1] -
             elem[0][1]*elem[1][0]*elem[2][2] + elem[0][0]*elem[1][1]*elem[2][2];
        if(den == (T)0.0)
            throw MatError::SingularMatrix;
        m.elem[0][0] = (-elem[1][2]*elem[2][1] + elem[1][1]*elem[2][2]) / den;
        m.elem[0][1] = ( elem[0][2]*elem[2][1] - elem[0][1]*elem[2][2]) / den;
        m.elem[0][2] = (-elem[0][2]*elem[1][1] + elem[0][1]*elem[1][2]) / den;
        m.elem[1][0] = ( elem[1][2]*elem[2][0] - elem[1][0]*elem[2][2]) / den;
        m.elem[1][1] = (-elem[0][2]*elem[2][0] + elem[0][0]*elem[2][2]) / den;
        m.elem[1][2] = ( elem[0][2]*elem[1][0] - elem[0][0]*elem[1][2]) / den;
        m.elem[2][0] = (-elem[1][1]*elem[2][0] + elem[1][0]*elem[2][1]) / den;
        m.elem[2][1] = ( elem[0][1]*elem[2][0] - elem[0][0]*elem[2][1]) / den;
        m.elem[2][2] = (-elem[0][1]*elem[1][0] + elem[0][0]*elem[1][1]) / den;
    }
    // inverse 4 order
    inline void inv4(Mat<4,4,T> &m)
    {
        T den =
            elem[0][3]*elem[1][2]*elem[2][1]*elem[3][0] - elem[0][2]*elem[1][3]*elem[2][1]*elem[3][0] -
            elem[0][3]*elem[1][1]*elem[2][2]*elem[3][0] + elem[0][1]*elem[1][3]*elem[2][2]*elem[3][0] +
            elem[0][2]*elem[1][1]*elem[2][3]*elem[3][0] - elem[0][1]*elem[1][2]*elem[2][3]*elem[3][0] -
            elem[0][3]*elem[1][2]*elem[2][0]*elem[3][1] + elem[0][2]*elem[1][3]*elem[2][0]*elem[3][1] +
            elem[0][3]*elem[1][0]*elem[2][2]*elem[3][1] - elem[0][0]*elem[1][3]*elem[2][2]*elem[3][1] -
            elem[0][2]*elem[1][0]*elem[2][3]*elem[3][1] + elem[0][0]*elem[1][2]*elem[2][3]*elem[3][1] +
            elem[0][3]*elem[1][1]*elem[2][0]*elem[3][2] - elem[0][1]*elem[1][3]*elem[2][0]*elem[3][2] -
            elem[0][3]*elem[1][0]*elem[2][1]*elem[3][2] + elem[0][0]*elem[1][3]*elem[2][1]*elem[3][2] +
            elem[0][1]*elem[1][0]*elem[2][3]*elem[3][2] - elem[0][0]*elem[1][1]*elem[2][3]*elem[3][2] -
            elem[0][2]*elem[1][1]*elem[2][0]*elem[3][3] + elem[0][1]*elem[1][2]*elem[2][0]*elem[3][3] +
            elem[0][2]*elem[1][0]*elem[2][1]*elem[3][3] - elem[0][0]*elem[1][2]*elem[2][1]*elem[3][3] -
            elem[0][1]*elem[1][0]*elem[2][2]*elem[3][3] + elem[0][0]*elem[1][1]*elem[2][2]*elem[3][3];
        if(den == (T)0.0)
            throw MatError::SingularMatrix;
        m.elem[0][0] =
            (-elem[1][3]*elem[2][2]*elem[3][1] + elem[1][2]*elem[2][3]*elem[3][1] +
              elem[1][3]*elem[2][1]*elem[3][2] - elem[1][1]*elem[2][3]*elem[3][2] -
              elem[1][2]*elem[2][1]*elem[3][3] + elem[1][1]*elem[2][2]*elem[3][3]) / den;
        m.elem[0][1] =
            ( elem[0][3]*elem[2][2]*elem[3][1] - elem[0][2]*elem[2][3]*elem[3][1] -
              elem[0][3]*elem[2][1]*elem[3][2] + elem[0][1]*elem[2][3]*elem[3][2] +
              elem[0][2]*elem[2][1]*elem[3][3] - elem[0][1]*elem[2][2]*elem[3][3]) / den;
        m.elem[0][2] =
            (-elem[0][3]*elem[1][2]*elem[3][1] + elem[0][2]*elem[1][3]*elem[3][1] +
              elem[0][3]*elem[1][1]*elem[3][2] - elem[0][1]*elem[1][3]*elem[3][2] -
              elem[0][2]*elem[1][1]*elem[3][3] + elem[0][1]*elem[1][2]*elem[3][3]) /den;
        m.elem[0][3] =
            ( elem[0][3]*elem[1][2]*elem[2][1] - elem[0][2]*elem[1][3]*elem[2][1] -
              elem[0][3]*elem[1][1]*elem[2][2] + elem[0][1]*elem[1][3]*elem[2][2] +
              elem[0][2]*elem[1][1]*elem[2][3] - elem[0][1]*elem[1][2]*elem[2][3]) / den;
        m.elem[1][0] =
            ( elem[1][3]*elem[2][2]*elem[3][0] - elem[1][2]*elem[2][3]*elem[3][0] -
              elem[1][3]*elem[2][0]*elem[3][2] + elem[1][0]*elem[2][3]*elem[3][2] +
              elem[1][2]*elem[2][0]*elem[3][3] - elem[1][0]*elem[2][2]*elem[3][3]) / den;
        m.elem[1][1] =
            (-elem[0][3]*elem[2][2]*elem[3][0] + elem[0][2]*elem[2][3]*elem[3][0] +
              elem[0][3]*elem[2][0]*elem[3][2] - elem[0][0]*elem[2][3]*elem[3][2] -
              elem[0][2]*elem[2][0]*elem[3][3] + elem[0][0]*elem[2][2]*elem[3][3]) / den;
        m.elem[1][2] =
            ( elem[0][3]*elem[1][2]*elem[3][0] - elem[0][2]*elem[1][3]*elem[3][0] -
              elem[0][3]*elem[1][0]*elem[3][2] + elem[0][0]*elem[1][3]*elem[3][2] +
              elem[0][2]*elem[1][0]*elem[3][3] - elem[0][0]*elem[1][2]*elem[3][3]) / den;
        m.elem[1][3] =
            (-elem[0][3]*elem[1][2]*elem[2][0] + elem[0][2]*elem[1][3]*elem[2][0] +
              elem[0][3]*elem[1][0]*elem[2][2] - elem[0][0]*elem[1][3]*elem[2][2] -
              elem[0][2]*elem[1][0]*elem[2][3] + elem[0][0]*elem[1][2]*elem[2][3]) / den;
        m.elem[2][0] =
            (-elem[1][3]*elem[2][1]*elem[3][0] + elem[1][1]*elem[2][3]*elem[3][0] +
              elem[1][3]*elem[2][0]*elem[3][1] - elem[1][0]*elem[2][3]*elem[3][1] -
              elem[1][1]*elem[2][0]*elem[3][3] + elem[1][0]*elem[2][1]*elem[3][3]) / den;
        m.elem[2][1] =
            ( elem[0][3]*elem[2][1]*elem[3][0] - elem[0][1]*elem[2][3]*elem[3][0] -
              elem[0][3]*elem[2][0]*elem[3][1] + elem[0][0]*elem[2][3]*elem[3][1] +
              elem[0][1]*elem[2][0]*elem[3][3] - elem[0][0]*elem[2][1]*elem[3][3]) / den;
        m.elem[2][2] =
            (-elem[0][3]*elem[1][1]*elem[3][0] + elem[0][1]*elem[1][3]*elem[3][0] +
              elem[0][3]*elem[1][0]*elem[3][1] - elem[0][0]*elem[1][3]*elem[3][1] -
              elem[0][1]*elem[1][0]*elem[3][3] + elem[0][0]*elem[1][1]*elem[3][3]) / den;
        m.elem[2][3] =
            ( elem[0][3]*elem[1][1]*elem[2][0] - elem[0][1]*elem[1][3]*elem[2][0] -
              elem[0][3]*elem[1][0]*elem[2][1] + elem[0][0]*elem[1][3]*elem[2][1] +
              elem[0][1]*elem[1][0]*elem[2][3] - elem[0][0]*elem[1][1]*elem[2][3]) / den;
        m.elem[3][0] =
            ( elem[1][2]*elem[2][1]*elem[3][0] - elem[1][1]*elem[2][2]*elem[3][0] -
              elem[1][2]*elem[2][0]*elem[3][1] + elem[1][0]*elem[2][2]*elem[3][1] +
              elem[1][1]*elem[2][0]*elem[3][2] - elem[1][0]*elem[2][1]*elem[3][2]) / den;
        m.elem[3][1] =
            (-elem[0][2]*elem[2][1]*elem[3][0] + elem[0][1]*elem[2][2]*elem[3][0] +
              elem[0][2]*elem[2][0]*elem[3][1] - elem[0][0]*elem[2][2]*elem[3][1] -
              elem[0][1]*elem[2][0]*elem[3][2] + elem[0][0]*elem[2][1]*elem[3][2]) / den;
        m.elem[3][2] =
            ( elem[0][2]*elem[1][1]*elem[3][0] - elem[0][1]*elem[1][2]*elem[3][0] -
              elem[0][2]*elem[1][0]*elem[3][1] + elem[0][0]*elem[1][2]*elem[3][1] +
              elem[0][1]*elem[1][0]*elem[3][2] - elem[0][0]*elem[1][1]*elem[3][2]) / den;
        m.elem[3][3] =
            (-elem[0][2]*elem[1][1]*elem[2][0] + elem[0][1]*elem[1][2]*elem[2][0] +
              elem[0][2]*elem[1][0]*elem[2][1] - elem[0][0]*elem[1][2]*elem[2][1] -
              elem[0][1]*elem[1][0]*elem[2][2] + elem[0][0]*elem[1][1]*elem[2][2]) / den;
    }
public:
    T elem[R][C];
    Mat()
    {
        memset(&(elem[0][0]), 0, R * C * sizeof(T));
    }
    Mat(T *elem_array)
    {
        memcpy(&(elem[0][0]), elem_array, R * C * sizeof(T));
    }
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
    void ResetAll(T v = (T)0)
    {
        for(int r = 0; r < R; r++)
        {
            for(int c = 0; c < C; c++)
            {
                elem[r][c] = v;
            }
        }
    }
    void ResetAll(T *elem_array)
    {
        memcpy(&(elem[0][0]), elem_array, R * C * sizeof(T));
    }
    void ResetToDiag(T diag = (T)0)
    {
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
                if(r == c)
                    elem[r][c] = diag;
                else
                    elem[r][c] = (T)0;
    }
    void ResetDiag(T diag = (T)0)
    {
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
                if(r == c)
                    elem[r][c] = diag;
    }
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
    Mat<R-1,C-1,T> CycAdjoint(int row, int col)
    {
        Mat<R-1,C-1,T> rtn;
        for(int r = row + 1, i = 0; i < R-1; r++, i++)
        {
            for(int c = col + 1, j = 0; j < C-1; c++, j++)
            {
                rtn.elem[i][j] = elem[r % R][c % C];
            }
        }
        return rtn;
    }
    Mat<R-1,C-1,T> Adjoint(int row, int col, int& sign)
    {
        Mat<R-1,C-1,T> rtn;
        sign = 1;
        for(int r = 0, i = 0; i < R-1; r++, i++)
        {
            if(r == row)
                r++;
            for(int c = 0, j = 0; j < C-1; c++, j++)
            {
                if(c == col)
                    c++;
                rtn.elem[i][j] = elem[r][c];
            }
        }
        for(int r = row + 1; r < R; r++)
            sign = -sign;
        for(int c = col + 1; c < R; c++)
            sign = -sign;
        return rtn;
    }
    template <int R1, int C1, typename T1>
    Mat<R,C,T> &CopyFrom(Mat<R1,C1,T1> &m, int srcR, int srcC, int dstR, int dstC, int rows = min(R, R1), int cols = min(C, C1))
    {
        for(int sr = srcR, dr = dstR; dr < dstR + rows; sr++, dr++)
            for(int sc = srcC, dc = dstC; dc < dstC + cols; sc++, dc++)
            {
                elem[dr % R][dc % C] = m.elem[sr % R1][sc % C1];
            }
        return *this;
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
    Mat<R,C,T> &RevMinus(const T &v)
    {
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                elem[r][c] = v - elem[r][c];
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
    
    Mat<R,C,T> operator*(T v) const
    {
        Mat<R,C,T> rtn;
        for(int r = 0; r < R; r++)
            for(int c = 0; c < C; c++)
            {
                rtn.elem[r][c] = elem[r][c] * v;
            }
        return rtn;
    }

    Mat<R,C,T> Cross(const Mat<R,C,T> &r) const
    {
        Mat<R,C,T> rtn;
        if(R != 3)
            throw MatError::RowNumNotEqual3;
        for(int c = 0; c < C; c++)
        {
            rtn.elem[0][c] = elem[1][c] * r.elem[2][c] - elem[2][c] * r.elem[1][c];
            rtn.elem[1][c] = elem[2][c] * r.elem[0][c] - elem[0][c] * r.elem[2][c];
            rtn.elem[2][c] = elem[0][c] * r.elem[1][c] - elem[1][c] * r.elem[0][c];
        }
        return rtn;
    }

    Mat<R,C,T> LCross(const Mat<R,C,T> &l) const
    {
        Mat<R,C,T> rtn;
        if(R != 3)
            throw MatError::RowNumNotEqual3;
        for(int c = 0; c < C; c++)
        {
            rtn.elem[0][c] = l.elem[1][c] * elem[2][c] - l.elem[2][c] * elem[1][c];
            rtn.elem[1][c] = l.elem[2][c] * elem[0][c] - l.elem[0][c] * elem[2][c];
            rtn.elem[2][c] = l.elem[0][c] * elem[1][c] - l.elem[1][c] * elem[0][c];
        }
        return rtn;
    }

    friend Mat<R,C,T> Cross(const Mat<R,C,T> &l, const Mat<R,C,T> &r)
    {
        Mat<R,C,T> rtn;
        if(R != 3)
            throw MatError::RowNumNotEqual3;
        for(int c = 0; c < C; c++)
        {
            rtn.elem[0][c] = l.elem[1][c] * r.elem[2][c] - l.elem[2][c] * r.elem[1][c];
            rtn.elem[1][c] = l.elem[2][c] * r.elem[0][c] - l.elem[0][c] * r.elem[2][c];
            rtn.elem[2][c] = l.elem[0][c] * r.elem[1][c] - l.elem[1][c] * r.elem[0][c];
        }
        return rtn;
    }
    
    Mat<1,C,T> Dots(const Mat<R,C,T> &r) const
    {
        Mat<1,C,T> rtn;
        for(int i = 0; i < C; i++)
        {
            rtn.elem[0][i] = (T)0;
            for(int j = 0; j < R; j++)
                rtn.elem[0][i] += elem[j][i] * r.elem[j][i];
        }
        return rtn;
    }
    
    friend Mat<1,C,T> Dots(const Mat<R,C,T> &l, const Mat<R,C,T> &r)
    {
        Mat<1,C,T> rtn;
        for(int i = 0; i < C; i++)
        {
            rtn.elem[0][i] = (T)0;
            for(int j = 0; j < R; j++)
                rtn.elem[0][i] += l.elem[j][i] * r.elem[j][i];
        }
        return rtn;
    }
    
    T Dot(const Mat<R,C,T> &r) const
    {
        T rtn;
        rtn = (T)0;
        for(int i = 0; i < C; i++)
            for(int j = 0; j < R; j++)
                rtn += elem[j][i] * r.elem[j][i];
        return rtn;
    }
    
    friend T Dot(const Mat<R,C,T> &l, const Mat<R,C,T> &r)
    {
        T rtn;
        rtn = (T)0;
        for(int i = 0; i < C; i++)
            for(int j = 0; j < R; j++)
                rtn += l.elem[j][i] * r.elem[j][i];
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
    T Det()
    {
        if(C != R)
            throw MatError::NotASquareMatrix;
#if(C == 1)
        else if(C == 1)
        {
            return elem[0][0];
        }
#elif(C == 2)
        else if(C == 2)
        {
            T rtn;
            rtn = elem[0][0]*elem[1][1] - elem[0][1]*elem[1][0];
            return rtn;
        }
#elif(C == 3)
        else if(C == 3)
        {
            T rtn;
            rtn =   elem[0][0]*elem[1][1]*elem[2][2] - elem[0][2]*elem[1][1]*elem[2][0]
                  + elem[0][1]*elem[1][2]*elem[2][0] - elem[0][0]*elem[1][2]*elem[2][1]
                  + elem[0][2]*elem[1][0]*elem[2][1] - elem[0][1]*elem[1][0]*elem[2][2]
            return rtn;
        }
#endif
        else
        {
            Mat<R,C,T> dm(*this);
            T coef;
            int swsign = 1;
            for(int r = 0; r < C - 1; r++)    // up triangle
            {
                int major = r;
                for(int i = r + 1; i < C; i++)    // find major row
                    if(abs(dm.elem[i][r]) > abs(dm.elem[major][r]))
                        major = i;
                if(major != r) // swap row
                {
                    for(int c = 0; c < R; c++)
                    {
                        swap(dm.elem[r][c], dm.elem[major][c]);
                    }
                    swsign = -swsign;
                }
                if(dm.elem[r][r] == (T)0)
                {
                    return (T)0.0;
                }
                for(int i = r + 1; i < C; i++)
                {
                    coef = -dm.elem[i][r] / dm.elem[r][r];
                    dm.elem[i][r] = (T)0.0;
                    for(int c = r + 1; c < R; c++)
                        dm.elem[i][c] += coef * dm.elem[r][c];
                }
            }
            coef = (T)swsign;
            for(int i = 0; i < R; i++)
            {
                coef *= dm.elem[i][i];
            }
            return coef;
        }
    }
    // inverse
    Mat<C,R,T> Inv()
    {
        if(C != R)
            throw MatError::NotASquareMatrix;
#if(C == 1)
        else if(C == 1)
        {
            return Mat<1,1,T>((T)1.0 / elem[0][0]);
        }
#elif(C == 2)
        else if(C == 2)
        {
            Mat<C,R,T> rtn;
            this->inv2(rtn);
            return rtn;
        }
#elif(C == 3)
        else if(C == 3)
        {
            Mat<C,R,T> rtn;
            this->inv3(rtn);
            return rtn;
        }
#elif(C == 4)
        else if(C == 4)
        {
            Mat<C,R,T> rtn;
            this->inv4(rtn);
            return rtn;
        }
#endif
        else
        {
            Mat<C,R+R,T> dm;
            Mat<C,R,T> rtn;
            T coef;
            for(int r = 0; r < C; r++)
                for(int c = 0; c < R+R; c++)
                    if(c < R)
                        dm.elem[r][c] = elem[r][c];
                    else if(c - R == r)
                        dm.elem[r][c] = 1;
                    else
                        dm.elem[r][c] = 0;
            for(int r = 0; r < C - 1; r++)    // up triangle
            {
                int major = r;
                for(int i = r + 1; i < C; i++)    // find major row
                    if(abs(dm.elem[i][r]) > abs(dm.elem[major][r]))
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
                    dm.elem[i][r] = (T)0.0;
                    for(int c = r + 1; c < R+R; c++)
                        dm.elem[i][c] += coef * dm.elem[r][c];
                }
            }
            for(int r = 0; r < C; r++)    // diag to 1
            {
                if(dm.elem[r][r] == (T)0.0)
                {
                    throw MatError::SingularMatrix;
                }
                coef = 1.f / dm.elem[r][r];
                dm.elem[r][r] = 1;
                for(int c = r + 1; c < R+R; c++)
                    dm.elem[r][c] *= coef;
            }
            for(int r = C - 1; r > 0; r--)    // to I
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
                a = abs(elem[r][c]);
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
        return std::sqrt(rtn);
    }

    // pseudo inverse use iteration algor, it's very slow,
    // you may use Inv(A.Trans() * A) * A.Trans() for left pinv,
    // and A.Trans() * Inv(A * A.Trans()) for right pinv.
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
    
    long PrintStr(char *str, const char *format = "%g")
    {
        char *s = str;
        for(int r = 0; r < R; r++)
        {
            for(int c = 0; c < C; c++)
            {
                s += sprintf(s, format, this->elem[r][c]);
                if(c < C - 1)
                {
                    *s++ = ',';
                    *s++ = ' ';
                }
            }
            *s++ = ';';
            *s++ = '\n';
        }
        *s = '\0';
        return s - str;
    }
    char *ToStr(char *buf, const char *format = "%g")
    {
        char *s = buf;
        for(int r = 0; r < R; r++)
        {
            for(int c = 0; c < C; c++)
            {
                s += sprintf(s, format, this->elem[r][c]);
                if(c < C - 1)
                    {
                        *s++ = ',';
                        *s++ = ' ';
                    }
                }
            *s++ = ';';
            *s++ = '\n';
        }
        *s = '\0';
        return buf;
    }
    std::string ToStr()
    {
        using namespace std;
        stringstream rtn;
        for(int r = 0; r < R; r++)
        {
            for(int c = 0; c < C; c++)
            {
                rtn << this->elem[r][c];
                if(c < C - 1) rtn << ", ";
            }
            rtn << ';' << endl;
        }
        return rtn.str();
    }
};

template<int R, typename T=float>
using Vec = Mat<R,1,T>;
    
}

#endif /* MATF_H_ */
