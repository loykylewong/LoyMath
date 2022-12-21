/*
 * Cordic.h
 *
 *  Created on: Dec 19, 2022
 *      Author: loywong
 */

#ifndef SRC_CORDIC_H__
#define SRC_CORDIC_H__

#include <type_traits>
#include <cstdint>
#include <array>
#include <cmath>

namespace LoyMath
{

// fixed point CORDIC for calculating triangular functions
// a are in [-pi, pi), mapped as [-128, 127) for int8_t, [-32768, 32767) for int16_t.
// \sqrt{x^2 + y^2} must be <= 77 for int8_t, <= 19897 for int16_t to prevent from overflow.
template<typename T = int16_t, typename = std::enable_if_t<
    std::is_same_v<T, int8_t> || std::is_same_v<T, int16_t> || std::is_same_v<T, int32_t>>>
struct Cordic
{
private:
    static constexpr size_t dw = 8 * sizeof(T);
    static constexpr size_t n_iter = dw - 2;
    static constexpr T lam = std::round(0.6072529350 * (double)(1UL << (dw - 1)));
    static constexpr auto thetas{ []() constexpr{
        std::array<T, n_iter> rtn{};
        double e2 = 1.0;
        for(size_t i = 0; i < n_iter; i++)
        {
            rtn[i] = std::round(std::atan(e2) / M_PI * (double)(1UL << (dw - 1)));
            e2 *= 0.5;
        }
        return rtn;
    }() };
    inline static void quad_trans(T &x, T &y, T &a)
    {
        x = -x;
        y = -y;
        if constexpr(std::is_same_v<T, int8_t>)
            a ^= 0x80;
        else if constexpr(std::is_same_v<T, int16_t>)
            a ^= 0x8000;
        else
            a ^= 0x80000000;
    }
    inline static void rot_quad_trans(T &x, T &y, T &a)
    {
        bool is_in_quad23;
        if constexpr(std::is_same_v<T, int8_t>)
            is_in_quad23 = (a & 0xc0      ) == 0x80       || (a & 0xc0      ) == 0x40      ;
        else if constexpr(std::is_same_v<T, int16_t>)
            is_in_quad23 = (a & 0xc000    ) == 0x8000     || (a & 0xc000    ) == 0x4000    ;
        else
            is_in_quad23 = (a & 0xc0000000) == 0x80000000 || (a & 0xc0000000) == 0x40000000;
        if(is_in_quad23)
            quad_trans(x, y, a);
    }
    inline static void vec_quad_trans(T &x, T &y, T &a)
    {
        bool is_in_quad23;
        if constexpr(std::is_same_v<T, int8_t>)
            is_in_quad23 = x & 0x80;
        else if constexpr(std::is_same_v<T, int16_t>)
            is_in_quad23 = x & 0x8000;
        else
            is_in_quad23 = x & 0x80000000;
        if(is_in_quad23)
        {
            quad_trans(x, y, a);
        }
    }
    inline static void scale(T &x, T &y)
    {
        if constexpr(std::is_same_v<T, int8_t>)
        {
            x = (T)(((int16_t)x * (int16_t)lam) >> 7);
            y = (T)(((int16_t)y * (int16_t)lam) >> 7);
        }
        else if constexpr(std::is_same_v<T, int16_t>)
        {
            x = (T)(((int32_t)x * (int32_t)lam) >> 15);
            y = (T)(((int32_t)y * (int32_t)lam) >> 15);
        }
        else
        {
            x = (T)(((int64_t)x * (int64_t)lam) >> 31);
            y = (T)(((int64_t)y * (int64_t)lam) >> 31);
        }
    }
public:
    static void Rotate(T &x, T &y, T &a)
    {
        T xs, ys;
        rot_quad_trans(x, y, a);
        for(size_t i = 0; i < n_iter; i++)
        {
            xs = x >> i;
            ys = y >> i;
            if(a > 0)
            {
                a -= thetas[i];
                x -= ys;
                y += xs;
            }
            else
            {
                a += thetas[i];
                x += ys;
                y -= xs;
            }
        }
        scale(x, y);
    }
    static void Vector(T &x, T &y, T &a)
    {
        T xs, ys;
        vec_quad_trans(x, y, a);
        for(size_t i = 0; i < n_iter; i++)
        {
            xs = x >> i;
            ys = y >> i;
            if(y < 0)
            {
                a -= thetas[i];
                x -= ys;
                y += xs;
            }
            else
            {
                a += thetas[i];
                x += ys;
                y -= xs;
            }
        }
        scale(x, y);
    }
    static T Sin(T ang, T one = 1 << (dw - 2))
    {
        T x = one, y = 0;
        Rotate(x, y, ang);
        return y;
    }
    static float Sin(float ang)
    {
        T x = 1 << (dw - 2);
        T y = 0;
        T a = std::round(ang / M_PI * (1UL << (dw - 1)));
        Rotate(x, y, a);
        return (float)y / (float)(1 << (dw - 2));
    }
    static T Cos(T ang, T one = 1 << (dw - 2))
    {
        T x = one, y = 0;
        Rotate(x, y, ang);
        return x;
    }
    static float Cos(float ang)
    {
        T x = 1 << (dw - 2);
        T y = 0;
        T a = std::round(ang / M_PI * (1UL << (dw - 1)));
        Rotate(x, y, a);
        return (float)x / (float)(1 << (dw - 2));
    }
    static void SinCos(T ang, T &sin, T &cos, T one = 1 << (dw - 2))
    {
        cos = one;
        sin = 0;
        Rotate(cos, sin, ang);
    }
    static void SinCos(float ang, float &s, float &c)
    {
        T x = 1 << (dw - 2);
        T y = 0;
        T a = std::round(ang / M_PI * (1UL << (dw - 1)));
        Rotate(x, y, a);
        s = (float)y / (float)(1 << (dw - 2));
        c = (float)x / (float)(1 << (dw - 2));
    }
    // \sqrt{i^2 + q^2} must be less than (1 << (dw - 1)) * 0.6072529350
    static T Atan2(T q, T i)
    {
        T a = 0;
        Vector(i, q, a);
        return a;
    }
    static float Atan2(float q, float i)
    {
        T x, y, a = 0;
        constexpr T one = std::round((1 << (dw - 2)) * std::sqrt(0.5f));
        x = i > q ? one : std::round(i / q * one);
        y = q > i ? one : std::round(q / i * one);
        Vector(x, y, a);
        return a * (float)M_PI / (float)(1 << (dw - 1));
    }
};

}

#endif
