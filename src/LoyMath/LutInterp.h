/*
 * LutInterp.h
 *
 *  Created on: Dec 20, 2022
 *      Author: loywong
 */

#ifndef SRC_LUTINTERP_H_
#define SRC_LUTINTERP_H_

#include <type_traits>
#include <cstdint>
#include <cassert>

namespace LoyMath
{

template<typename XT = uint16_t, typename YT = int16_t, typename = std::enable_if_t< std::is_arithmetic_v<YT> &&
    (std::is_same_v<XT, uint8_t> || std::is_same_v<XT, uint16_t> || std::is_same_v<XT, uint32_t>) &&
    sizeof(YT) >= sizeof(XT) >>
class LutInterp
{
private:
    uint8_t lut_aw;
    uint8_t lut_fw;
    XT lut_fmask;
    size_t lut_size;
    YT *lut;
public:
    LutInterp() = delete;
    LutInterp(const LutInterp &) = delete;
    LutInterp(size_t lut_aw) : lut_aw(lut_aw)
    {
        assert(lut_aw >= 4 && lut_aw <= sizeof(XT) * 8);
        lut_fw = 8 * sizeof(XT) - lut_aw;
        lut_fmask = (1UL << lut_fw) - 1;
        lut_size = 1 + ((size_t)1 << lut_aw);  // extend one for the last interpolation domain
        lut = new YT[lut_size];
    }
    virtual ~LutInterp()
    {
        delete[] lut;
    }
    // use for custom initialize the lut
    // CAUTION: actual LUT size is (1 << lut_aw) + 1, the last one
    //          is extended for the last interpolation domain, and
    //          MUST be initialized.
    YT &Lut(size_t idx)
    {
        assert(idx < lut_size);
        return lut[idx];
    }
    YT operator()(XT x) const
    {
        if(!lut_fw)
        {
            return lut[x];
        }
        else
        {
            XT idx = x >> lut_fw;
            YT y0 = lut[idx];
            YT y1 = lut[idx + 1];
            x &= lut_fmask;
            YT w1 = (YT)x;
            YT w0 = (YT)(lut_fmask - x + 1);
            if constexpr(std::is_same_v<YT, uint8_t> || std::is_same_v<YT, int8_t>)
            {
                return (YT)(((int16_t)y0 * w0 + (int16_t)y1 * w1) >> lut_fw);
            }
            else if constexpr(std::is_same_v<YT, uint16_t> || std::is_same_v<YT, int16_t>)
            {
                return (YT)(((int32_t)y0 * w0 + (int32_t)y1 * w1) >> lut_fw);
            }
            else if constexpr(std::is_same_v<YT, uint32_t> || std::is_same_v<YT, int32_t>)
            {
                return (YT)(((int64_t)y0 * w0 + (int64_t)y1 * w1) >> lut_fw);
            }
            else if constexpr(std::is_floating_point_v<YT>)
            {
                return (y0 * w0 + y1 * w1) / (YT)(1 << lut_fw);
            }
        }
    }
};

}

#endif /* SRC_LUTINTERP_H_ */
