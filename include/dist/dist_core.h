#ifndef DIST_CORE_H__
#define DIST_CORE_H__
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include <cstdint>
#include <cmath>

namespace dist {
using std::size_t;
using std::uint64_t;



struct NotImplementedError: public std::runtime_error {
     template<typename...Args>NotImplementedError(Args &&...args):
        std::runtime_error(std::forward<Args>(args)...) {}
};

// See http://www.iiisci.org/journal/CV$/sci/pdfs/GS315JG.pdf

template<typename FT=double>
struct SetDistInterface {
    FT operator()(size_t a, size_t b, size_t c) const {
        throw NotImplementedError("Not written");
    }
};

#define DEFINE_SETSIM_FUNC(NAME, RETVAL) \
    template<typename FT=double>\
    struct NAME##Similarity: public SetDistInterface<FT> {\
        constexpr FT operator()(size_t a, size_t b, size_t c) const {\
            return (RETVAL);\
        }\
    };\
    template<typename FT=double>\
    struct NAME##Distance: public SetDistInterface<FT> {\
        constexpr FT operator()(size_t a, size_t b, size_t c) const {\
            return 1. - (RETVAL);\
        }\
    };\
    template<typename FT> using NAME##Sim = NAME##Similarity<FT>;\
    template<typename FT> using NAME##Dist = NAME##Distance<FT>;
    
#define DEFINE_SETDIST_FUNC(NAME, RETVAL) \
    template<typename FT=double>\
    struct NAME##Distance: public SetDistInterface<FT> {\
        constexpr FT operator()(size_t a, size_t b, size_t c) const {\
            return (RETVAL);\
        }\
    };\
    template<typename FT> using NAME##Dist = NAME##Distance<FT>;

DEFINE_SETSIM_FUNC(Jaccard, static_cast<FT>(a) / (a+b+c)) // Also equivalent to Tanimoto
DEFINE_SETSIM_FUNC(Dice, static_cast<FT>(2*a) / (2 * a + b + c))
DEFINE_SETSIM_FUNC(ThreeW, static_cast<FT>(3*a) / (3 * a + b + c))
DEFINE_SETSIM_FUNC(SokalSneath, static_cast<FT>(a) / (a + 2*(b + c)))
DEFINE_SETSIM_FUNC(Cosine, static_cast<FT>(a) / std::sqrt((a+b) * (a+c)))
DEFINE_SETSIM_FUNC(OchiaiI, FT(a) / ((a+b)*(a+c)))
DEFINE_SETSIM_FUNC(SorgenFrei, a*a / static_cast<FT>((a+b)*(a+c)))
DEFINE_SETSIM_FUNC(MountFord, FT(a) / (.5 * (a * b + a * c) + b * c))
DEFINE_SETSIM_FUNC(McConnaughey, (a*a - b * c) / ((a + b) * (b + c)))
DEFINE_SETSIM_FUNC(Otsuka, a / std::sqrt((a+b)*(a+c)))
DEFINE_SETSIM_FUNC(KulczynskiII, .5 * a * (2*a + b + c)  / ((a + b) * (a + c)))
DEFINE_SETSIM_FUNC(DriverKroeber, .5 * a / (1./(a+b) + 1./(a + c)))
DEFINE_SETSIM_FUNC(Johnson, FT(a) / (1./(a+b) + 1./(a + c)))
DEFINE_SETSIM_FUNC(Simpson, FT(a) / std::min(a+b, a+c))
DEFINE_SETSIM_FUNC(BraunBauquet, FT(a) / std::max(a+b, a+c))


DEFINE_SETDIST_FUNC(Hamming, b+c)
DEFINE_SETDIST_FUNC(Euclid, std::sqrt(b+c))
DEFINE_SETDIST_FUNC(LanceWilliams, (b+c)/FT(2*a+b+c)) // Also Bray & Curtis
DEFINE_SETDIST_FUNC(Hellinger, 2. * std::sqrt(1. - FT(a) / ((a+b) * (a+c))))
DEFINE_SETDIST_FUNC(Chord, std::sqrt(2. * (1. - FT(a) / ((a+b) * (a+c)))))

template<typename FT=double>
struct MinkowskiDistance: public SetDistInterface<FT> {
    const FT pinv_;
    constexpr MinkowskiDistance(FT p=1.): pinv_(1./p) {}
    constexpr FT operator()(size_t a, size_t b, size_t c) const {
        return std::pow(b+c, pinv_);
    }
};
template<typename FT=double>
struct SoftJaccardSimilarity: public SetDistInterface<FT> {
    const FT smooth_;
    constexpr SoftJaccardSimilarity(FT p=1e-6): smooth_(p) {}
    constexpr FT operator()(size_t a, size_t b, size_t c) const {
        return FT(a) / (a+b+c + smooth_);
    }
};




template<typename FT> using CzekanowskiSimilarity = DiceSimilarity<FT>;
template<typename FT> using CzekanowskiSim = CzekanowskiSimilarity<FT>;
template<typename FT> using SquaredEuclidDist = HammingDist<FT>;
template<typename FT> using Canberra = HammingDist<FT>;
} // dist

#endif  /* DIST_CORE_H__ */
