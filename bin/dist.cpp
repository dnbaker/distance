#include "dist/dist_core.h"
#include <cstdio>

#define XTR(x) STRINGI(x)
#define STRINGI(x) #x

int main() {
    size_t a = 14, b = 12, c = 15;
    double v = dist::MinkowskiDistance<>()(a,b,c);
    std::fprintf(stderr, "v: %lf\n", v);
#define DO_THING(x)\
    std::fprintf(stderr, STRINGI(x) " distance is %lf\n", dist::x<>()(a,b,c))
    DO_THING(JaccardSimilarity);
    DO_THING(DiceSimilarity);
    DO_THING(SoftJaccardSimilarity);
    DO_THING(CzekanowskiSim);
    DO_THING(McConnaugheySim);
    DO_THING(MountFordSimilarity);
    DO_THING(HammingDist);
}
