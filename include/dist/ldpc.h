#ifndef LDPC_H__
#define LDPC_H__
#include <cassert>
#include <set>
#include <climits>
#include <random>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <numeric>

namespace ldpc {

using std::size_t;


template<typename It1, typename It2, typename F>
void fisher_yates_shuffle(It1 i1, It2 i2, F gen=F());
template<typename Con, typename F>
void fisher_yates_shuffle(Con &c, F gen=F());

template<typename It1, typename It2, typename F>
void fisher_yates_shuffle(It1 i1, It2 i2, F gen) {
    using std::swap;
    size_t dist = std::distance(i1, i2) - 1;
    for(;dist > 1; --dist) {
        swap(i1[dist], i1[gen()%dist]);
    }
}


template<typename Con, typename F>
void fisher_yates_shuffle(Con &c, F gen) {
    fisher_yates_shuffle<decltype(std::begin(c)), decltype(std::end(c)),F>(std::begin(c), std::end(c), gen);
}

template<typename T, typename Alloc, typename=std::enable_if_t<std::is_arithmetic<T>::value>>
void dump_binary_matrix(const std::vector<T, Alloc> &c, int rowlen, std::FILE *fp=stdout) {
    for(size_t i = 0; i < c.size() / rowlen; ++i) {
        auto p = &c[i * rowlen];
        static constexpr size_t BPE = sizeof(T) * CHAR_BIT;
        for(size_t j = 0; j < (rowlen + (BPE - 1)) / BPE; ++j) {
            auto v = p[j];
            for(size_t k = 0; k < BPE; )
                std::fputc('0' + ((v>>k++)&1), fp);
        }
        std::fputc('\n', fp);
    }
}

template<typename T>void show_vec(const T &x) {for(const auto i: x) std::fprintf(stderr, "%zu\n", size_t(i)); std::fputc('\n', stderr);}

// Based on https://github.com/yunwilliamyu/opal code

template<typename T, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
std::vector<T, Alloc> generate_ldpc(int rowlen, int ones_per_row, int height, bool unaltered_include=false) {
    if(rowlen % ones_per_row) throw "up";
    size_t rat = rowlen / ones_per_row;
    static constexpr size_t BPE = sizeof(T) * CHAR_BIT;
    const size_t items_per_row = (rowlen + (BPE - 1)) / BPE;
    size_t nrows = ((height + 1) * rat);
    std::vector<char> base(rowlen * rat);
    for(size_t i = 0; i < rat; ++i) // For each row
        for(size_t j = i * ones_per_row; j != (i + 1) * ones_per_row; ++j)
            base[i * rowlen + j] = 1;
    std::fprintf(stderr, "BPE: %zu. items per row: %zu\n", BPE, items_per_row);
    std::vector<int> copy(base.size());
    std::iota(copy.begin(), copy.end(), 0u);
    size_t offset = unaltered_include ? items_per_row * rat: size_t(0);
    std::vector<T, Alloc> ret((unaltered_include ? nrows + rat: nrows)  * items_per_row);
    if(unaltered_include) {
        for(size_t i = 0; i < rat; ++i) {
            for(size_t j = 0; j < rowlen; ++j) {
                size_t bit_index = i * rowlen + j;
                ret[bit_index / BPE] |= base[i * rowlen + j] << (bit_index % BPE);
            }
        }
        
    }
#if !NDEBUG
    std::set<int> lset(copy.begin(), copy.end());
#endif
    std::mt19937_64 mt;
    auto retptr = ret.data();
    if(unaltered_include) retptr += items_per_row * rat;
    for(size_t i = 0; i < nrows; ++i) {
        assert(std::set<int>(copy.begin(), copy.end()) == lset);
        fisher_yates_shuffle(copy, mt);
        for(size_t j = 0; j < rowlen; ++j) {
            size_t idx = items_per_row * i + j / BPE;
            assert(j < copy.size());
            auto bidx = copy[j];
            assert(idx < ret.size());
            assert(bidx < base.size());
            retptr[idx] |= T(base[bidx]) << (j % BPE);
        }
    }
    return ret;
}

}

#endif /* LDPC_H__ */
