#ifndef LDPC_H__
#define LDPC_H__
#include <climits>
#include <cstdlib>

namespace ldpc {

using std::size_t;

template<typename It, typename F>
void fisher_yates_shuffle(It1 i1, It2 i2, F gen=F()) {
    using std::swap;
    size_t dist = std::distance(i1, i2);
    for(;dist > 1; --dist;) {
        swap(i1[dist], il[gen()%dist]);
    }
}
template<typename Con, typename F>
void fisher_yates_shuffle(Con &c, F gen=F()) {
    fisher_yates_shuffle(std::begin(c), std::end(c));
}

// Based on https://github.com/yunwilliamyu/opal code

template<typename T, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
std::vector<T, Alloc> generate_ldpc(int rowlen, int ones_per_row, int height, bool unaltered_include=false) {
    if(rowlen % ones_per_row) throw "up";
    size_t rat = rowlen / ones_per_row;
    static constexpr size_t BPE = sizeof(T) * CHAR_BIT;
    const size_t items_per_row = (rowlen + (BPE - 1)) / BPE;
    size_t nrows = ((height + 1) * rat);
    std::vector<char> base(rowlen * rat),;
    for(size_t i = 0; i < rat; ++i) // For each row
        for(size_t j = i * ones_per_row; j != (i + 1) * ones_per_row; ++j)
            base[i * rowlen + j] = 1;
    auto copy(base.size());
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
    for(size_t i = 0; i < nrows; ++i) {
        fisher_yates_shuffle(copy); {
        for(size_t j = 0; j < rowlen; ++j) {
            ret[offset + items_per_row * i + j / BPE] |= T(base[copy[j]]) << (j % BPE);
        }
    }
    return ret;
}

}

#endif /* LDPC_H__ */