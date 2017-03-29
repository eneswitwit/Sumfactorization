#ifndef __CONSTEXPR_ARRAY_H__
#define __CONSTEXPR_ARRAY_H__

#include <cstddef>

template<typename Number, std::size_t N>
class constexpr_array{
private:
    constexpr static std::size_t size_ = N;
    Number data_[N] {}; // Number default constructor is essential!
public:
    constexpr std::size_t size() const { return N;};

    constexpr Number& operator[](std::size_t n){
        return data_[n];
    }

    constexpr const Number& operator[](std::size_t n) const {
        return data_[n];
    }

    using iterator = Number*;

    constexpr iterator begin() { return &data_[0];};
    constexpr iterator end() { return &data_[N];};


};


#endif