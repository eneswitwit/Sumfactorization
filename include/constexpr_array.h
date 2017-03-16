#ifndef __ARRAY_CONSTEXPR_H__
#define __ARRAY_CONSTEXPR_H__

template<typename T, size_t N>
class constexpr_array{
private:
    constexpr static size_t size_ = N;
    T data_[N] {}; // T default constructor is essential!
public:
    constexpr size_t size() const { return N;};

    constexpr T& operator[](size_t n){
        return data_[n];
    }

    constexpr const T& operator[](size_t n) const {
        return data_[n];
    }

    using iterator = T*;

    constexpr iterator begin() { return &data_[0];};
    constexpr iterator end() { return &data_[N];};


};


#endif