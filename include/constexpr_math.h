#ifndef __CONSTEXPR_MATH_H_
#define __CONSTEXPR_MATH_H_

#include <array>
#include <cstdint>


// Template constexpr function for calculating fabs as constexpr
template<typename Number>
constexpr Number fabs(Number val) {
    if (val < 0) {
        return -val;
    }
    else {
        return val;
    }
}

template<class T, class dcy = std::decay_t<T>>
constexpr inline std::enable_if_t<std::is_floating_point<T>::value, dcy> inverse(T value) {
    return (value == 0) ? 0.0 : 1.0 / value;
}
template <class T>
constexpr inline std::decay_t<T> sign(T value) {
    return value < 0 ? -1 : 1;
}
template <typename T>
constexpr inline std::decay_t<T> abs(T value) {
    return value * sign(value);
}
template<class T>
constexpr inline std::decay_t<T> power(T const& base, std::size_t const& pow) {
    if (pow == 0) {return 1;}
    else if (pow == 1) {return base;}
    else {
        T result = base;
        for (std::size_t i = 1; i < pow; ++i) {
            result *= base;
        }
        return result;
    }
}
constexpr std::intmax_t factorial(std::intmax_t const& n) {
    if (n == 0) {return 1;}
    std::intmax_t result = n;
    for (std::intmax_t i = n - 1; i > 0; --i) {
        result *= i;
    }
    return result;
}
constexpr static std::size_t max_factorial = 10;//replace with calculated version later




template<std::size_t... Is> struct seq {};

template<std::size_t N, std::size_t... Is>
struct gen_seq: gen_seq < N - 1, N, Is... > {};

template<std::size_t... Is>
struct gen_seq<0, Is...> : seq<Is...> {};

template<class T, std::size_t N>
struct sin_coeffs {
    using array_type = std::array<T, N>;
    constexpr static inline T coeff(std::size_t n) {
        return (n % 2 ? 1 : -1) * inverse((T)factorial((2 * n) - 1));
    }
    template<std::size_t...NS>
    constexpr static array_type _coeffs(seq<NS...>) {
        return {{coeff(NS)...}};
    }
    constexpr static array_type coeffs = _coeffs(gen_seq<N> {});
};


template<class T, std::size_t N = max_factorial, class dcy = std::decay_t<T> >
constexpr std::enable_if_t<std::is_floating_point<dcy>::value, dcy> sin(T x) noexcept {
    constexpr const std::array<dcy, N>& coeffs = sin_coeffs<dcy, N>::coeffs;
    const dcy x_2 = x * x;
    dcy pow = x;
    dcy result = 0;
    for (std::size_t i = 0; i < N; ++i) {
        result += coeffs[i] * pow;
        pow *= x_2;
    }
    return result;
}


#endif