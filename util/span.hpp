// For license details see LICENSE.

#pragma once

#include <cassert>
#include <cstddef>

template <typename T, std::size_t N> struct static_span {
    typedef T value_type;

    constexpr static_span(const T *ptr) : _data(ptr) {}
    constexpr static_span(const static_span &other) = default;
    constexpr static_span(static_span &&other) = default;

    template <typename Cont>
    constexpr static_span(const Cont &c) : _data(c.data()) {
        static_assert(c.size() == N);
    }

    constexpr const T *cbegin() const { return _data; }
    constexpr const T *cend() const { return _data + N; }
    constexpr const T *begin() const { return _data; }
    constexpr const T *end() const { return _data + N; }

    constexpr const T &operator[](size_t i) const {
        assert(i < N);
        return _data[i];
    }

    constexpr const T *data() const { return _data; }
    constexpr std::size_t size() const { return N; }

  private:
    const T *_data;
};

template <typename T> using span3 = static_span<T, 3>;

typedef span3<double> span3d;

template <typename T> inline span3<T> make_span3(const T *ptr) {
    return span3<T>{ptr};
}

template <typename Cont>
inline span3<typename Cont::value_type> make_span3(const Cont &c) {
    return span3<typename Cont::value_type>{c.data()};
}
