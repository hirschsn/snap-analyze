#pragma once

#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>

/** Stores a running statistic of max, min, average, sum of a scalar varibale.
 * Used by sana.cpp.
 */
template <typename T>
struct ScalarStatistics {
    T _max = T(0);
    T _min = T(0);
    T _sum = T(0);
    size_t _n = 0;

    struct ScalarStatistics& sample(const T& sample) {
        if (_n == 0) {
            _max = _min = sample;
        } else {
            _max = std::max(_max, sample);
            _min = std::min(_min, sample);
        }
        _sum += sample;
        _n++;
        return *this;
    }

    size_t nsamples() const { return _n; }
    T max() const { return _max; }
    T min() const { return _min; }
    T sum() const { return _sum; }
    T avg() const {
        return nsamples() == static_cast<size_t>(0)
                ? T(0)
                : sum() / nsamples();
    }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const ScalarStatistics<T>& ss)
{
    return os << ss.min() << " " << ss.avg() << " " << ss.max();
}


template <typename T>
struct Histogram {
    std::vector<T> data;

    Histogram& sample(const T& sample) {
        data.push_back(sample);
        return *this;
    }

    size_t nsamples() const { return data.size(); }
    T max() const {
        if (nsamples() == 0)
            return T(0);
        else
            return *std::max_element(std::begin(data), std::end(data));
    }
    T min() const {
        if (nsamples() == 0)
            return T(0);
        else
            return *std::min_element(std::begin(data), std::end(data));
    }
    T sum() const { return std::accumulate(std::begin(data), std::end(data), T(0)); }
    T avg() const {
        return nsamples() == static_cast<size_t>(0)
                ? T(0)
                : sum() / nsamples();
    }

    //std::vector<size_t> hist(size_t nbins) {
    //    T incr = 1.01 * max() / nbins;
    //    std::vector<size_t> h(nbins, 0);

    //    for (const auto& d: data)
    //        h[d / incr]++;
    //    return h;
    //}

    const auto begin() const { return data.begin(); }
    const auto end() const { return data.end(); }

    //struct boxplot {
    //    T median;
    //    T q1, q3;
    //    T min, max;
    //};
    //boxpot box() {
    //    std::sort(std::begin(data), std::end(data));

    //    boxplot bp;
    //    bp.median = data[nsamples() / 2];
    //    bp.q1 = data[nsamples() / 4];
    //    bp.q3 = data[3 * nsamples() / 4];
    //    bp.min = min();
    //    bp.max = max();

    //    return bp;
    //}
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Histogram<T>& ss)
{
    return os << ss.min() << " " << ss.avg() << " " << ss.max();
}
