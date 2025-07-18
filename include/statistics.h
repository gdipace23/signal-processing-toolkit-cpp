#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

/* *************************************************************************************************

    statistics.h

    A collection of statistical functions for signal processing.
    These functions compute common statistical indicators such as mean, standard deviation,
    RMS, skewness, kurtosis, crest factor, etc., and are useful in condition monitoring
    and vibration analysis.

************************************************************************************************* */
template <typename T> T mean(const std::vector<T>& signal);
template <typename T> T standardDeviation(const std::vector<T>& signal);
template <typename T> T variance(const std::vector<T>& signal);
template <typename T> T rms(const std::vector<T>& signal);
template <typename T> T kurtosis(const std::vector<T>& signal);
template <typename T> T skewness(const std::vector<T>& signal);
template <typename T> T averageAmplitude(const std::vector<T>& signal);
template <typename T> T peakToPeak(const std::vector<T>& signal);
template <typename T> T crestFactor(const std::vector<T>& signal);
template <typename T> T clearanceFactor(const std::vector<T>& signal);

/**
 * @brief Compute the arithmetic mean of a signal vector.
 *
 * @tparam T          Numeric type (e.g., float, double)
 * @param signal      Input signal vector
 * @return T          Mean (average) value of the signal
 */
template <typename T>
T mean(const std::vector<T>& signal) {
    T sum = static_cast<T>(0.0);
    for (const auto& x : signal)
        sum += x;
    return sum / static_cast<T>(signal.size());
}

/**
 * @brief Compute the standard deviation of a signal vector.
 *
 * @tparam T              Numeric type (e.g., float, double)
 * @param signal          Input signal vector
 * @return T              Standard deviation (unbiased estimate)
 */
template <typename T>
T standardDeviation(const std::vector<T>& signal) {
    T mu = mean(signal);
    T accum = static_cast<T>(0.0);
    for (const auto& x : signal)
        accum += std::pow(std::abs(x - mu), 2.0);
    return std::sqrt(accum / static_cast<T>(signal.size() - 1));
}

/**
 * @brief Compute the variance of a signal vector.
 *
 * @tparam T              Numeric type (e.g., float, double)
 * @param signal          Input signal vector
 * @return T              Variance (unbiased estimate)
 */
template <typename T>
T variance(const std::vector<T>& signal) {
    T mu = mean(signal);
    T accum = static_cast<T>(0.0);
    for (const auto& x : signal)
        accum += std::pow(std::abs(x - mu), 2.0);
    return accum / static_cast<T>(signal.size() - 1);
}

/**
 * @brief Compute the Root Mean Square (RMS) of a signal vector.
 *
 * @tparam T          Numeric type (e.g., float, double)
 * @param signal      Input signal vector
 * @return T          RMS value
 */
template <typename T>
T rms(const std::vector<T>& signal) {
    T accum = static_cast<T>(0.0);
    for (const auto& x : signal)
        accum += x * x;
    return std::sqrt(accum / static_cast<T>(signal.size()));
}

/**
 * @brief Compute the kurtosis of the signal.
 *
 * Kurtosis measures the "tailedness" of the signal distribution.
 *
 * @tparam T        Numeric type (e.g., float, double)
 * @param signal    Input signal vector
 * @return T        Kurtosis value (unitless)
 */
template <typename T>
T kurtosis(const std::vector<T>& signal) {
    T mu = mean(signal);
    T num = static_cast<T>(0.0);
    T denom = static_cast<T>(0.0);

    for (const auto& x : signal) {
        T diff = x - mu;
        num += std::pow(diff, 4.0);
        denom += std::pow(diff, 2.0);
    }

    num /= signal.size();
    denom = std::pow(denom / signal.size(), 2.0);

    return num / denom;
}

/**
 * @brief Compute the skewness of the signal.
 *
 * Skewness measures the asymmetry of the signal distribution.
 *
 * @tparam T        Numeric type (e.g., float, double)
 * @param signal    Input signal vector
 * @return T        Skewness value (unitless)
 */
template <typename T>
T skewness(const std::vector<T>& signal) {
    T mu = mean(signal);
    T num = static_cast<T>(0.0);
    T denom = static_cast<T>(0.0);

    for (const auto& x : signal) {
        T diff = x - mu;
        num += std::pow(diff, 3.0);
        denom += std::pow(diff, 2.0);
    }

    num /= signal.size();
    denom = std::pow(denom / signal.size(), 1.5);

    return num / denom;
}

/**
 * @brief Compute the average amplitude of the signal.
 *
 * @tparam T        Numeric type
 * @param signal    Input signal vector
 * @return T        Mean absolute value
 */
template <typename T>
T averageAmplitude(const std::vector<T>& signal) {
    T sum = static_cast<T>(0.0);
    for (const auto& x : signal)
        sum += std::abs(x);
    return sum / static_cast<T>(signal.size());
}

/**
 * @brief Compute the peak-to-peak amplitude of the signal.
 *
 * @tparam T        Numeric type
 * @param signal    Input signal vector
 * @return T        Peak-to-peak value (max - min)
 */
template <typename T>
T peakToPeak(const std::vector<T>& signal) {
    auto [minIt, maxIt] = std::minmax_element(signal.begin(), signal.end());
    return *maxIt - *minIt;
}

/**
 * @brief Compute the crest factor of the signal.
 *
 * Crest Factor = (Peak amplitude) / RMS
 * Here, peak is approximated as half the peak-to-peak.
 *
 * @tparam T        Numeric type
 * @param signal    Input signal vector
 * @return T        Crest factor (unitless)
 */
template <typename T>
T crestFactor(const std::vector<T>& signal) {
    T peak = peakToPeak(signal) / static_cast<T>(2.0);
    T rmsVal = rms(signal);
    return peak / rmsVal;
}

/**
 * @brief Compute the clearance factor of the signal.
 *
 * Clearance Factor = (Peak amplitude) / (mean absolute value)^2
 *
 * @tparam T        Numeric type
 * @param signal    Input signal vector
 * @return T        Clearance factor (unitless)
 */
template <typename T>
T clearanceFactor(const std::vector<T>& signal) {
    T peak = peakToPeak(signal) / static_cast<T>(2.0);
    T avgAbs = averageAmplitude(signal);
    return peak / (avgAbs * avgAbs);
}
