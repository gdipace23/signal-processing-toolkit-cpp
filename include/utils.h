#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <cstdlib> // for std::abs()

constexpr double pi = 3.14159265358979323846;

/* *************************************************************************************************

    utils.h

    Utility functions for signal processing:

    - Vector absolute value (ABS) for real and complex inputs
    - Element-wise squaring (alQuadro)
    - Phase angle extraction (in [−π, π] or [0, 2π])
    - Time and frequency axis generators

************************************************************************************************* */

// Compute absolute value for real-valued vector
template <typename T> std::vector<T> ABS(const std::vector<T>& signal);

// Compute magnitude of complex-valued vector
template <typename T> std::vector<T> ABS(const std::vector<std::complex<T>>& signal);

// Square each element of a real vector
template <typename T> std::vector<T> alQuadro(const std::vector<T>& signal);

// Compute phase angle in range [-π, π]
template <typename T> std::vector<T> phasePI(const std::vector<std::complex<T>>& signal);

// Compute phase angle in range [0, 2π]
template <typename T> std::vector<T> phase2PI(const std::vector<std::complex<T>>& signal);

// Generate a time vector given sampling frequency and number of samples
template <typename T> std::vector<T> timeVector(const size_t N, const T Fs);

// Generate a frequency vector (linear) given Fs and N
template <typename T> std::vector<T> freqVector(const size_t N, const T Fs);


/**
 * @brief Compute the absolute value of each element in a real-valued vector.
 *
 * @tparam T       Numeric type (float, double, etc.)
 * @param signal   Input signal
 * @return         Output vector with absolute values
 */
template <typename T>
std::vector<T> ABS(const std::vector<T>& signal) {
    std::vector<T> absVector(signal.size());
    for (size_t i = 0; i < signal.size(); ++i)
        absVector[i] = std::abs(signal[i]);
    return absVector;
}


/**
 * @brief Compute the magnitude of each element in a complex-valued vector.
 *
 * @tparam T       Numeric type (float, double, etc.)
 * @param signal   Input complex-valued vector
 * @return         Output vector of magnitudes
 */
template <typename T>
std::vector<T> ABS(const std::vector<std::complex<T>>& signal) {
    std::vector<T> absVector(signal.size());
    for (size_t i = 0; i < signal.size(); ++i)
        absVector[i] = std::abs(signal[i]);
    return absVector;
}


/**
 * @brief Square each element in a real-valued vector.
 *
 * @tparam T       Numeric type (float, double, etc.)
 * @param signal   Input vector
 * @return         Output vector with each element squared
 */
template <typename T>
std::vector<T> alQuadro(const std::vector<T>& signal) {
    std::vector<T> squared(signal.size());
    for (size_t i = 0; i < signal.size(); ++i)
        squared[i] = std::pow(signal[i], 2.0);
    return squared;
}


/**
 * @brief Compute phase angles of a complex vector in range [−π, π].
 *
 * @tparam T       Numeric type
 * @param signal   Complex input vector
 * @return         Vector of phase angles in radians
 */
template <typename T>
std::vector<T> phasePI(const std::vector<std::complex<T>>& signal) {
    std::vector<T> phase(signal.size());
    for (size_t i = 0; i < signal.size(); ++i)
        phase[i] = std::arg(signal[i]);
    return phase;
}


/**
 * @brief Compute phase angles of a complex vector in range [0, 2π].
 *
 * @tparam T       Numeric type
 * @param signal   Complex input vector
 * @return         Vector of phase angles in radians, wrapped to [0, 2π]
 */
template <typename T>
std::vector<T> phase2PI(const std::vector<std::complex<T>>& signal) {
    std::vector<T> phase(signal.size());
    for (size_t i = 0; i < signal.size(); ++i) {
        T angle = std::arg(signal[i]);
        phase[i] = (angle < static_cast<T>(0.0)) ? angle + 2.0 * pi : angle;
    }
    return phase;
}


/**
 * @brief Generate a time vector of size N with sampling frequency Fs.
 *
 * @tparam T       Numeric type
 * @param N        Number of samples
 * @param Fs       Sampling frequency
 * @return         Vector of time values
 */
template <typename T>
std::vector<T> timeVector(const size_t N, const T Fs) {
    std::vector<T> t(N);
    for (size_t i = 0; i < N; ++i)
        t[i] = static_cast<T>(i) / Fs;
    return t;
}


/**
 * @brief Generate a frequency vector of size N with sampling frequency Fs.
 *
 * @tparam T       Numeric type
 * @param N        Number of frequency bins
 * @param Fs       Sampling frequency
 * @return         Vector of frequency values (0 to Fs)
 */
template <typename T>
std::vector<T> freqVector(const size_t N, const T Fs) {
    std::vector<T> f(N);
    for (size_t i = 0; i < N; ++i)
        f[i] = static_cast<T>(i) * Fs / static_cast<T>(N);
    return f;
}
