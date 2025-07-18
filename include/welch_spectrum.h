#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include "utils.h"

/* *************************************************************************************************

    welch_estimator.h

    Provides spectral analysis using the Welch method:

    - spectrumEstimator: computes the averaged magnitude spectrum (|FFT|) using overlapping windows
    - psdEstimator: estimates the Power Spectral Density (PSD) using Welch’s method

    These estimators are suitable for analyzing noisy or non-stationary signals with improved stability.

************************************************************************************************* */

// Estimate signal spectrum (magnitude) using Welch method
template <typename T>
std::vector<T> spectrumEstimator(const std::vector<T>& signal, const int Nw, const int Noverlap);

// Estimate Power Spectral Density (PSD) using Welch method
template <typename T>
std::vector<T> psdEstimator(const std::vector<T>& signal, const int Nw, const int Noverlap, const T Fs);


/**
 * @brief Estimate the averaged magnitude spectrum using Welch's method.
 *
 * @tparam T            Numeric type (float, double)
 * @param signal        Input time-domain signal
 * @param Nw            Window length (samples)
 * @param Noverlap      Number of overlapping samples between segments
 * @return vector<T>    Averaged amplitude spectrum (same length as window)
 */
template <typename T>
std::vector<T> spectrumEstimator(const std::vector<T>& signal, const int Nw, const int Noverlap) {
    std::vector<T> w(Nw);
    T Aw = static_cast<T>(0.0);

    // Hanning window and amplitude correction factor
    for (int i = 0; i < Nw; i++) {
        w[i] = 0.5 * (1.0 - std::cos(2.0 * pi * i / (Nw - 1)));
        Aw += w[i];
    }
    Aw = Nw / Aw;

    std::vector<T> spectrum(Nw, 0.0);
    int N = static_cast<int>(signal.size());
    int M = (N - Noverlap) / (Nw - Noverlap);
    int shift = 0;

    for (int m = 0; m < M; ++m) {
        for (int k = 0; k < Nw; ++k) {
            std::complex<T> sum(0.0, 0.0);
            for (int n = 0; n < Nw; ++n) {
                T sample = w[n] * signal[n + shift] * Aw / Nw;
                T angle = -2.0 * pi * k * n / Nw;
                sum += sample * std::polar(static_cast<T>(1.0), angle);
            }
            spectrum[k] += std::pow(std::abs(sum), 2.0);
        }
        shift += (Nw - Noverlap);
    }

    for (int i = 0; i < Nw; ++i)
        spectrum[i] = std::sqrt(spectrum[i] / M);

    return spectrum;
}


/**
 * @brief Estimate the Power Spectral Density (PSD) using Welch's method.
 *
 * @tparam T            Numeric type (float, double)
 * @param signal        Input time-domain signal
 * @param Nw            Window length (samples)
 * @param Noverlap      Number of overlapping samples between segments
 * @param Fs            Sampling frequency [Hz]
 * @return vector<T>    Power Spectral Density estimate (length = Nw)
 */
template <typename T>
std::vector<T> psdEstimator(const std::vector<T>& signal, const int Nw, const int Noverlap, const T Fs) {
    std::vector<T> w(Nw);
    T sumW = 0.0, numW = 0.0;

    // Hanning window and correction factors
    for (int i = 0; i < Nw; ++i) {
        w[i] = 0.5 * (1.0 - std::cos(2 * pi * i / (Nw - 1)));
        sumW += w[i];
        numW += std::pow(w[i], 2.0);
    }

    T Aw = Nw / sumW;                             // Amplitude correction
    T Be = Fs * numW / std::pow(sumW, 2.0);       // Power correction
    T Sp = std::pow(Aw, 2.0) / (std::pow(Nw, 2) * Be);

    std::vector<T> PSD(Nw, 0.0);
    int N = static_cast<int>(signal.size());
    int M = (N - Noverlap) / (Nw - Noverlap);
    int shift = 0;

    for (int m = 0; m < M; ++m) {
        for (int k = 0; k < Nw; ++k) {
            std::complex<T> sum(0.0, 0.0);
            for (int n = 0; n < Nw; ++n) {
                T angle = -2.0 * pi * k * n / Nw;
                sum += w[n] * signal[n + shift] * std::polar(static_cast<T>(1.0), angle);
            }
            PSD[k] += std::pow(std::abs(sum), 2.0) * Sp / M;
        }
        shift += (Nw - Noverlap);
    }

    return PSD;
}
