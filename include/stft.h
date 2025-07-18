#pragma once

#include <vector>
#include <cmath>
#include <complex>
#include "utils.h"

/* *************************************************************************************************

    stft.h

    Provides Short-Time Fourier Transform (STFT) functionality for time-frequency analysis of signals.

    Includes:
    - STFT with Hanning window and customizable overlap
    - Returns the time-frequency magnitude map

    Useful for analyzing non-stationary signals and detecting transient or speed-dependent phenomena.

************************************************************************************************* */


/**
 * @brief Compute the Short-Time Fourier Transform (STFT) of a signal using a Hanning window.
 *
 * @tparam T              Numeric type (float, double)
 * @param signal          Input time-domain signal
 * @param Fs              Sampling frequency [Hz]
 * @param Nw              Window size (number of samples per segment)
 * @param Noverlap        Number of overlapping samples between adjacent windows
 * @param stft            Output time-frequency matrix of magnitudes [time][frequency]
 * @param time            Output time vector corresponding to each window center [s]
 * @param frequency       Output frequency vector [Hz]
 */
template <typename T>
void STFT(
    const std::vector<T>& signal,
    const T Fs,
    const int Nw,
    const int Noverlap,
    std::vector<std::vector<T>>& stft,
    std::vector<T>& time,
    std::vector<T>& frequency
) {
    const int N = static_cast<int>(signal.size());

    // Compute Hanning window
    std::vector<T> window(Nw);
    T amplitudeCorrection = static_cast<T>(0.0);
    for (int i = 0; i < Nw; ++i) {
        window[i] = static_cast<T>(0.5) * (1.0 - std::cos(2 * pi * i / (Nw - 1)));
        amplitudeCorrection += window[i];
    }
    amplitudeCorrection = Nw / amplitudeCorrection;

    // Number of time segments
    int M = (N - Noverlap) / (Nw - Noverlap);
    int shift = 0;

    // Resize output containers
    stft.resize(M, std::vector<T>(Nw));
    frequency.resize(Nw);
    time.resize(M);

    for (int m = 0; m < M; ++m) {
        for (int k = 0; k < Nw; ++k) {
            std::complex<T> sum(0.0, 0.0);
            for (int n = 0; n < Nw; ++n) {
                T sample = window[n] * signal[n + shift] * amplitudeCorrection / static_cast<T>(Nw);
                T angle = -2.0 * pi * k * n / static_cast<T>(Nw);
                sum += sample * std::polar(static_cast<T>(1.0), angle);
            }
            stft[m][k] = std::abs(sum);
            frequency[k] = k * Fs / Nw;
        }

        time[m] = m * (Nw - Noverlap) / Fs;
        shift += (Nw - Noverlap);
    }
}
