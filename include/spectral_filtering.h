#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include "fft.h"
#include "utils.h"

/* *************************************************************************************************

    spectral_filtering.h

    Provides spectral-domain filtering utilities for time-domain signal preprocessing.

    Includes:
    - Band-pass filtering via DFT with custom frequency bounds (bandpassFilterByFFT)
    - Separation of periodic and residual components based on shaft frequency (extractPeriodicComponent)

    These functions are useful for applications such as envelope analysis, residual extraction,
    and pre-processing for fault detection in rotating machinery.

************************************************************************************************* */


/**
 * @brief Applies a band-pass filter in the frequency domain using DFT.
 *
 * @tparam T         Numeric type (e.g. float, double)
 * @param signal     Input time-domain signal
 * @param Fs         Sampling frequency [Hz]
 * @param Fstart     Lower frequency bound of the filter [Hz]
 * @param Fend       Upper frequency bound of the filter [Hz]
 * @return vector<T> Time-domain signal after filtering
 *
 * Filters the input signal by zeroing all frequency components outside
 * the range [Fstart, Fend] and applying the inverse DFT.
 * Note: Fstart and Fend must satisfy 0 <= Fstart < Fend < Fs/2
 */
template <typename T>
std::vector<T> bandpassFilterByFFT(const std::vector<T>& signal, const T Fs, const T Fstart, const T Fend) {
    int N = static_cast<int>(signal.size());
    std::vector<std::complex<T>> spectrum = DFTslow(signal);  // or FFT if applicable
    std::vector<std::complex<T>> filteredSpectrum(N, std::complex<T>(0.0, 0.0));

    int index_start = static_cast<int>(Fstart * N / Fs);
    int index_end = static_cast<int>(Fend * N / Fs);

    if (index_start == 0) {
        filteredSpectrum[0] = spectrum[0];
    }
    else {
        filteredSpectrum[index_start] = spectrum[index_start];
        filteredSpectrum[N - index_start] = spectrum[N - index_start];
    }

    for (int i = index_start + 1; i <= index_end; ++i) {
        filteredSpectrum[i] = spectrum[i];
        filteredSpectrum[N - i] = spectrum[N - i];
    }

    return iDFTslow(filteredSpectrum);  // or iFFT if applicable
}


/**
 * @brief Decomposes a signal into periodic and residual components based on shaft frequency.
 *
 * @tparam T                  Numeric type (e.g. float, double)
 * @param signal              Input time-domain signal
 * @param signal_periodic     Output periodic signal containing only harmonics of RotationFrequency
 * @param signal_residual     Output residual signal (all non-periodic content)
 * @param Fs                  Sampling frequency [Hz]
 * @param RotationFrequency   Fundamental shaft rotation frequency [Hz]
 * @param epsilonFrequency    Frequency tolerance to identify harmonics [Hz]
 *
 * The function identifies all harmonics of RotationFrequency in the spectrum
 * and isolates them in `signal_periodic`. The remaining content is returned
 * as `signal_residual`.
 */
template <typename T>
void extractPeriodicComponent(
    const std::vector<T>& signal,
    std::vector<T>& signal_periodic,
    std::vector<T>& signal_residual,
    const T Fs,
    const T RotationFrequency,
    const T epsilonFrequency
) {
    std::vector<std::complex<T>> spectrum = DFTslow(signal);  // or FFT
    std::vector<std::complex<T>> periodicSpectrum(spectrum.size(), std::complex<T>(0.0, 0.0));
    std::vector<std::complex<T>> residualSpectrum = spectrum;

    residualSpectrum[0] = std::complex<T>(0.0, 0.0);  // Remove DC if needed

    std::vector<T> freq = freqVector<T>(spectrum.size(), Fs);
    int count = 1;

    for (size_t i = 0; i < spectrum.size(); ++i) {
        if (std::abs(count * RotationFrequency - freq[i]) <= epsilonFrequency) {
            periodicSpectrum[i] = spectrum[i];
            residualSpectrum[i] = std::complex<T>(0.0, 0.0);
        }

        if (freq[i] > count * RotationFrequency + epsilonFrequency) {
            ++count;
        }
    }

    signal_periodic = iDFTslow(periodicSpectrum);  // or iFFT
    signal_residual = iDFTslow(residualSpectrum);  // or iFFT
}
