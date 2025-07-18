#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include "utils.h"

/* *************************************************************************************************

    envelope.h

    Envelope analysis via frequency band demodulation (manual DFT approach).
    Suitable for vibration signals acquired from accelerometers or microphones,
    where fault signatures are hidden in amplitude modulation.

************************************************************************************************* */

/**
 * @brief Compute the envelope of a signal via frequency band demodulation.
 *
 * The envelope is extracted by isolating a frequency band using DFT, then performing
 * an inverse transform and computing the modulus.
 *
 * @tparam T           Numeric type (e.g., float, double)
 * @param signal       Input time-domain signal
 * @param Fs           Sampling frequency [Hz]
 * @param Fstart       Lower bound of band-pass filter [Hz]
 * @param Fend         Upper bound of band-pass filter [Hz]
 * @return vector<T>   Envelope signal (real-valued, same length as input)
 */
template <typename T>
std::vector<T> envelope(const std::vector<T>& signal, T Fs, T Fstart, T Fend)
{
    // Number of samples
    int N = static_cast<int>(signal.size());

    // Define band-pass limits in frequency bin indices
    int index_start = static_cast<int>(Fstart * N / Fs);
    int index_end = static_cast<int>(Fend * N / Fs) + 1;

    // Prepare filtered spectrum (zero outside the band)
    std::vector<std::complex<T>> spectrum_band(index_end - index_start, { 0.0, 0.0 });

    // DFT: compute only for bins within the band
    for (int k = index_start; k < index_end; ++k) {
        std::complex<T> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            T angle = -2.0 * pi * k * n / N;
            sum += signal[n] * std::polar(static_cast<T>(1.0), angle);
        }
        spectrum_band[k - index_start] = sum / static_cast<T>(N);
    }

    // iDFT: reconstruct time-domain signal from band-limited spectrum
    std::vector<T> envelope_signal(N);
    for (int n = 0; n < N; ++n) {
        std::complex<T> sum(0.0, 0.0);
        for (int k = 0; k < static_cast<int>(spectrum_band.size()); ++k) {
            T angle = 2.0 * pi * (k + index_start) * n / N;
            sum += spectrum_band[k] * std::polar(static_cast<T>(1.0), angle);
        }
        envelope_signal[n] = static_cast<T>(2.0) * std::abs(sum);  // real envelope
    }

    return envelope_signal;
}
