#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <cstring>     // for memcpy
#include <stdexcept>   // for exceptions
#include "utils.h"

/* *************************************************************************************************

    fft.h

    Contains functions for spectral analysis, including:

    - Direct Discrete Fourier Transform (DFT) and its inverse (iDFT)
    - Fast Fourier Transform (FFT) and its inverse (IFFT) using FFTW3

    The direct implementations (DFTslow/iDFTslow) are provided for educational and testing purposes.
    FFTreal2complex and iFFTreal2complex leverage FFTW for high performance on real signals.

************************************************************************************************* */

// Compute DFT using a direct implementation (slow, generic)
template <typename T> std::vector<std::complex<T>> DFTslow(const std::vector<T>& signal);

// Compute inverse DFT using a direct implementation (slow, generic)
template <typename T> std::vector<T> iDFTslow(const std::vector<std::complex<T>>& spectrum);

// Compute FFT using FFTW3 (real-to-complex) for double
std::vector<std::complex<double>> FFTreal2complex(const std::vector<double>& signal);

// Compute IFFT using FFTW3 (complex-to-real) for double
std::vector<double> iFFTreal2complex(const std::vector<std::complex<double>>& spectrum);


/**
 * @brief Compute the Discrete Fourier Transform (DFT) of a real signal.
 *
 * @tparam T              Numeric type (float, double, etc.)
 * @param signal          Input signal vector
 * @return vector<complex<T>> Complex frequency-domain spectrum
 * 
 * This function uses the mathematical definition of DFT and has O(N^2) complexity.
 */
template <typename T>
std::vector<std::complex<T>> DFTslow(const std::vector<T>& signal) {
    size_t N = signal.size();
    std::vector<std::complex<T>> spectrum(N);

    // Loop over output frequencies k
    for (size_t k = 0; k < N; ++k) {
        std::complex<T> sum(0.0, 0.0);
        for (size_t n = 0; n < N; ++n) {
            T angle = -2.0 * pi * static_cast<T>(k * n) / static_cast<T>(N);
            sum += signal[n] * std::polar(static_cast<T>(1.0), angle);
        }
        spectrum[k] = sum / static_cast<T>(N); // Normalize
    }

    return spectrum;
}


/**
 * @brief Compute the inverse Discrete Fourier Transform (iDFT).
 *
 * @tparam T              Numeric type (float, double, etc.)
 * @param spectrum        Complex frequency-domain input
 * @return vector<T>      Reconstructed time-domain signal
 *
 * This function uses the mathematical definition of iDFT and has O(N^2) complexity.
 */
template <typename T>
std::vector<T> iDFTslow(const std::vector<std::complex<T>>& spectrum) {
    size_t N = spectrum.size();
    std::vector<T> signal(N);

    // Loop over output time indices n
    for (size_t n = 0; n < N; ++n) {
        std::complex<T> sum(0.0, 0.0);
        for (size_t k = 0; k < N; ++k) {
            T angle = 2.0 * pi * static_cast<T>(k * n) / static_cast<T>(N);
            sum += spectrum[k] * std::polar(static_cast<T>(1.0), angle);
        }
        signal[n] = std::real(sum); 
    }

    return signal;
}


/**
 * @brief Compute the FFT using FFTW (real to complex).
 *
 * @param signal      Input signal (real)
 * @return vector<complex<double>> Half-complex spectrum (size = N/2 + 1)
 *
 * FFTreal2complex is significantly faster than DFTslow. The output contains N/2 + 1 values.
 */
inline std::vector<std::complex<double>> FFTreal2complex(const std::vector<double>& signal) {
    int N = static_cast<int>(signal.size());
    std::vector<std::complex<double>> spectrum(N / 2 + 1);

    // Allocate FFTW buffers
    double* in = (double*)fftw_malloc(sizeof(double) * N);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    
    // Copy data and execute plan
    std::memcpy(in, signal.data(), sizeof(double) * N);
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Copy result
    std::memcpy(spectrum.data(), out, sizeof(std::complex<double>) * (N / 2 + 1));

    // Normalize
    for (auto& c : spectrum)
        c /= static_cast<double>(N);

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return spectrum;
}


/**
 * @brief Compute the inverse FFT using FFTW (complex-to-real).
 *
 * @param spectrum      Input complex spectrum (size = N/2 + 1)
 * @return vector<double> Real reconstructed time-domain signal (size = N)
 *
 * iFFTreal2complex is significantly faster than iDFTslow.
 */
inline std::vector<double> iFFTreal2complex(const std::vector<std::complex<double>>& spectrum) {
    int N = static_cast<int>((spectrum.size() - 1) * 2);
    std::vector<double> signal(N);

    // Allocate FFTW buffers
    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    double* out = (double*)fftw_malloc(sizeof(double) * N);

    // Copy data and execute plan
    std::memcpy(in, spectrum.data(), sizeof(std::complex<double>) * (N / 2 + 1));
    fftw_plan plan = fftw_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Copy result
    std::memcpy(signal.data(), out, sizeof(double) * N);

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return signal;
}
