#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "fft.h"  // Includes DFT/FFT utilities and wrappers for FFTW

int main() {
    // --------------------------------------------------------------------------------------------
    // 1. Signal Parameters
    // --------------------------------------------------------------------------------------------

    const double fs = 1000.0;       // Sampling frequency [Hz]
    const double T = 1.0;           // Total duration [s]
    const size_t N = static_cast<size_t>(fs * T); // Number of samples
    const double f1 = 50.0;         // First sinusoidal component [Hz]
    const double f2 = 120.0;        // Second sinusoidal component [Hz]

    // --------------------------------------------------------------------------------------------
    // 2. Generate Test Signal (2 tones + noise)
    // --------------------------------------------------------------------------------------------

    std::vector<double> signal(N);
    for (size_t n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / fs;
        signal[n] = std::sin(2 * pi * f1 * t)                        // Tone 1
            + 0.5 * std::sin(2 * pi * f2 * t)                  // Tone 2 (lower amplitude)
            + 0.1 * ((rand() / static_cast<double>(RAND_MAX)) - 0.5); // Small noise
    }

    // --------------------------------------------------------------------------------------------
    // 3. Compute FFT Using FFTW
    // --------------------------------------------------------------------------------------------

    std::vector<std::complex<double>> spectrum = FFTreal2complex(signal);

    // Display magnitude spectrum for first 20 bins
    std::cout << "Frequency [Hz]   |   Magnitude\n";
    std::cout << "-----------------|------------\n";
    for (size_t k = 0; k < 20; ++k) {
        double freq = static_cast<double>(k) * fs / N;
        std::cout << freq << "               |   " << std::abs(spectrum[k]) << "\n";
    }

    // --------------------------------------------------------------------------------------------
    // 4. Reconstruct Time-Domain Signal Using Inverse FFT
    // --------------------------------------------------------------------------------------------

    std::vector<double> reconstructed = iFFTreal2complex(spectrum);

    // Print first 10 samples for visual comparison
    std::cout << "\nOriginal vs Reconstructed (first 10 samples):\n";
    for (size_t i = 0; i < 10; ++i) {
        std::cout << "original[" << i << "] = " << signal[i]
            << ",  reconstructed[" << i << "] = " << reconstructed[i] << "\n";
    }

    // --------------------------------------------------------------------------------------------
    // End
    // --------------------------------------------------------------------------------------------

    std::cout << "\nPress ENTER to exit...";
    std::cin.ignore();

    return 0;
}
