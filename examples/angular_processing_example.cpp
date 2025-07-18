#include <iostream>
#include <vector>
#include <cmath>
#include "angular_processing.h"  // Main header that includes all modules

int main() {
    using T = double;

    // --------------------------------------------------------------------------------------------
    // 1. Simulation Setup
    // --------------------------------------------------------------------------------------------

    T Fs = 10000.0;            // Sampling frequency [Hz]
    T T_total = 2.0;           // Total signal duration [s]
    int N = static_cast<int>(Fs * T_total);  // Total number of samples
    std::vector<T> t = timeVector<T>(N, Fs); // Time vector

    // --------------------------------------------------------------------------------------------
    // 2. Simulate Shaft Speed Fluctuations
    // --------------------------------------------------------------------------------------------

    T f_rot_mean = 30.0;       // Mean rotational speed [Hz]
    T modAmp = 2.0;            // Amplitude of modulation [Hz]
    T modFreq = 2.0;           // Modulation frequency [Hz]
    int pulsesPerRev = 1;      // One tacho pulse per revolution

    // Generate instantaneous rotation frequency (sinusoidal modulation)
    std::vector<T> f_inst(N);
    for (int i = 0; i < N; ++i)
        f_inst[i] = f_rot_mean + modAmp * std::sin(2 * pi * modFreq * t[i]);

    // Compute instantaneous phase (integral of frequency)
    std::vector<T> phase(N);
    T cumulative = 0.0;
    for (int i = 0; i < N; ++i) {
        cumulative += f_inst[i] / Fs;
        phase[i] = cumulative * 2.0 * pi * pulsesPerRev;
    }

    // --------------------------------------------------------------------------------------------
    // 3. Generate Simulated Tacho Signal
    // --------------------------------------------------------------------------------------------

    // Binary square wave based on instantaneous phase (0 or 1)
    std::vector<T> tacho(N);
    for (int i = 0; i < N; ++i)
        tacho[i] = (std::fmod(phase[i], 2 * pi) < pi) ? 0.0 : 1.0;

    // --------------------------------------------------------------------------------------------
    // 4. Simulate Vibration Signal
    // --------------------------------------------------------------------------------------------

    // Consists of a fixed tone and a speed-dependent 2x-order component
    std::vector<T> signal(N, 0.0);
    for (int i = 0; i < N; ++i) {
        T freq1 = 50.0;               // Constant tone
        T freq2 = f_inst[i] * 2.0;    // Second-order component
        signal[i] = std::sin(2 * pi * freq1 * t[i]) +
            0.5 * std::sin(2 * pi * freq2 * t[i] + 0.5);
    }

    // --------------------------------------------------------------------------------------------
    // 5. Estimate Rotation Speed from Tacho
    // --------------------------------------------------------------------------------------------

    // Automatically suggest a trigger threshold for edge detection
    T trigger = triggerSuggest(tacho);

    // Estimate rotation frequency [Hz] and shaft revolutions
    std::vector<T> speedHz, revs;
    RotationFrequencyFromTacho(tacho, Fs, pulsesPerRev, trigger, speedHz, revs);

    std::cout << "Estimated shaft speed [Hz] at first 5 points:\n";
    for (int i = 0; i < 5 && i < speedHz.size(); ++i)
        std::cout << "  " << speedHz[i] << " Hz\n";

    // --------------------------------------------------------------------------------------------
    // 6. Angular Resampling
    // --------------------------------------------------------------------------------------------

    std::vector<T> signal_resampled, angle_resampled;
    int samplesPerRev = 0, totalRevs = 0;

    AngularResamplingWithTacho(signal, Fs, tacho, pulsesPerRev, trigger,
        signal_resampled, angle_resampled,
        samplesPerRev, totalRevs);

    std::cout << "\nAngular resampling complete: " << totalRevs << " revolutions, "
        << samplesPerRev << " samples per revolution.\n";

    // --------------------------------------------------------------------------------------------
    // 7. Time Synchronous Averaging (TSA)
    // --------------------------------------------------------------------------------------------

    std::vector<T> tsa_signal, tsa_angle;

    TSA(signal_resampled, angle_resampled, samplesPerRev, totalRevs,
        tsa_signal, tsa_angle);

    std::cout << "\nTSA complete. First 5 averaged values:\n";
    for (int i = 0; i < 5 && i < tsa_signal.size(); ++i)
        std::cout << "  angle " << tsa_angle[i] << " rad => " << tsa_signal[i] << "\n";

    // --------------------------------------------------------------------------------------------
    // End
    // --------------------------------------------------------------------------------------------

    std::cout << "\nPress ENTER to exit...";
    std::cin.ignore();

    return 0;
}
