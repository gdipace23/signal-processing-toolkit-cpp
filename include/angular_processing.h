#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "spline.h"  // Required for cubicSpline interpolation
#include "utils.h"   // Contains timeVector()

/* *************************************************************************************************

    angular_processing.h

    Provides tools for angular signal processing and resampling in rotating machinery analysis.

    Includes:
    - Shaft rotation frequency estimation from tacho signals
    - Trigger threshold estimation
    - Angular resampling using tacho signals
    - Angular resampling via transmission ratio (gear ratio)
    - Time Synchronous Averaging (TSA) for periodic component extraction

    These tools are essential for applications such as order tracking, envelope analysis,
    and fault diagnostics in rotating systems.

************************************************************************************************* */

// Estimate rotation frequency from a tachometer signal
template <typename T>
void RotationFrequencyFromTacho(
    const std::vector<T>& tachoSignal,
    const T samplingFreq,
    const int pulsesPerRevolution,
    const T triggerLevel,
    std::vector<T>& rotationFrequency,
    std::vector<T>& shaftRevolutions
);

// Automatically suggest a trigger level from tacho signal
template <typename T>
T triggerSuggest(const std::vector<T>& tachoSignal);

// Perform angular resampling based on tacho signal
template <typename T>
void AngularResamplingWithTacho(
    const std::vector<T>& timeSignal,
    const T samplingFreq,
    const std::vector<T>& tachoSignal,
    const int pulsesPerRevolution,
    const T triggerLevel,
    std::vector<T>& resampledSignal,
    std::vector<T>& angularPositions,
    int& samplesPerRevolution,
    int& totalRevolutions
);

// Perform angular resampling based on gear transmission ratio
template <typename T>
void AngularResamplingWithTau(
    const std::vector<T>& inputSignalResampled,
    const std::vector<T>& inputAngle,
    const int& inputSamplesPerRev,
    const int& inputRevs,
    const T transmissionRatio,
    std::vector<T>& outputSignalResampled,
    std::vector<T>& outputAngle,
    int& outputSamplesPerRev,
    int& outputRevs
);

// Perform time synchronous averaging on an angular-resampled signal
template <typename T>
void TSA(
    const std::vector<T>& signal_resampled,
    const std::vector<T>& angle,
    const int& samplesPerRev,
    const int& numRevolutions,
    std::vector<T>& signal_averaged,
    std::vector<T>& angle_averaged
);


/**
 * @brief Estimate the instantaneous rotation frequency from a tacho signal.
 *
 * @tparam T Numeric type (e.g. float, double)
 * @param tachoSignal         Input tachometer signal
 * @param samplingFreq        Sampling frequency in Hz
 * @param pulsesPerRevolution Number of pulses per shaft revolution
 * @param triggerLevel        Trigger threshold for rising edge detection
 * @param rotationFrequency   Output vector of rotation frequencies (Hz)
 * @param shaftRevolutions    Output vector of revolution counts
 */
template <typename T>
void RotationFrequencyFromTacho(
    const std::vector<T>& tachoSignal,
    const T samplingFreq,
    const int pulsesPerRevolution,
    const T triggerLevel,
    std::vector<T>& rotationFrequency,
    std::vector<T>& shaftRevolutions)
{
    int N = static_cast<int>(tachoSignal.size());
    std::vector<T> time = timeVector(N, samplingFreq);
    std::vector<T> triggerTimes;

    for (int i = 0; i < N - 1; i++) {
        if (tachoSignal[i] <= triggerLevel && tachoSignal[i + 1] >= triggerLevel) {
            T crossing = (triggerLevel - tachoSignal[i + 1]) / (tachoSignal[i] - tachoSignal[i + 1]) * (time[i] - time[i + 1]) + time[i + 1];
            triggerTimes.push_back(crossing);
        }
    }

    int count = static_cast<int>(triggerTimes.size());

    if (count < 2) {
        rotationFrequency.assign(1, 0);
        shaftRevolutions.assign(1, 0);
        return;
    }

    rotationFrequency.resize(count - 1);
    shaftRevolutions.resize(count - 1);
    for (int i = 0; i < count - 1; i++) {
        rotationFrequency[i] = static_cast<T>(1.0) / ((triggerTimes[i + 1] - triggerTimes[i]) * pulsesPerRevolution);
        shaftRevolutions[i] = static_cast<T>(i) / pulsesPerRevolution;
    }
}


/**
 * @brief Suggest a trigger level based on 60% of the tacho signal amplitude range.
 *
 * @tparam T Numeric type
 * @param tachoSignal  Input tachometer signal
 * @return Suggested trigger level
 */
template <typename T>
T triggerSuggest(const std::vector<T>& tachoSignal) {
    auto [minVal, maxVal] = std::minmax_element(tachoSignal.begin(), tachoSignal.end());
    return (*maxVal - *minVal) * static_cast<T>(0.6);
}


/**
 * @brief Perform angular resampling using a tacho signal.
 *
 * @tparam T Numeric type
 * @param timeSignal           Input time-domain signal
 * @param samplingFreq         Sampling frequency
 * @param tachoSignal          Tachometer signal
 * @param pulsesPerRevolution  Number of pulses per revolution
 * @param triggerLevel         Trigger level for edge detection
 * @param resampledSignal      Output signal resampled over angle
 * @param angularPositions     Angular vector (in radians)
 * @param samplesPerRevolution Number of samples per revolution
 * @param totalRevolutions     Total number of revolutions detected
 */
template <typename T>
void AngularResamplingWithTacho(
    const std::vector<T>& timeSignal,
    const T samplingFreq,
    const std::vector<T>& tachoSignal,
    const int pulsesPerRevolution,
    const T triggerLevel,
    std::vector<T>& resampledSignal,
    std::vector<T>& angularPositions,
    int& samplesPerRevolution,
    int& totalRevolutions)
{
    int N = static_cast<int>(timeSignal.size());
    std::vector<T> time = timeVector(N, samplingFreq);
    std::vector<T> triggerTimes;

    for (int i = 0; i < N - 1; i++) {
        if (tachoSignal[i] <= triggerLevel && tachoSignal[i + 1] >= triggerLevel) {
            T crossing = (triggerLevel - tachoSignal[i + 1]) / (tachoSignal[i] - tachoSignal[i + 1]) * (time[i] - time[i + 1]) + time[i + 1];
            triggerTimes.push_back(crossing);
        }
    }

    totalRevolutions = static_cast<int>(triggerTimes.size()) - 2;

    T meanFreq = 0.0;
    for (int i = 0; i < totalRevolutions + 1; i++) {
        meanFreq += static_cast<T>(1.0) / (triggerTimes[i + 1] - triggerTimes[i]);
    }
    meanFreq /= ((totalRevolutions + 1) * pulsesPerRevolution);

    int samplesPerPulse = static_cast<int>(samplingFreq / (meanFreq * pulsesPerRevolution));
    samplesPerRevolution = samplesPerPulse * pulsesPerRevolution;

    std::vector<T> timeTheta(samplesPerPulse * totalRevolutions);
    angularPositions.resize(samplesPerPulse * totalRevolutions);

    T dtheta = static_cast<T>(2.0 * pi / samplesPerRevolution);
    T delta = static_cast<T>(2.0 * pi / pulsesPerRevolution);
    int shift = 0;

    for (int i = 0; i < totalRevolutions; i++) {
        T t1 = triggerTimes[i];
        T t2 = triggerTimes[i + 1];
        T t3 = triggerTimes[i + 2];

        T Db0 = 2 * delta * t1 * t2 * t2 + delta * t1 * t1 * t3 - 2 * delta * t1 * t1 * t2 - delta * t3 * t3 * t1;
        T Db1 = delta * t3 * t3 + 2 * delta * t1 * t1 - delta * t1 * t1 - 2 * delta * t2 * t2;
        T Db2 = 2 * delta * t2 + delta * t1 - delta * t3 - 2 * delta * t1;
        T D = t3 * t3 * t2 + t1 * t2 * t2 + t1 * t1 * t3 - t1 * t1 * t2 - t2 * t2 * t3 - t3 * t3 * t1;

        T b0 = Db0 / D;
        T b1 = Db1 / D;
        T b2 = Db2 / D;

        for (int k = 0; k < samplesPerPulse; k++) {
            T theta = k * dtheta;
            timeTheta[k + shift] = (1 / (2 * b2)) * (std::sqrt(4 * b2 * (theta - b0) + b1 * b1) - b1);
            angularPositions[k + shift] = (k + shift) * dtheta;
        }

        shift += samplesPerPulse;
    }

    resampledSignal = cubicSpline(time, timeSignal, timeTheta);
}


/**
 * @brief Resample a signal in the angular domain from input shaft to output shaft using gear ratio.
 *
 * @tparam T Numeric type
 * @param inputSignalResampled  Input signal already resampled by angle (shaft 1)
 * @param inputAngle            Angle vector of input signal (shaft 1)
 * @param inputSamplesPerRev    Samples per revolution of input
 * @param inputRevs             Total revolutions of input
 * @param transmissionRatio     Transmission ratio (tau = shaft2/shaft1)
 * @param outputSignalResampled Output signal resampled for shaft 2
 * @param outputAngle           Output angular positions (shaft 2)
 * @param outputSamplesPerRev   Samples per revolution of output
 * @param outputRevs            Total revolutions of output
 */
template <typename T>
void AngularResamplingWithTau(
    const std::vector<T>& inputSignalResampled,
    const std::vector<T>& inputAngle,
    const int& inputSamplesPerRev,
    const int& inputRevs,
    const T transmissionRatio,
    std::vector<T>& outputSignalResampled,
    std::vector<T>& outputAngle,
    int& outputSamplesPerRev,
    int& outputRevs)
{
    int N_in = inputSamplesPerRev * inputRevs;
    std::vector<T> scaledAngle(N_in);
    for (int i = 0; i < N_in; i++)
        scaledAngle[i] = inputAngle[i] * transmissionRatio;

    outputSamplesPerRev = static_cast<int>(inputSamplesPerRev / transmissionRatio);
    outputRevs = static_cast<int>(inputRevs * transmissionRatio);
    int N_out = outputSamplesPerRev * outputRevs;

    outputAngle.resize(N_out);
    T dtheta = static_cast<T>(2.0 * pi / outputSamplesPerRev);
    for (int i = 0; i < N_out; i++)
        outputAngle[i] = i * dtheta;

    outputSignalResampled = cubicSpline(scaledAngle, inputSignalResampled, outputAngle);
}


/**
 * @brief Compute the Time Synchronous Average (TSA) of a signal.
 *
 * This function averages multiple revolutions of an angular-resampled signal
 * to extract its periodic content and reduce noise or asynchronous disturbances.
 *
 * @tparam T                    Numeric type (e.g., float, double)
 * @param signal_resampled      Input signal resampled by angle (length = samplesPerRev * numRevolutions)
 * @param angle                 Input angular vector (same length as signal_resampled)
 * @param samplesPerRev         Number of samples per revolution
 * @param numRevolutions        Number of revolutions in the input
 * @param signal_averaged       Output: Time-synchronous averaged signal (length = samplesPerRev)
 * @param angle_averaged        Output: Averaged angular vector (length = samplesPerRev)
 */
template <typename T>
void TSA(
    const std::vector<T>& signal_resampled,
    const std::vector<T>& angle,
    const int& samplesPerRev,
    const int& numRevolutions,
    std::vector<T>& signal_averaged,
    std::vector<T>& angle_averaged
) {
    signal_averaged.assign(samplesPerRev, static_cast<T>(0.0));  // Initialize output
    int shift = 0;

    // Sum each sample across all revolutions
    for (int rev = 0; rev < numRevolutions; ++rev) {
        for (int i = 0; i < samplesPerRev; ++i)
            signal_averaged[i] += signal_resampled[shift + i];
        shift += samplesPerRev;
    }

    // Divide to compute the mean
    for (int i = 0; i < samplesPerRev; ++i)
        signal_averaged[i] /= static_cast<T>(numRevolutions);

    // Reconstruct angle vector
    angle_averaged.resize(samplesPerRev);
    T delta_angle = static_cast<T>(2.0 * pi / samplesPerRev);
    for (int i = 0; i < samplesPerRev; ++i)
        angle_averaged[i] = i * delta_angle;
}