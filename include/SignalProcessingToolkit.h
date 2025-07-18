#pragma once

/* *************************************************************************************************

    SignalProcessingToolkit.h

    Centralized header for core signal processing modules used in vibration analysis,
    rotating machinery diagnostics, and general signal analytics.

    Includes modules for:

    - Statistical analysis
    - Spectral estimation (FFT, Welch)
    - Time-frequency transformations (STFT)
    - Angular domain resampling and synchronization
    - Envelope analysis and spectral filtering
    - Interpolation and general utilities

    This header enables easy inclusion of the complete signal processing toolkit
    in a single import.

************************************************************************************************* */

// Utilities & Core Math
#include "utils.h"
#include "statistics.h"

// Spectral Analysis
#include "fft.h"
#include "welch_spectrum.h"     
#include "stft.h"

// Filtering & Modulation
#include "envelope.h"
#include "spectral_filtering.h"

// Angular Domain Processing
#include "angular_processing.h"

// Interpolation
#include "spline.h"
