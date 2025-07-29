# 🔧 Signal Processing toolkit (C++)

**Signal Processing Toolkit** is a modular C++ library designed for advanced signal analysis, with a special focus on rotating machinery signals.  
It offers tools for time, frequency, and angular domain analysis, making it ideal for vibration diagnostics, condition monitoring, and academic projects.

---

## 📌 Key Features

✅ Custom implementations of:
- Fourier Transforms (DFT/FFT and inverse)
- Short-Time Fourier Transform (STFT)
- Spectrum and PSD estimation (Welch method)
- Natural spline interpolation
- Angular resampling using tacho signals or transmission ratios
- Envelope detection and spectral filtering
- Basic statistics and signal operations

✅ Clean, portable and well-documented code  
✅ High-performance FFT using **FFTW3**  
✅ Modular design via `SignalProcessingToolkit.h`  
✅ Practical examples in the `examples/` folder  

---

## 📂 Project Structure

```
SignalProcessingToolkit/
│
├── include/
│ ├── SignalProcessingToolkit.h # Main header (includes all modules)
│ ├── fft.h
│ ├── stft.h
│ ├── welch_spectrum.h
│ ├── envelope.h
│ ├── angular_processing.h
│ ├── spectral_filtering.h
│ ├── statistics.h
│ ├── utils.h
│ └── spline.h
│
├── examples/
│ ├── fft_example.cpp
│ └── angular_processing_example.cpp
│
├── README.md
└── LICENSE (MIT)

```

---

## ⚙️ Requirements

- C++17 or later
- [FFTW3](http://www.fftw.org/) (Fastest Fourier Transform in the West)

---

### ➕ Installing FFTW3

#### On Debian/Ubuntu (Linux):
```bash
sudo apt-get install libfftw3-dev
```

#### On macOS (using Homebrew):
```bash
brew install fftw
```

#### On Windows:
Download precompiled binaries from the official site:
https://fftw.org/download.html

---

### 🔁 How to Use

1. Include the main header:
```cpp
#include "SignalProcessingToolkit.h"
```

2. Compile your code with FFTW3 linked:
```bash
g++ your_program.cpp -lfftw3 -o your_program
```

3. Run the example:
```bash
./fft_example
```

---
## 📘 Module Overview

### ✔️ fft.h
- DFTslow() / iDFTslow() – Educational implementations of Discrete Fourier Transform and its inverse (O(N²) complexity)

- FFTreal2complex() / iFFTreal2complex() – High-performance FFT and IFFT using FFTW3 for real-valued signals

### ✔️ welch_spectrum.h
- spectrumEstimator() – Computes averaged magnitude spectrum via Welch’s method

- psdEstimator() – Estimates Power Spectral Density (PSD) with windowing and overlap

### ✔️ stft.h
- STFT() – Computes Short-Time Fourier Transform for time-frequency analysis using Hanning windows

### ✔️ envelope.h
- envelope() – Extracts the signal envelope using frequency-domain band demodulation

### ✔️ spectral_filtering.h
- bandpassFilterByFFT() – Filters a signal in a given frequency band via FFT manipulation

- extractPeriodicComponent() – Isolates periodic and residual signal components based on rotational frequency

### ✔️ angular_processing.h
- RotationFrequencyFromTacho() – Computes shaft speed from tachometric pulses

- AngularResamplingWithTacho() – Resamples signal according to angular position from tacho signal

- AngularResamplingWithTau() – Resamples signal between shafts using a known transmission ratio

- TSA() – Performs Time Synchronous Averaging across multiple revolutions

### ✔️ statistics.h
- Basic statistical operations on signals: mean(), variance(), standardDeviation(), RMS()

### ✔️ utils.h
- General-purpose utilities: ABS(), alQuadro(), phasePI(), phase2PI(), timeVector(), freqVector()

### ✔️ spline.h
- cubicSpline() – Natural cubic spline interpolation for smooth signal resampling

---

## 💡 Included Examples

| File                    | Description                                           |
|-----------------------------|-------------------------------------------------------|
| `fft_example.cpp`  | Demonstrates DFT/FFT and their inverses on a synthetic sinusoidal signal         |
| `angular_processing_example.cpp`  | Performs angular resampling from a tachometric signal and applies Time Synchronous Averaging (TSA)      |

---

## 🚀 Suggestions
Include SignalProcessingToolkit.h to access all modules easily

Adapt each module to your application: rotating machinery, audio analysis, vibration monitoring, etc.

The library is modular and extensible — feel free to extend it for your needs!

---

## 👤 Author

Developed by **Giuseppe Dipace**.

---

## 📄 License

This project is licensed under the MIT License.  
See `LICENSE.md` for details.

---

## ⭐ Contributing

Pull requests and suggestions are welcome!  
Feel free to open an [issue](https://github.com/gdipace23/signal-processing-toolkit-cpp/issues) or create a PR.

---
