#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

/* *************************************************************************************************

    @file spline.h

    Implements cubic spline interpolation using a natural spline formulation (second derivatives at
    endpoints equal zero). This is useful for estimating intermediate values between known data
    points with smooth transitions.

    - cubicSpline(): Performs 1D interpolation using a tridiagonal solver.

************************************************************************************************* */

// Cubic spline interpolation (natural boundary conditions)
template <typename T>
std::vector<T> cubicSpline(const std::vector<T>& x_sampled, const std::vector<T>& y_sampled, const std::vector<T>& x_new);

/**
 * @brief Natural cubic spline interpolation.
 *
 * @tparam T               Numeric type (float, double, etc.)
 * @param x_sampled        Sampled x-values (must be strictly increasing)
 * @param y_sampled        Sampled y-values (same size as x_sampled)
 * @param x_new            New x-values at which to interpolate
 * @return std::vector<T>  Interpolated y-values at x_new
 *
 * Uses "natural" boundary conditions: second derivative at both ends is zero.
 * Internally solves a tridiagonal system to obtain the second derivatives.
 */
template <typename T>
std::vector<T> cubicSpline(const std::vector<T>& x_sampled, const std::vector<T>& y_sampled, const std::vector<T>& x_new) {

    int N = static_cast<int>(x_sampled.size());
    int N_new = static_cast<int>(x_new.size());
    std::vector<T> y_new(N_new);

    // Tridiagonal system setup
    std::vector<T> diagA(N, 0.0);      // Main diagonal
    std::vector<T> diagInfA(N, 0.0);   // Lower diagonal
    std::vector<T> diagSupA(N, 0.0);   // Upper diagonal
    std::vector<T> v(N, 0.0);          // Right-hand side vector

    diagA[0] = static_cast<T>(1.0);    // Natural boundary: second derivative = 0
    v[0] = static_cast<T>(0.0);

    for (int i = 1; i < N - 1; i++) {
        T h_prev = x_sampled[i] - x_sampled[i - 1];
        T h_next = x_sampled[i + 1] - x_sampled[i];
        T slope_prev = y_sampled[i] - y_sampled[i - 1];
        T slope_next = y_sampled[i + 1] - y_sampled[i];

        diagA[i] = static_cast<T>(2.0 / 3.0) * (h_prev + h_next);
        diagSupA[i] = static_cast<T>(1.0 / 3.0) * h_next;
        diagInfA[i + 1] = diagSupA[i];
        v[i] = slope_next / h_next - slope_prev / h_prev;
    }

    diagA[N - 1] = static_cast<T>(1.0);  // Natural boundary

    // Solve tridiagonal system (Forward substitution)
    for (int j = 1; j < N; j++) {
        T factor = diagInfA[j] / diagA[j - 1];
        diagA[j] -= factor * diagSupA[j - 1];
        v[j] -= factor * v[j - 1];

        if (diagA[j] == static_cast<T>(0.0))
            throw std::runtime_error("cubicSpline: Division by zero in tridiagonal solver.");
    }

    // Back substitution
    std::vector<T> b(N);
    b[N - 1] = v[N - 1] / diagA[N - 1];
    for (int j = N - 2; j >= 0; j--) {
        b[j] = (v[j] - diagSupA[j] * b[j + 1]) / diagA[j];
    }

    // Interpolation loop
    int index = 0;
    for (int s = 0; s < N_new; s++) {
        // Find the interval x_sampled[index] <= x_new[s] < x_sampled[index+1]
        for (int i = index; i < N - 1; i++) {
            if (x_new[s] >= x_sampled[i] && x_new[s] < x_sampled[i + 1]) {
                index = i;
                break;
            }
        }

        T h = x_sampled[index + 1] - x_sampled[index];
        T deltaY = y_sampled[index + 1] - y_sampled[index];
        T di = y_sampled[index];
        T ai = (b[index + 1] - b[index]) / (3.0 * h);
        T ci = deltaY / h - (1.0 / 3.0) * h * (b[index + 1] + 2.0 * b[index]);
        T dx = x_new[s] - x_sampled[index];

        y_new[s] = ai * std::pow(dx, 3) + b[index] * std::pow(dx, 2) + ci * dx + di;
    }

    return y_new;
}
