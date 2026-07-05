#ifndef FFT_H
#define FFT_H

#include <stdint.h>

#define FFT_SIZE 256
#define FFT_BITS 8


void compute_fft(const float *input_real, float *output_magnitude, uint8_t window_type);

#endif // FFT_H
