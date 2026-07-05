#include "fft.h"

typedef struct {
    float real;
    float imag;
} complex_t;


// taylor series trig approx
static float local_sin(float x) {
    // Normalize angle to (-PI, PI]
    while (x > 3.14159265f)  x -= 6.2831853f;
    while (x <= -3.14159265f) x += 6.2831853f;
    
    // Taylor series expansion for sin(x) = x - x^3/3! + x^5/5! - x^7/7!
    float x2 = x * x;
    float term = x;
    float sum = x;
    
    term = (term * x2) / 6.0f;    // x^3 / 3!
    sum -= term;
    term = (term * x2) / 20.0f;   // x^5 / 5!
    sum += term;
    term = (term * x2) / 42.0f;   // x^7 / 7!
    sum -= term;
    
    return sum;
}

static float local_cos(float x) {
    // cos(x) = sin(x + PI/2)
    return local_sin(x + 1.57079632f);
}

// Newton-Raphson sqrt. appx.
static float local_sqrt(float x) {
    if (x <= 0.0f) return 0.0f;
    
    // Initial rough guess
    float guess = x;
    if (x > 1.0f) guess = x * 0.5f;
    
    // Run 4 quick iterations for clean float accuracy
    for (int i = 0; i < 4; i++) {
        guess = 0.5f * (guess + x / guess);
    }
    return guess;
}

// bit-reversal
static uint8_t fft_bit_reverse(uint8_t index) {
    uint8_t rev = 0;
    for (uint8_t i = 0; i < FFT_BITS; i++) {
        if (index & (1 << i)) {
            rev |= (1u << (FFT_BITS - 1u - i));
        }
    }
    return rev;
}

// Coley-Tukey FFT
void compute_fft(const float *input_real, float *output_magnitude, uint8_t window_type) {
    static complex_t data[FFT_SIZE];
    
    // Bit-reversal sorting step combined with Windowing Function
    for (int i = 0; i < FFT_SIZE; i++) {
        uint8_t rev_idx = fft_bit_reverse((uint8_t)i);
        
        // Default multiplier is 1.0f (Rectangle Window)
        float w = 1.0f; 
        
        // Calculate common angle for Hann, Hamming, and Blackman
        float angle = (2.0f * 3.14159265f * (float)i) / (float)(FFT_SIZE - 1);

        if (window_type == 0) { // HANN
            w = 0.5f * (1.0f - local_cos(angle));
        } 
        else if (window_type == 1) { // HAMMING
            w = 0.54f - 0.46f * local_cos(angle);
        } 
        else if (window_type == 2) { // BLACKMAN
            w = 0.42f - 0.5f * local_cos(angle) + 0.08f * local_cos(2.0f * angle);
        }
        // window_type == 3 is RECTANGLE, leaves w = 1.0f

        // Apply window shaping curve to incoming real voltage trace
        data[rev_idx].real = input_real[i] * w;
        data[rev_idx].imag = 0.0f;
    }

    // Cooley-Tukey Radix-2 processing loop
    for (int size = 2; size <= FFT_SIZE; size <<= 1) {
        int half_size = size >> 1;
        float tab_step = -2.0f * 3.14159265f / (float)size;

        for (int i = 0; i < FFT_SIZE; i += size) {
            for (int j = 0; j < half_size; j++) {
                float angle = (float)j * tab_step;
                
                float twiddle_real = local_cos(angle);
                float twiddle_imag = local_sin(angle);

                int t_idx = i + j + half_size;
                int u_idx = i + j;

                float t_real = data[t_idx].real * twiddle_real - data[t_idx].imag * twiddle_imag;
                float t_imag = data[t_idx].real * twiddle_imag + data[t_idx].imag * twiddle_real;

                data[t_idx].real = data[u_idx].real - t_real;
                data[t_idx].imag = data[u_idx].imag - t_imag;
                data[u_idx].real += t_real;
                data[u_idx].imag += t_imag;
            }
        }
    }

    // Absolute complex magnitude calculation pass
    for (int i = 0; i < FFT_SIZE / 2; i++) {
        float r = data[i].real;
        float im = data[i].imag;
        output_magnitude[i] = local_sqrt(r * r + im * im);
    }
}

// Manually satisfy the bare-metal __aeabi_memclr4 symbol requirement
__attribute__((used)) 
void __aeabi_memclr4(void *dest, unsigned int bytes) {
    uint32_t *d = (uint32_t *)dest;
    unsigned int words = bytes >> 2; // Divide by 4 to get word count
    
    while (words--) {
        *d++ = 0;
    }
}
