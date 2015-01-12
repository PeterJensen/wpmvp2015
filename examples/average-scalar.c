#include <stdint.h>

__attribute__ ((noinline))
float averageScalar(float *a, uint32_t length) {
  float  sum = 0.0f;
  for (uint32_t j = 0, l = length; j < l; j = j + 4) {
    sum = sum + *(a++);
  }
  return sum/length;
}

float a[1000];
int main() {
  return (int)averageScalar(a, 1000);
}
