#include <stdint.h>
#include <xmmintrin.h>

typedef float floatx4 __attribute__ ((vector_size(16)));

__attribute__ ((noinline))
float averageVectorSize(float *a, uint32_t length) {
  floatx4 sumx4 = {0.0f, 0.0f, 0.0f, 0.0f};
  floatx4 *ax4  = (floatx4 *)a;
  for (uint32_t j = 0, l = length; j < l; j = j + 4) {
    sumx4 = sumx4 + *(ax4++);
  }
  return (sumx4[0] + sumx4[1] + sumx4[2] + sumx4[3])/length;
}

float a[1000];
int main() {
  return (int)averageVectorSize(a, 1000);
}