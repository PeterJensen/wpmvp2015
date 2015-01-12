#include <stdint.h>
#include <xmmintrin.h>
__attribute__ ((noinline))
float averageIntrin(float *a, uint32_t length) {
  __m128 sumx4 = _mm_set_ps1(0.0);
  for (uint32_t j = 0, l = length; j < l; j = j + 4) {
    sumx4 = _mm_add_ps(sumx4, _mm_loadu_ps(&(a[j])));
  }
  float mSumx4[4];
  _mm_storeu_ps(mSumx4, sumx4);
  return (mSumx4[0] + mSumx4[1] + mSumx4[2] + mSumx4[3])/length;
}

float a[1000];
int main() {
  return (int)averageIntrin(a, 1000);
}