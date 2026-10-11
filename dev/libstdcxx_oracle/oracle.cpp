// Oracle for the Rust port of libstdc++ random machinery. Compile with the GCC under test.
#include <random>
#include <cstdio>
#include <vector>
#include <algorithm>

int main() {
  std::printf("# compiler %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
  // mt19937 first outputs
  for (unsigned s : {5489u, 1u, 123u, 124u, 0u, 4294967295u}) {
    std::mt19937 mt(s);
    std::printf("mt %u", s);
    for (int i = 0; i < 5; ++i) std::printf(" %lu", (unsigned long)mt());
    std::printf("\n");
  }
  {
    std::mt19937 mt(5489u);
    mt.discard(9999);
    std::printf("mt10000 5489 %lu\n", (unsigned long)mt());
  }
  // seed narrowing (unsigned long long seed)
  {
    unsigned long long big = 4294967296ULL + 123ULL;
    std::mt19937 mt(big);
    std::printf("mtbig %lu\n", (unsigned long)mt());
  }
  // generate_canonical
  {
    std::mt19937 mt(123u);
    std::printf("canon 123");
    for (int i = 0; i < 4; ++i) std::printf(" %.17g", std::generate_canonical<double, 53>(mt));
    std::printf("\n");
  }
  // normal: fresh distribution each call (FAUST's cppRnorm)
  {
    std::mt19937 mt(123u);
    std::printf("normfresh 123");
    for (int i = 0; i < 4; ++i) { std::normal_distribution<> d(0.5, 2.0); std::printf(" %.17g", d(mt)); }
    std::printf("\n");
  }
  // normal: reused distribution (saved second deviate)
  {
    std::mt19937 mt(123u);
    std::normal_distribution<> d(0.0, 1.0);
    std::printf("normreuse 123");
    for (int i = 0; i < 4; ++i) std::printf(" %.17g", d(mt));
    std::printf("\n");
  }
  // uniform_real
  {
    std::mt19937 mt(124u);
    std::uniform_real_distribution<double> d(-1.5, 2.25);
    std::printf("unif 124");
    for (int i = 0; i < 4; ++i) std::printf(" %.17g", d(mt));
    std::printf("\n");
  }
  // bernoulli
  {
    std::mt19937 mt(7u);
    std::bernoulli_distribution d(0.3);
    std::printf("bern 7");
    for (int i = 0; i < 16; ++i) std::printf(" %d", (int)d(mt));
    std::printf("\n");
  }
  // uniform_int<int>
  for (int hi : {0, 1, 6, 99, 999, 65535, 1000000, 2147483647}) {
    std::mt19937 mt(42u);
    std::uniform_int_distribution<> d(0, hi);
    std::printf("uint 42 0 %d", hi);
    for (int i = 0; i < 8; ++i) std::printf(" %d", d(mt));
    std::printf("\n");
  }
  {
    std::mt19937 mt(42u);
    std::uniform_int_distribution<> d(-5, 5);
    std::printf("uint 42 -5 5");
    for (int i = 0; i < 8; ++i) std::printf(" %d", d(mt));
    std::printf("\n");
  }
  // shuffle of vector<long>
  for (int n : {1, 2, 3, 7, 8, 100, 1001}) {
    std::mt19937 mt(99u);
    std::vector<long> v(n);
    for (int i = 0; i < n; ++i) v[i] = i;
    std::shuffle(v.begin(), v.end(), mt);
    std::printf("shuf 99 %d", n);
    int m = std::min(n, 12);
    for (int i = 0; i < m; ++i) std::printf(" %ld", v[i]);
    std::printf(" next %lu\n", (unsigned long)mt());
  }
  // shuffle large (above the two-per-draw threshold: n*n > 2^32-1 -> n >= 65536)
  {
    std::mt19937 mt(99u);
    std::vector<long> v(70000);
    for (int i = 0; i < 70000; ++i) v[i] = i;
    std::shuffle(v.begin(), v.end(), mt);
    std::printf("shuf 99 70000");
    for (int i = 0; i < 12; ++i) std::printf(" %ld", v[i]);
    std::printf(" next %lu\n", (unsigned long)mt());
  }
  return 0;
}
