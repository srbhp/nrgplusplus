#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <stdio.h>
#include <string>
#include <tuple>
#include <vector>
int someFunction(int a, int b) {
  if (b == 0)
    return 1;
  int value = someFunction(a, b / 2);
  if (b % 2 != 0)
    return value * value * a;
  else
    return value * value;
}
int main() {
  int a = 3;
  int b = 5;
  int c = someFunction(a, b);
  printf("%d", c);
  return 0;
}
