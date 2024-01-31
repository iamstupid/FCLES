# FCLES
Funny CPP Linear Equation Solver
Inspired from one zhihu question: https://www.zhihu.com/question/641902455

# Usage

```cpp
#include "FCLES.h"
using var = linearExpr::var;
int main(){
  var a, b;
  a + b = 12;
  a - b = 6;
  std::cout << a*b << endl; // 27
  return 0;
}
```
