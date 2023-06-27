
#include <semt/Semt.h>
#include <semt/Shortcuts.h>
#include <iostream>

using namespace std;
using namespace SEMT;

int main()
{
 DVAR(x,0);
 auto f = pow(x, INT(3)) + sin(x);
 cout << f << endl;

 // We have one variable and want four derivatives.
 DifferentiableVectorExpr<1, 4>DVE;
 // Insert the expression, the actual differentiation has taken place
 // long before the whole program is even executed ^^
 // But this line issues the compiler to instantiate the derivatives.
 DVE.push_back(f);

 // Print derivatives and evaluate.
 vector<double>x0(1, 3.1415);
 for (int i = 0; i < 5; ++i)
  cout << DVE.get_derivative(i) << " == "
          << (DVE.get_derivative(i)(x0))
          << " for x0 = " << (x0[0]) << endl;
 
 return 0;
}
