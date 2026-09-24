#include "ReverseProblem.h"
#include <ctime>

int main()
{
   srand(static_cast<unsigned>(time(nullptr)));

   ReverseProblem problem(3, 1.0);
   problem.lnSigma = std::log(0.1);
   problem.sigma = 0.1;
   problem.eps = 1e-15;
   try
   {
      problem.Solve();
   }
   catch (const std::runtime_error& e)
   {
      std::cout << "Caught a runtime error: " << e.what() << std::endl;
   }
   
   std::cout << "Final sigma = " << problem.sigma << std::endl;

   return 0;
}