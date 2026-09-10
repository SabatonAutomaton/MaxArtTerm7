#include "ReverseProblem.h"

int main( )
{
   std::vector<double> trueI = { 2,3,4 };
   ReverseProblem problem( 3, trueI);
   for ( double alpha = 1e-20; alpha < 1e-1; alpha *= 10 )
   {
      std::cout << "alpha = " << alpha << std::endl;
      try
      {
         problem.Solve( );
      }
      catch ( const std::runtime_error &e )
      {
         std::cout << "Caught a runtime error: " << e.what( ) << std::endl;
      }
      for ( int i = 0; i < problem.n; i++ )
      {
         std::cout << problem.I[i] << std::endl;
      }
   }


   return 0;
}