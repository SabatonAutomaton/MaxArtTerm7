#include "ReverseProblem.h"

int main( )
{
   std::vector<double> trueI = { 4.3,10.1,1.1 };
   ReverseProblem problem( 3, trueI);

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

   return 0;
}