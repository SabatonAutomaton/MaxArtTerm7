#include "ReverseProblem.h"
#include "Gauss.h"

double ReverseProblem::Distance( xyz p1, xyz p2 )
{
   return sqrt( ( p2.x - p1.x ) * ( p2.x - p1.x ) +
      ( p2.y - p1.y ) * ( p2.y - p1.y ) +
      ( p2.z - p1.z ) * ( p2.z - p1.z ) );
}

void ReverseProblem::ComputeTrueV( )
{
   double piSigma = 0.5 * 1 / ( std::numbers::pi * sigma );
   for ( int i = 0; i < n; i++ )
   {
      trueV[i] = piSigma * trueI[i] * ( 1 / Distance( sourceBPos[i], receiverMPos[i] ) - 1 / Distance( sourceAPos[i], receiverMPos[i] ) -
         ( 1 / Distance( sourceBPos[i], receiverNPos[i] ) - 1 / Distance( sourceAPos[i], receiverNPos[i] ) ) );
   }
}

void ReverseProblem::ComputeCurrentV( )
{
   double piSigma = 0.5 * 1 / ( std::numbers::pi * sigma );
   for ( int i = 0; i < n; i++ )
   {
      V[i] = piSigma * I[i] * ( 1 / Distance( sourceBPos[i], receiverMPos[i] ) - 1 / Distance( sourceAPos[i], receiverMPos[i] ) -
         ( 1 / Distance( sourceBPos[i], receiverNPos[i] ) - 1 / Distance( sourceAPos[i], receiverNPos[i] ) ) );
   }
}

void ReverseProblem::ComputeJacobian( )
{
   double piSigma = 0.5 * 1 / ( std::numbers::pi * sigma );
   for ( int i = 0; i < n; i++ )
      for ( int j = 0; j < n; j++ )
      {
         if ( i == j )
            dVdI[i][j] = piSigma * ( 1 / Distance( sourceBPos[j], receiverMPos[i] ) - 1 / Distance( sourceAPos[j], receiverMPos[i] ) -
               ( 1 / Distance( sourceBPos[j], receiverNPos[i] ) - 1 / Distance( sourceAPos[j], receiverNPos[i] ) ) );
         else
            dVdI[i][j] = 0;
      }
}

void ReverseProblem::WeightsInit( )
{
   for ( int i = 0; i < n; i++ )
      w2[i] = 1 / ( trueV[i] * trueV[i] );
}

void ReverseProblem::BuildMatrix( )
{
   for ( int q = 0; q < n; q++ )
      for ( int s = 0; s < n; s++ )
         for ( int i = 0; i < n; i++ )
            A[q][s] += w2[i] * dVdI[i][q] * dVdI[i][s];
}

void ReverseProblem::BuildRightPart( )
{
   for ( int q = 0; q < n; q++ )
      for ( int i = 0; i < n; i++ )
         b[q] += w2[i] * ( trueV[i] - V[i] ) * dVdI[i][q];
}

void ReverseProblem::Regularisation( )
{
   for ( int i = 0; i < n; i++ )
      A[i][i] += alpha;
}

double ReverseProblem::CalcResidual( )
{
   double res = 0;
   for ( int i = 0; i < n; i++ )
   {
      double diff = V[i] - trueV[i];
      res += w2[i] * diff * diff;
   }
   return res;
}

void ReverseProblem::CalcIndividualResiduals( )
{
   for ( int i = 0; i < n; i++ )
   {
      double diff = V[i] - trueV[i];
      residuals[i] = w2[i] * diff * diff;
   }
}

void ReverseProblem::Solve( )
{
   double residual = 1;
   Gauss gauss;
   ComputeTrueV( );
   WeightsInit( );
   ComputeJacobian( );

   std::cout << "Iter\tI1\tI2\tI3\tTotal_Residual\tPhi1\tPhi2\tPhi3" << std::endl;

   ComputeCurrentV( );
   residual = CalcResidual( );
   CalcIndividualResiduals( );
   std::cout << 0 << "\t" << I[0] << "\t" << I[1] << "\t" << I[2] << "\t" << residual << "\t";
   for ( int i = 0; i < n; i++ )
   {
      std::cout << residuals[i];
      if ( i < n - 1 ) std::cout << "\t";
   }
   std::cout << std::endl;

   for ( int iter = 0; iter < maxIter && residual > eps; iter++ )
   {
      for ( int i = 0; i < n; i++ )
         for ( int j = 0; j < n; j++ )
            A[i][j] = 0;

      for ( int i = 0; i < n; i++ )
         b[i] = 0;

      ComputeCurrentV( );
      BuildMatrix( );
      BuildRightPart( );
      Regularisation( );

      std::vector<double> deltaI = gauss.solve( A, b );

      for ( int i = 0; i < n; i++ )
      {
         I[i] += deltaI[i];
      }

      ComputeCurrentV( );
      residual = CalcResidual( );
      CalcIndividualResiduals( );
      std::cout << iter + 1 << "\t" << I[0] << "\t" << I[1] << "\t" << I[2] << "\t" <<  residual << "\t";
      for ( int i = 0; i < n; i++ )
      {
         std::cout << residuals[i];
         if ( i < n - 1 ) std::cout << "\t";
      }
      std::cout << std::endl;
   }
}

