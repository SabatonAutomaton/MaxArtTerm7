#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <numbers>
class ReverseProblem
{
public:
   ReverseProblem( int _n, double _trueSigma)
   {
      n = _n;
      trueSigma = _trueSigma;
      lnSigma = std::log(0.01);
      sigma = 0.01;
      A.resize( 1, std::vector<double>( 1, 0 ) );
      b.resize( 1 );
      dVdSigma.resize( n, 0.0 );
      trueV.resize( n );
      V.resize( n );
      w2.resize( n );
      residuals.resize( n );
      receiverMPos = { xyz( 200,0,0 ), xyz( 500,0,0 ), xyz( 1000,0,0 ) };
      receiverNPos = { xyz( 300,0,0 ), xyz( 600,0,0 ), xyz( 1100,0,0 ) };
      sourceAPos = { xyz( 0,-500,0 ), xyz( 0,0,0 ), xyz( 0,500,0 ) };
      sourceBPos = { xyz( 100,-500,0 ), xyz( 100,0,0 ), xyz( 100,500,0 ) };
   }

   struct xyz {

      xyz( ) : x( 0 ), y( 0 ), z( 0 ) {}

      xyz( double _x, double _y, double _z ) : x( _x ), y( _y ), z( _z ) {}

      double x = 0;
      double y = 0;
      double z = 0;
   };

   int n = 0;
   std::vector<std::vector<double>> A;
   std::vector<double> b;
   std::vector<double> trueV;
   std::vector<double> V;
   std::vector<xyz> receiverMPos;
   std::vector<xyz> receiverNPos;
   std::vector<xyz> sourceAPos;
   std::vector<xyz> sourceBPos;
   double lnSigma;
   double sigma;
   double trueSigma;
   std::vector<double> dVdSigma;
   std::vector<double> w2;
   std::vector<double> residuals;
   double alpha = 0.0;
   double eps = 1e-10;
   int maxIter = 10;

   double Distance( xyz p1, xyz p2 );
   void ComputeTrueV( );
   void ComputeCurrentV( );
   void ComputeJacobian( );
   void WeightsInit( );
   void BuildMatrix( );
   void BuildRightPart( );
   void Regularisation( );
   double CalcResidual( );
   void CalcIndividualResiduals( );
   void Solve( );

};

