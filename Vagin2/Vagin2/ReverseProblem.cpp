#include "ReverseProblem.h"
#include "Gauss.h"

double ReverseProblem::Distance(xyz p1, xyz p2)
{
   return sqrt((p2.x - p1.x) * (p2.x - p1.x) +
      (p2.y - p1.y) * (p2.y - p1.y) +
      (p2.z - p1.z) * (p2.z - p1.z));
}

void ReverseProblem::ComputeTrueV()
{
   double piSigma = 0.5 * 1 / (std::numbers::pi * trueSigma);
   for (int i = 0; i < n; i++)
   {
      double geom = 1.0 / Distance(sourceBPos[1], receiverMPos[i])
         - 1.0 / Distance(sourceAPos[1], receiverMPos[i])
         - (1.0 / Distance(sourceBPos[1], receiverNPos[i])
            - 1.0 / Distance(sourceAPos[1], receiverNPos[i]));
      trueV[i] = piSigma * geom;
      double noise = 0.1 * (2.0 * rand() / RAND_MAX - 1.0);// случайное в [-0.1, 0.1]
      //trueV[i] *= (1.0 + noise);
   }
}

void ReverseProblem::ComputeCurrentV()
{
   sigma = std::exp(lnSigma);
   double piSigma = 0.5 * 1 / (std::numbers::pi * sigma);
   for (int i = 0; i < n; i++)
   {
      double geom = 1.0 / Distance(sourceBPos[1], receiverMPos[i])
         - 1.0 / Distance(sourceAPos[1], receiverMPos[i])
         - (1.0 / Distance(sourceBPos[1], receiverNPos[i])
            - 1.0 / Distance(sourceAPos[1], receiverNPos[i]));

      V[i] = piSigma * geom;
   }
}

void ReverseProblem::ComputeJacobian()
{
   sigma = std::exp(lnSigma);
   ComputeCurrentV();
   for (int i = 0; i < n; i++)
      dVdSigma[i] = -V[i];
   double piSigma = 0.5 * 1 / (std::numbers::pi * sigma);
   for (int i = 0; i < n; i++)
   {
      double geom = 1.0 / Distance(sourceBPos[1], receiverMPos[i])
         - 1.0 / Distance(sourceAPos[1], receiverMPos[i])
         - (1.0 / Distance(sourceBPos[1], receiverNPos[i])
            - 1.0 / Distance(sourceAPos[1], receiverNPos[i]));

      // dV_i/dσ = - geom / (2π σ²)
   }
}

void ReverseProblem::WeightsInit()
{
   for (int i = 0; i < n; i++)
      w2[i] = 1;
}

void ReverseProblem::BuildMatrix()
{
   A[0][0] = 0.0;
   for (int i = 0; i < n; i++)
      A[0][0] += w2[i] * dVdSigma[i] * dVdSigma[i];
}

void ReverseProblem::BuildRightPart()
{
   b[0] = 0.0;
   for (int i = 0; i < n; i++)
      b[0] += w2[i] * (trueV[i] - V[i]) * dVdSigma[i];
}

//void ReverseProblem::Regularisation()
//{
//   for (int i = 0; i < n; i++)
//   {
//      A[i][i] += alpha;
//      b[i] -= isApriori[i] * alpha * (I[i] - aprioriI[i]);
//   }
//}

double ReverseProblem::CalcResidual()
{
   double res = 0;
   for (int i = 0; i < n; i++)
   {
      double diff = V[i] - trueV[i];
      res += w2[i] * diff * diff;
   }

   //// Добавляем регуляризационный член
   //for ( int j = 0; j < n; j++ )
   //{
   //   res += alpha * I[j] * I[j];
   //}

   return res;
}

void ReverseProblem::CalcIndividualResiduals()
{
   for (int i = 0; i < n; i++)
   {
      double diff = V[i] - trueV[i];
      residuals[i] = w2[i] * diff * diff;
   }
}

void ReverseProblem::Solve()
{
   double residual = 1;
   Gauss gauss;
   
   ComputeTrueV();
   WeightsInit();
   ComputeJacobian();

   std::cout << "Iter\tsigma\tTotal_Residual\tPhi1\tPhi2\tPhi3" << std::endl;

   // Начальное приближение sigma задано в конструкторе (например, 0.01)
   ComputeCurrentV();
   residual = CalcResidual();
   CalcIndividualResiduals();

   std::cout << 0 << "\t" << sigma << "\t" << residual << "\t";
   for (int i = 0; i < n; i++)
   {
      std::cout << residuals[i];
      if (i < n - 1) std::cout << "\t";
   }
   std::cout << std::endl;

   for (int iter = 0; iter < maxIter && residual > eps; iter++)
   {
      A[0][0] = 0.0;
      b[0] = 0.0;

      ComputeCurrentV();
      ComputeJacobian();
      BuildMatrix();
      BuildRightPart();
      //gauss.printMatrix( A, b );
      //Regularisation();

      std::vector<double> delta = gauss.solve(A, b);
      if (std::abs(delta[0]) < 1e-6) break;
      lnSigma += delta[0];
      sigma = std::exp(lnSigma);

      ComputeCurrentV();
      residual = CalcResidual();
      CalcIndividualResiduals();
      std::cout << iter + 1 << "\t" << sigma << "\t" << residual << "\t";
      for (int i = 0; i < n; i++)
      {
         std::cout << residuals[i];
         if (i < n - 1) std::cout << "\t";
      }
      std::cout << std::endl;
   }
}

