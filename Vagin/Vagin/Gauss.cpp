#include "Gauss.h"

vector<double> Gauss::solve( vector<vector<double>> A, vector<double> b ) {
   int n = A.size( );

   if ( n == 0 || A[0].size( ) != n || b.size( ) != n ) {
      throw invalid_argument( "Неверные размеры матрицы или вектора" );
   }

   for ( int k = 0; k < n; ++k ) {
      int maxRow = k;
      double maxVal = fabs( A[k][k] );

      for ( int i = k + 1; i < n; ++i ) {
         if ( fabs( A[i][k] ) > maxVal ) {
            maxVal = fabs( A[i][k] );
            maxRow = i;
         }
      }

      if ( maxVal < 1e-10 ) {
         throw runtime_error( "Division by zero" );
      }

      if ( maxRow != k ) {
         swap( A[k], A[maxRow] );
         swap( b[k], b[maxRow] );
      }

      double pivot = A[k][k];
      for ( int j = k; j < n; ++j ) {
         A[k][j] /= pivot;
      }
      b[k] /= pivot;

      for ( int i = k + 1; i < n; ++i ) {
         double factor = A[i][k];
         for ( int j = k; j < n; ++j ) {
            A[i][j] -= factor * A[k][j];
         }
         b[i] -= factor * b[k];
      }
   }

   vector<double> x( n, 0.0 );
   for ( int i = n - 1; i >= 0; --i ) {
      x[i] = b[i];
      for ( int j = i + 1; j < n; ++j ) {
         x[i] -= A[i][j] * x[j];
      }
   }

   return x;
}

void Gauss::printMatrix( const vector<vector<double>> &A, const vector<double> &b ) {
   int n = A.size( );
   for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < n; ++j ) {
         cout << A[i][j] << "\t";
      }
      cout << "| " << b[i] << endl;
   }
   cout << endl;
}
