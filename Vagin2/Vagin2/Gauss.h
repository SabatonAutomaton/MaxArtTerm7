#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

class Gauss {
public:
   vector<double> solve( vector<vector<double>> A, vector<double> b );
   void printMatrix( const vector<vector<double>> &A, const vector<double> &b );
};