#include <iostream>
#include <Eigen/Dense> // Include the Eigen library for matrix operations
#include <cmath>// Include cmath for mathematical functions
#include <iomanip> // Includeiomanip for setting precision of output

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Function to perform the Gauss-Seidel method
VectorXd gaussSeidel(const MatrixXd&, const VectorXd&, int, double);



int main() {
    
    int p = 3;


    
    MatrixXd A(p, p);
    VectorXd b(p);

    /*
  
    - Static initialization of the matrix
    The matrix m is a 5x5 matrix with the following values with the last column being the result of the sum of each row:
    -7, -4, -1, -1, = -10
    -1,  8, -5,  0, = 35
    -1, -2, -4,  0, = 0
    2,  0,  0, -6; = 36    
    */

    // initialize the matrix with values
    // A << -7, -4, -1, -1,
    //      -1,  8, -5,  0,
    //      -1, -2, -4,  0,
    //       2,  0,  0, -6;

    // // initialize the vector with values
    // b << -10, 35, 0, -36;

    A << 5, 2,-1,
         2, 3, 4,
         3, 4, 2;

    b<<3, -1, 8;


    // Print the matrix A
    cout << "Matriz A (" << A.rows() << "x" << A.cols() << "):\n" << A << "\n\n";
    // Print the vector  b
    cout << "Vector b:\n" << b.transpose() << "\n\n";

    int iteraciones = 100; // Number of iterations for the Gauss-Seidel method

    // Execute  Gauss-Seidel Mehtod
    VectorXd result = gaussSeidel(A, b, iteraciones, 0.0001); // Call the Gauss-Seidel method
    cout << "Result:" << result.transpose() << std::endl; // Print the result


  





    return 0;
}


/*
    The Gauss-Seidel method is an iterative method for solving a system of linear equations of the form Ax = b.
    It is particularly useful for large sparse systems where direct methods may be inefficient.
    The method works by iteratively updating the solution vector x until it converges to the true solution.
    The convergence criteria can be based on a maximum number of iterations or a specified tolerance level.

    Parameters:
    - A: The coefficient matrix of the system of equations.
    - b: Vector of constants.
    - maxIter: The maximum number of iterations to perform.
    - tol: The convergence tolerance. The algorithm stops if the error is less than this value.

    Return:
        - x: The solution vector.

    */



   VectorXd gaussSeidel(const MatrixXd& A, const VectorXd& b, int maxIter, double tol) {

    int n = A.rows();
    VectorXd x = VectorXd::Ones(n);     // Actual solution
    VectorXd xOld = x;                  // previus Solution
    VectorXd error(n);

    cout << fixed << setprecision(6);
    cout << "Iter\t";
    
    // Print the header of the table
    for (int i = 0; i < n; ++i)
        cout << "x" << i+1 << "\t\t";
    for (int i = 0; i < n; ++i)
        cout << "err(x" << i+1 << ")\t\t";
    cout << "\n";

    // Main loop for Gauss-Seidel method
    for (int k = 1; k <= maxIter; ++k) {
        xOld = x;

        for (int i = 0; i < n; ++i) {
            double suma = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i)
                    suma += A(i,j) * x(j);
            }
            x(i) = (b(i) - suma) / A(i,i);
        }

        // Calculate relative error
        for (int i = 0; i < n; ++i)
            error(i) = fabs((x(i) - xOld(i)) / x(i));


        //print a row from the table
        cout << k << "\t";
        for (int i = 0; i < n; ++i) cout << x(i) << "\t";
        for (int i = 0; i < n; ++i) cout << error(i) << "\t";
        cout << "\n";

        // Check if the error is less than the tolerance level
        if (error.maxCoeff() <= tol) {
            cout << "\nThere is convergence after " << k << " iteration.\n";
            return x;
        }
    }

    cout << "\nThere was no convergence after " << maxIter << " iteration.\n";
    return x;
}



