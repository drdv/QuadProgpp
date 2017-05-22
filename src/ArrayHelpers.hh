// Copyright (C) 2001--2008 Andrea Schaerf, Luca Di Gaspero.
// Copyright (C) 2017 Alexander Sherikov
//
// This software may be modified and distributed under the terms
// of the MIT license.  See the LICENSE file for details.

#pragma once


#include <sstream> // print_matrix, print_vector

namespace QuadProgpp
{

template<typename T>
void print_matrix(const char* name, const Matrix<T>& A, int n = -1, int m = -1)
{
    std::ostringstream s;
    std::string t;
    if (n == -1)
        n = A.rows();
    if (m == -1)
        m = A.cols();

    s << name << ": " << std::endl;
    for (int i = 0; i < n; i++)
    {
        s << " ";
        for (int j = 0; j < m; j++)
            s << A[i][j] << ", ";
        s << std::endl;
    }
    t = s.str();
    t = t.substr(0, t.size() - 3); // To remove the trailing space, comma and newline

    std::cout << t << std::endl;
}


template<typename T>
void print_vector(const char* name, const Vector<T>& v, int n = -1)
{
    std::ostringstream s;
    std::string t;
    if (n == -1)
        n = v.size();

    s << name << ": " << std::endl << " ";
    for (int i = 0; i < n; i++)
    {
        s << v[i] << ", ";
    }
    t = s.str();
    t = t.substr(0, t.size() - 2); // To remove the trailing space and comma

    std::cout << t << std::endl;
}


// Utility functions for computing the Cholesky decomposition and solving
// linear systems
template<typename T>
class CholeskyDecomposition
{
    public:
        void compute (QPPP_MATRIX(T)& A)
        {
            register int i, j, k, n = A.rows();
            register double sum;

            for (i = 0; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    sum = A(i, j);
                    for (k = i - 1; k >= 0; k--)
                        sum -= A(i, k)*A(j, k);
                    if (i == j)
                    {
                        if (sum <= 0.0)
                        {
                            std::ostringstream os;
                            // raise error
#ifdef QUADPROGPP_ENABLE_TRACING
                            print_matrix("A", A);
#endif
                            os << "Error in cholesky decomposition, sum: " << sum;
                            throw std::logic_error(os.str());
                            exit(-1);
                        }
                        A(i, i) = std::sqrt(sum);
                    }
                    else
                        A(j, i) = sum / A(i, i);
                }
                for (k = i + 1; k < n; k++)
                    A(i, k) = A(k, i);
            }
        }


        void invert_upper(QPPP_MATRIX(T)& A, QPPP_MATRIX(T)& J, QPPP_VECTOR(T)& z, QPPP_VECTOR(T)& d)
        {
            int n = A.rows();

            for (int i = 0; i < n; i++)
            {
                d[i] = 1.0;
                forward_elimination(A, z, d);
                for (int j = 0; j < n; j++)
                {
                    J(i, j) = z[j];
                }
                d[i] = 0.0;
            }
        }


        void solve(QPPP_MATRIX(T)& A, QPPP_VECTOR(T)& x, const QPPP_VECTOR(T)& b)
        {
            int n = A.rows();
            QPPP_VECTOR(T) y(n);

            /* Solve L * y = b */
            forward_elimination(A, y, b);
            /* Solve L^T * x = y */
            backward_elimination(A, x, y);
        }


    private:
        inline void forward_elimination(QPPP_MATRIX(T)& A, QPPP_VECTOR(T)& y, const QPPP_VECTOR(T)& b)
        {
            register int i, j, n = A.rows();

            y[0] = b[0] / A(0, 0);
            for (i = 1; i < n; i++)
            {
                y[i] = b[i];
                for (j = 0; j < i; j++)
                    y[i] -= A(i, j) * y[j];
                y[i] = y[i] / A(i, i);
            }
        }

        inline void backward_elimination(QPPP_MATRIX(T)& A, QPPP_VECTOR(T)& x, const QPPP_VECTOR(T)& y)
        {
            register int i, j, n = A.rows();

            x[n - 1] = y[n - 1] / A(n - 1, n - 1);
            for (i = n - 2; i >= 0; i--)
            {
                x[i] = y[i];
                for (j = i + 1; j < n; j++)
                    x[i] -= A(i, j) * x[j];
                x[i] = x[i] / A(i, i);
            }
        }
};


template<typename T>
void multiply_and_add(QPPP_VECTOR(T)& y, const QPPP_MATRIX(T)& A, const QPPP_VECTOR(T)& x, const QPPP_VECTOR(T)& b)
{
    for (int i = 0; i < A.cols(); i++)
    {
        double sum = 0.0;
        for (int j = 0; j < A.rows(); j++)
            sum += A(j, i) * x[j];
        sum += b[i];
        y[i] = sum;
    }
}
} //namespace QuadProgpp
