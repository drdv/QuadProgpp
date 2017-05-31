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


        void invert_upper(  const QPPP_MATRIX(T)& A,
                            QPPP_MATRIX(T)& J)
        {
            int n = A.rows();

            // in-place inversion of upper triangular A
            for (int i = n-1; i >=0; --i)
            {
                // zeros below the main diagonal
                for(int k = n-1; k > i; --k)
                {
                    J(k, i) = 0.0;
                }

                // main diagonal
                J(i,i) = 1.0 / A(i, i);

                // elements above the main diagonal
                for (int j = i-1; j >= 0; --j)
                {
                    double tmp = 0.0;
                    for (int k = j+1; k <= i; ++k)
                    {
                        tmp -= A(k,j) * J(k,i);
                    }
                    J(j, i) = tmp / A(j,j);
                }
            }


            // Original version -- requires extra temporary vectors and should
            // be slower
            /*
            int n = A.rows();
            d.setZero(n);

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
            */
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

        void solve_negative(QPPP_MATRIX(T)& A, QPPP_VECTOR(T)& x, const QPPP_VECTOR(T)& b)
        {
            solve(A, x, b);
            for (int i = 0; i < x.size(); ++i)
            {
                x[i] = -x[i];
            }
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
inline void multiply_and_add(   QPPP_VECTOR(T)& y,
                                const QPPP_MATRIX(T)& A,
                                const QPPP_VECTOR(T)& x,
                                const QPPP_VECTOR(T)& b)
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


template<typename T>
inline void multiply_and_add_i( QPPP_VECTOR(T)& y,
                                const QPPP_MATRIX(T)& A,
                                const QPPP_VECTOR(T)& x,
                                const QPPP_VECTOR(T)& b,
                                const int index)
{
    double sum = 0.0;
    for (int k = 0; k < x.size(); ++k)
        sum += A(k, index) * x[k];
    y[index] = sum + b[index];
}


//-----------------------------------------------------------------------
// Utility functions for updating some data needed by the solution method
//-----------------------------------------------------------------------
template<typename T>
inline void compute_d(QPPP_VECTOR(T)& d, const QPPP_MATRIX(T)& J, const QPPP_VECTOR(T)& np)
{
    register int i, j, n = d.size();
    register double sum;

    /* compute d = H^T * np */
    for (i = 0; i < n; i++)
    {
        sum = 0.0;
        for (j = 0; j < n; j++)
            sum += J(j, i) * np[j];
        d[i] = sum;
    }
}


template<typename T>
inline void compute_primal_step_direction(QPPP_VECTOR(T)& z, const QPPP_MATRIX(T)& J, const QPPP_VECTOR(T)& d, const int iq)
{
    register int i, j, n = z.size();

    /* setting of z = H * d */
    for (i = 0; i < n; i++)
    {
        z[i] = 0.0;
        for (j = iq; j < n; j++)
            z[i] += J(i, j) * d[j];
    }
}


template<typename T>
inline void compute_dual_step_direction(const QPPP_MATRIX(T)& R, QPPP_VECTOR(T)& r, const QPPP_VECTOR(T)& d, const int iq)
{
    register int i, j, n = d.size();
    register double sum;

    /* setting of r = R^-1 d */
    for (i = iq - 1; i >= 0; i--)
    {
        sum = 0.0;
        for (j = i + 1; j < iq; j++)
            sum += R(i, j) * r[j];
        r[i] = (d[i] - sum) / R(i, i);
    }
}
} //namespace QuadProgpp
