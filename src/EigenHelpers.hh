// Copyright (C) 2017 Alexander Sherikov
//
// This software may be modified and distributed under the terms
// of the MIT license.  See the LICENSE file for details.

#pragma once

namespace QuadProgpp
{

template<typename T>
void print_matrix(const char* name, const QPPP_MATRIX(T)& A, int n = -1, int m = -1)
{
    if (n == -1)
    {
        n = A.rows();
    }
    if (m == -1)
    {
        m = A.cols();
    }

    std::cout << name << ": " << std::endl << " ";
    std::cout << A.block(0, 0, n, m) << std::endl;
}


template<typename T>
void print_vector(const char* name, const QPPP_VECTOR(T)& v, int n = -1)
{
    if (n == -1)
    {
        n = v.rows();
    }

    std::cout << name << ": " << std::endl << " ";
    std::cout << v.head(n) << std::endl;
}


// Utility functions for computing the Cholesky decomposition and solving
// linear systems
template<typename T>
class CholeskyDecomposition
{
    public:
        template<class t_Derived>
        void compute (Eigen::PlainObjectBase<t_Derived>& A)
        {
            lltOfA.compute(A);
        }


        template<class t_Derived>
        void solve(Eigen::PlainObjectBase<t_Derived>& A, QPPP_VECTOR(T)& x, const QPPP_VECTOR(T)& b)
        {
            x = lltOfA.solve(b);
        }


        template<   class t_Derived_J,
                    class t_Derived_A>
        void invert_upper(  const Eigen::PlainObjectBase<t_Derived_A> & A,
                            Eigen::PlainObjectBase<t_Derived_J> & J)
        {
            // this is slightly slower
            /*
            J.setIdentity();
            lltOfA.matrixU().solveInPlace( J );
            */


            int n = A.rows();
            // a copy is needed since row-wise access to triangular matrix
            // matrixU is not supported
            J = lltOfA.matrixU();

            // in-place inversion of upper triangular J
            for (int i = n-1; i >=0; --i)
            {
                J(i,i) = 1.0 / J(i, i);
                for (int j = i-1; j >= 0; --j)
                {
                    double tmp = J.row(j).segment(j+1, i-j)*J.col(i).segment(j+1, i-j);
                    J(j, i) = - tmp / J(j,j);
                }
            }
        }


    private:
        Eigen::LLT< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > lltOfA;
};


template<typename T>
void multiply_and_add(QPPP_VECTOR(T)& y, const QPPP_MATRIX(T)& A, const QPPP_VECTOR(T)& x, const QPPP_VECTOR(T)& b)
{
    y.noalias() = A.transpose()*x + b;
}
} // namesace QuadProgpp
