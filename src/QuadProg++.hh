/*
 File $Id: QuadProg++.hh 232 2007-06-21 12:29:00Z digasper $

 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani
 for the solution of a (convex) Quadratic Programming problem
 by means of an active-set dual method.

The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^T x + ce0 = 0
    CI^T x + ci0 >= 0

 The matrix and vectors dimensions are as follows:
     G: n * n
		g0: n

		CE: n * p
	 ce0: p

	  CI: n * m
   ci0: m

     x: n

 The function will return the cost of the solution written in the x vector or
 std::numeric_limits::infinity() if the problem is infeasible. In the latter case
 the value of the x vector is not correct.

 References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

 Notes:
  1. pay attention in setting up the vectors ce0 and ci0.
	   If the constraints of your problem are specified in the form
	   A^T x = b and C^T x >= d, then you should set ce0 = -b and ci0 = -d.

TODO deprecated?
  2. The matrices have column dimension equal to MATRIX_DIM,
     a constant set to 20 in this file (by means of a #define macro).
     If the matrices are bigger than 20 x 20 the limit could be
		 increased by means of a -DMATRIX_DIM=n on the compiler command line.

  3. The matrix G is modified within the function since it is used to compute
     the G = L^T L cholesky factorization for further computations inside the function.
     If you need the original matrix G you should make a copy of it and pass the copy
     to the function.

 Author: Luca Di Gaspero
  			 DIEGM - University of Udine, Italy
				 luca.digaspero@uniud.it
				 http://www.diegm.uniud.it/digaspero/

 The author will be grateful if the researchers using this software will
 acknowledge the contribution of this function in their research papers.

 Copyright (c) 2007-2016 Luca Di Gaspero

 This software may be modified and distributed under the terms
 of the MIT license.  See the LICENSE file for details.
*/


#ifndef _QUADPROGPP
#define _QUADPROGPP

#include "config.hh"


#ifdef QUADPROGPP_ENABLE_EIGEN

#include <Eigen/Dense>

#define QPPP_VECTOR(t_Scalar)   Eigen::Matrix<t_Scalar, Eigen::Dynamic, 1, Eigen::ColMajor>
#define QPPP_MATRIX(t_Scalar)   Eigen::Matrix<t_Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>

#else

#include "Array.hh"

#define QPPP_VECTOR(t_Scalar)   QuadProgpp::Vector<t_Scalar>
#define QPPP_MATRIX(t_Scalar)   QuadProgpp::Matrix<t_Scalar>

#endif


namespace QuadProgpp
{
    class Status
    {
        public:
            enum Value
            {
                OK = 0,
                FAILURE = 1
            };
    };


    class Solver
    {
        public:
            Solver();
            ~Solver();

            Status::Value solve(   QPPP_MATRIX(double)& G,
                            QPPP_VECTOR(double)& g0,
                            const QPPP_MATRIX(double)& CE,
                            const QPPP_VECTOR(double)& ce0,
                            const QPPP_MATRIX(double)& CI,
                            const QPPP_VECTOR(double)& ci0,
                            QPPP_VECTOR(double)& x);

        private:
            class Implementation;
            Implementation *impl;
    };
} // namespace QuadProgpp

#endif // #define _QUADPROGPP
