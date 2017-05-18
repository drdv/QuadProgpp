// Copyright (C) 2017 Alexander Sherikov
//
// This software may be modified and distributed under the terms
// of the MIT license.  See the LICENSE file for details.

#pragma once

// Utility function for computing the scalar product
template<typename T>
inline double scalar_product(const QPPP_VECTOR(T)& x, const QPPP_VECTOR(T)& y)
{
    return (x.transpose()*y);
}


template<typename T>
void print_matrix(const char* name, const QPPP_MATRIX(T)& A, int n = -1, int m = -1)
{
  if (n == -1)
    n = A.rows();
  if (m == -1)
    m = A.cols();

  std::cout << name << ": " << std::endl << " ";
  std::cout << A.block(0, 0, n, m) << std::endl;
}


template<typename T>
void print_vector(const char* name, const QPPP_VECTOR(T)& v, int n = -1)
{
  if (n == -1)
    n = v.rows();

  std::cout << name << ": " << std::endl << " ";
  std::cout << v.head(n) << std::endl;
}
