// Copyright (C) 2001--2008 Andrea Schaerf, Luca Di Gaspero.
// Copyright (C) 2017 Alexander Sherikov
//
// This software may be modified and distributed under the terms
// of the MIT license.  See the LICENSE file for details.

#pragma once


#include <sstream> // print_matrix, print_vector

namespace QuadProgpp
{

// Utility function for computing the scalar product
template<typename T>
inline double scalar_product(const Vector<T>& x, const Vector<T>& y)
{
    register int i, n = x.size();
    register double sum;

    sum = 0.0;
    for (i = 0; i < n; i++)
        sum += x[i] * y[i];
    return sum;
}


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

} //namespace QuadProgpp
