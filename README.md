This repository contains QuadProg++ solver with my modifications:

* added optional support for Eigen, which results in a better performance
  almost identical to eiQuadProg (which is under more restrictive GPL license);

* added Octave interface (it should work with MATLAB too);

* a number of minor optimizations and refactoring.

You are free to use this version of the solver under the terms of the MIT
license.

I am not planning to support or develop QuadProg++ any further - I am going to
focus on my own implementation of the algorithm instead:
> https://github.com/asherikov/qpmad

You can find the contents of the original README file below.

-----

# QuadProg++

A C++ library for Quadratic Programming which implements the
[Goldfarb-Idnani active-set dual method](http://www.javaquant.net/papers/goldfarbidnani.pdf).

At present it is limited to the solution of strictly convex quadratic programs.

Previous versions of the project were hosted on
[sourceforge](https://sourceforge.net/projects/quadprog/?source=directory).

## Install

To build the library simply go through the `cmake .; make; make install` cycle.

In order to use it, you will be required to include in your code file the
`Array.hh` header, which contains a handy C++ implementation of Vectors and
Matrices.

## Contribution

Contributions and bug fixes are welcome.

Copyright (C) 2007-2016 Luca Di Gaspero, [MIT License](LICENSE).
