% ------------------------------------------------------------
% compilation
% ------------------------------------------------------------

copyfile('config.hh', '../src');
mex -v -I../src/   quadprogpp_interface.cpp  ../src/Array.cc ../src/QuadProg++.cc
%mex -v -I../src/ -I/usr/local/include/eigen3/  quadprogpp_interface.cpp  ../src/Array.cc ../src/QuadProg++.cc

% note: when I define QUADPROGPP_ENABLE_EIGEN, I get unexpected results
%       on the test qpgi/matlab/other_test_problems/test_problem_005.m

%%%EOF
