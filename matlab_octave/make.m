% ------------------------------------------------------------
% compilation
% ------------------------------------------------------------

copyfile('config.hh', '../src');
mex -v -I../src/   quadprogpp_interface.cpp  ../src/Array.cc ../src/QuadProg++.cc
%mex -v -I../src/ -I/usr/local/include/eigen3/  quadprogpp_interface.cpp  ../src/Array.cc ../src/QuadProg++.cc

%%%EOF
