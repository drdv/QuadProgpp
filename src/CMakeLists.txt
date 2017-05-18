# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

if(QUADPROGPP_ENABLE_EIGEN)
    set(sources     QuadProg++.cc)
    set(headers     EigenHelpers.hh QuadProg++.hh config.hh)
else(QUADPROGPP_ENABLE_EIGEN)
    set(sources     Array.cc QuadProg++.cc)
    set(headers     Array.hh ArrayHelpers.hh QuadProg++.hh config.hh)
endif(QUADPROGPP_ENABLE_EIGEN)

set(library_sources ${sources})
list(REMOVE_ITEM library_sources ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_library(quadprog STATIC ${library_sources} ${headers})
set_property(TARGET quadprog PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET quadprog PROPERTY PUBLIC_HEADER ${headers})

if(QUADPROGPP_BUILD_EXAMPLE)
    add_executable(main main.cc ${headers})
    target_link_libraries(main quadprog)
endif(QUADPROGPP_BUILD_EXAMPLE)

install(TARGETS quadprog
  EXPORT quadprog-targets
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/QuadProg++
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(EXPORT quadprog-targets DESTINATION cmake)
