BUILD_DIR=build

test_build:
	mkdir -p build;
	cd build; cmake -DQUADPROGPP_BUILD_EXAMPLE=ON -DQUADPROGPP_ENABLE_TRACING=OFF -DQUADPROGPP_ENABLE_EIGEN=OFF ..;
	cd build; ${MAKE}
	rm -Rf build/*
	cd build; cmake -DQUADPROGPP_BUILD_EXAMPLE=ON -DQUADPROGPP_ENABLE_TRACING=ON -DQUADPROGPP_ENABLE_EIGEN=OFF ..;
	cd build; ${MAKE}
	rm -Rf build/*
	cd build; cmake -DQUADPROGPP_BUILD_EXAMPLE=ON -DQUADPROGPP_ENABLE_TRACING=ON -DQUADPROGPP_ENABLE_EIGEN=ON ..;
	cd build; ${MAKE}
	rm -Rf build/*
	cd build; cmake -DQUADPROGPP_BUILD_EXAMPLE=ON -DQUADPROGPP_ENABLE_TRACING=OFF -DQUADPROGPP_ENABLE_EIGEN=ON ..;
	cd build; ${MAKE}

clean:
	rm -Rf build;
	cd matlab_octave; ${MAKE} clean
