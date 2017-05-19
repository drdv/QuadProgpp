octave: forceclean
	octave --no-history --silent make.m

matlab: forceclean
	matlab -nosplash -nodesktop -nojvm -r 'make; quit' | tail +10

octave_test: octave
	octave --no-history --silent qppp_test.m

clean:
	rm -f *.o

forceclean: clean
	rm -f *.mex*
