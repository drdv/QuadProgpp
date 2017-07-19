all: matlab install

matlab:
	matlab -nosplash -nodesktop -nojvm -r 'make; quit' | tail +10

octave:
	octave --no-history --silent make.m

octave_test: octave
	octave --no-history --silent qppp_test.m

install:
	mv *.mex* ~/local/bin/qp_mex
	cp quadprogpp.m ~/local/bin/qp_mex

clean:
	rm -f *.o
	rm -f *.mex*
