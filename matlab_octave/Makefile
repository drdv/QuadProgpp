octave: forceclean
	octave --no-history --silent make.m

matlab: forceclean
	matlab -nosplash -nodesktop -nojvm -r 'make; quit' | tail +10

clean:
	rm -f *.o

forceclean: clean
	rm -f *.mex*
