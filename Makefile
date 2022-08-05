FLAGS=-DRSB_CONST_MAX_SUPPORTED_THREADS=128 -DRSB_WANT_VERBOSE_MESSAGES=0
LINKING=-I

build:
	nvcc ${FLAGS} librsb/examples/hello_cuda.cu librsb/*.c librsb/*.cu -o librsb/examples/hello_cuda ${LINKING}

run:
	./librsb/examples/hello_cuda
