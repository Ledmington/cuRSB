FILE=matrix_conversions

build:
	gcc `librsb-config --I_opts` -c src/${FILE}.c -o bin/${FILE}.o
	gcc -o bin/${FILE} ${FILE}.o `librsb-config --static --ldflags --extra_libs`

run:
	./bin/${FILE}

clean:
	rm bin/*
