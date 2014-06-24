default: parallel-intel

clean:
	rm -f *.i *.o *.s test_AS_ALN_bruteforcedp_parallel test_AS_ALN_bruteforcedp_sse test_AS_ALN_bruteforcedp

ref-gcc:
	g++ -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer -o test_AS_ALN_bruteforcedp AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp.C test_AS_ALN_bruteforcedp.cpp

ref-intel:
	icc -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer -o test_AS_ALN_bruteforcedp AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp.C test_AS_ALN_bruteforcedp.cpp

ref-profile:
	icc -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -g -O0 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer -o test_AS_ALN_bruteforcedp AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp.C test_AS_ALN_bruteforcedp.cpp

sse-gcc:
	g++ -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer -o test_AS_ALN_bruteforcedp_sse AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_sse.C test_AS_ALN_bruteforcedp_sse.cpp

sse-intel:
	icc -m64 -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer -o test_AS_ALN_bruteforcedp_sse AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_sse.C test_AS_ALN_bruteforcedp_sse.cpp

sse-intel-debug:
	icc -m64 -pthread -g3 -O0 -save-temps -o test_AS_ALN_bruteforcedp_sse AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_sse.C test_AS_ALN_bruteforcedp_sse.cpp

sse-profile:
	icc -m64 -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -g -O0 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer -o test_AS_ALN_bruteforcedp_sse AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_sse.C test_AS_ALN_bruteforcedp_sse.cpp

parallel-gcc:
	g++ -fopenmp -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer -o test_AS_ALN_bruteforcedp_parallel AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_parallel.C test_AS_ALN_bruteforcedp_parallel.cpp

parallel-intel:
	icc -m64 -pthread -fopenmp -Wall -Wno-write-strings -Wno-unused  -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -funroll-loops -fomit-frame-pointer -finline-functions -o test_AS_ALN_bruteforcedp_parallel AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_parallel.C test_AS_ALN_bruteforcedp_parallel.cpp

parallel-intel-debug:
	icc -m64 -pthread -fopenmp -g3 -O0 -save-temps -o test_AS_ALN_bruteforcedp_parallel AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_parallel.C test_AS_ALN_bruteforcedp_parallel.cpp

parallel-profile:
	icc -m64 -pthread -fopenmp -Wall -Wno-write-strings -Wno-unused -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -g -O0 -funroll-loops  -o test_AS_ALN_bruteforcedp_parallel AS_UTL_alloc.C AS_UTL_reverseComplement.C AS_ALN_bruteforcedp_parallel.C test_AS_ALN_bruteforcedp_parallel.cpp
