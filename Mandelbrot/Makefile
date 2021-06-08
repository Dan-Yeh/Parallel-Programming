CXX = mpicxx
CXX_FLAGS = --std=c++17 -Wall -Wextra -mavx -O3 -g -DOMPI_SKIP_MPICXX -fopenmp
# this compiler definition is needed to silence warnings caused by the openmpi CXX
# bindings that are deprecated. This is needed on gcc 8 forward.
# see: https://github.com/open-mpi/ompi/issues/5157

all: env_file student_submission sequential_implementation

reference_implementation:
	$(CXX) $(CXX_FLAGS) -o reference_implementation reference_implementation.cpp

sequential_implementation:
	$(CXX) $(CXX_FLAGS) -o sequential_implementation sequential_implementation.cpp

student_submission:
	$(CXX) $(CXX_FLAGS) -o student_submission student_submission.cpp

env_file: student_submission.cpp
	sed -n -e 's/.*!submission_env \(.*\)/\1/p' student_submission.cpp > student_submission.env

clean:
	rm -f sequential_implementation student_submission sequential_implementation student_submission.env mandelbrot.ppm
