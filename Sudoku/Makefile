CXX = g++
CXX_FLAGS = --std=c++17 -Wall -Wextra -fopenmp -O3

all: env_file sequential_implementation student_submission

sequential_implementation: sequential_implementation.cpp 
	$(CXX) $(CXX_FLAGS) -o sequential_implementation sequential_implementation.cpp

reference_implementation: reference_implementation.cpp
	$(CXX) $(CXX_FLAGS) -o reference_implementation reference_implementation.cpp

student_submission: student_submission.cpp
	$(CXX) $(CXX_FLAGS) -o student_submission student_submission.cpp

env_file: student_submission.cpp
	sed -n -e 's/.*!submission_env \(.*\)/\1/p' student_submission.cpp > student_submission.env

clean:
	rm -f sequential_implementation reference_implementation student_submission student_submission.env
