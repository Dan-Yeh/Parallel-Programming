all: sequential_implementation student_submission

env_file: student_submission.cpp
	sed -n -e 's/.*!submission_env \(.*\)/\1/p' student_submission.cpp > student_submission.env

reference_implementation: reference_implementation.cpp
	g++ -std=c++17 -Wall -Wextra -o reference_implementation -O3 -g reference_implementation.cpp -lpthread -lm

student_submission: student_submission.cpp
<<<<<<< HEAD
	g++ -std=c++17 -Wall -Wextra -o student_submission -O3 -g student_submission.cpp -lpthread -lm

sequential_implementation: sequential_implementation.cpp
	g++ -std=c++17 -Wall -Wextra -o sequential_implementation -O3 -g sequential_implementation.cpp -lpthread -lm
=======
	g++ -std=c++17 -Wall -Wextra -o student_submission -O3 -pg -g student_submission.cpp -lpthread -lm

sequential_implementation: sequential_implementation.cpp
	g++ -std=c++17 -Wall -Wextra -o sequential_implementation -O3 -pg -g sequential_implementation.cpp -lpthread -lm
>>>>>>> f8b52302b29d8ea022ffda4d4834c084e219aa41

clean:
	rm -f sequential_implementation reference_implementation student_submission

