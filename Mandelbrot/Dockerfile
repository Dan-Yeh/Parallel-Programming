FROM nlknguyen/alpine-mpich

RUN sudo apk add --update --no-cache libgomp

USER root

RUN mkdir /parallel
WORKDIR /parallel
COPY ./ /parallel

# docker build -t parallel . && docker run --rm parallel sh -c "make && mpirun -np4 /parallel/student_submission"