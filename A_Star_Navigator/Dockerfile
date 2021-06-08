FROM alpine:latest

RUN apk add --update --no-cache g++ libgomp make ffmpeg

RUN mkdir /parallel
WORKDIR /parallel
COPY ./ /parallel

# docker build -t parallel . && docker run --rm parallel sh -c "make && /parallel/student_submission"