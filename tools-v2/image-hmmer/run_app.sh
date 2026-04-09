#!/bin/bash

### build 
docker build -t hmmer -f application/Dockerfile application
docker run -it -v $(pwd):/home --entrypoint bash hmmer

### run container
docker run -it -w /home -u $(id -u):$(id -g) \
    -v $(pwd):/home \
    -e INPUT_FASTA=examples/sequences.fasta \
    -e PROFILE_HMM=examples/profile.hmm \
    -e OUTPUT_DIR=examples/output \
    -e BUFFER_SIZE=300000 \
    -e E_VALUE=1e-4 \
    -e DATABASE_SIZE=0 \
    -e N_THREADS=3 \
    hmmer /tmp/my_app/run.py

