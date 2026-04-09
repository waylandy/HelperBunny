# build image
docker build -t cdd-search -f application/Dockerfile application
docker run -it -v $(pwd):/home --entrypoint bash cdd-search

# download database if not there
./download_database.sh

# run container
docker run -it -w /home -u $(id -u):$(id -g) \
    -v $(pwd):/home \
    -v ./database/data:/tmp/rpsbproc/data \
    -v ./database/db:/tmp/rpsbproc/db \
    -e INPUT_FASTA=examples/sequences.fasta \
    -e OUTPUT_DIR=examples/cdd_output \
    -e NUM_THREADS=4 \
    cdd-search

