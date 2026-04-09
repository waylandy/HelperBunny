# build 
docker build -t ccmpred -f application/Dockerfile application
docker run -it -v $(pwd):/home --entrypoint bash --gpus all ccmpred
# sudo chown -R $USER:$USER . && sudo chmod -R 755 .

# run ccmpred
docker run -it -w /home -u $(id -u):$(id -g) \
    --gpus all \
    -v $(pwd):/home \
    -e INPUT_A2M=examples/PF00017_gap20_head.a2m \
    -e OUTPUT_NPZ=examples/model.npz \
    ccmpred /tmp/my_app/get_model.py

# score
docker run -it -w /home -u $(id -u):$(id -g) \
    --gpus all \
    -v $(pwd):/home \
    -e INPUT_A2M=examples/PF00017_gap20_head.a2m \
    -e POTTS_NPZ=examples/model.npz \
    -e OUTPUT_FILE=examples/scores.parquet \
    ccmpred /tmp/my_app/get_scores.py


