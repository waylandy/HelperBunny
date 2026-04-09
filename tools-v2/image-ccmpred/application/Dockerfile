FROM nvidia/cuda:12.4.0-devel-ubuntu22.04

# sudo chown -R $USER:$USER . && sudo chmod -R 755 .
# docker build -t ccmpred . && docker run -it -w /home -v $(pwd):/home --gpus all --entrypoint bash ccmpred

# docker build -t ccmpred -f application/Dockerfile application && docker run -it -w /home -v $(pwd):/home --gpus all --entrypoint bash ccmpred

##################

RUN apt-get update && apt-get install -y git wget unzip htop
RUN apt-get update && apt-get install -y python3.10 python3-pip 
RUN pip3 install --upgrade pip
RUN pip3 install numpy==1.26.3 pandas==2.2.3 polars==1.23.0 pyarrow fastparquet

RUN mkdir /tmp/my_app
ADD . /tmp/my_app

# RUN apt-get update && apt-get install -y cmake
# RUN unzip /tmp/my_app/ccmpred/ccmpred.zip -d /tmp \
#  && cd /tmp/CCMpred \
#  && cmake . \
#  && make \
#  && make install
RUN unzip /tmp/my_app/ccmpred/ccmpred-binary.zip -d /usr/local/bin/

ENTRYPOINT ["python3"]
