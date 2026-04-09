FROM ubuntu:22.04

RUN apt-get update && apt-get install -y git wget unzip htop
RUN apt-get update && apt-get install -y python3.10 python3-pip 
RUN pip3 install --upgrade pip
RUN pip3 install numpy==1.26.3 pandas==2.2.3 polars==1.23.0 pyarrow fastparquet

RUN apt-get update && apt-get install -y ncbi-blast+ libdw-dev

RUN mkdir /tmp/my_app
ADD . /tmp/my_app

RUN tar -xzvf /tmp/my_app/src/RpsbProc-x64-linux.tar.gz -C /tmp \
 && mv /tmp/RpsbProc-x64-linux /tmp/rpsbproc \
 && ln -s /tmp/rpsbproc/rpsbproc /usr/bin/rpsbproc \
 && mkdir /tmp/rpsbproc/data /tmp/rpsbproc/db

ENTRYPOINT ["python3", "/tmp/my_app/process.py"]
