FROM ubuntu:22.04

USER root

RUN apt-get update && apt-get install -y git wget unzip htop
RUN apt-get update && apt-get install -y python3.10 python3-pip 
RUN pip3 install --upgrade pip

RUN pip3 install numpy==1.26.3 pandas==2.2.3 polars==1.23.0 pyarrow fastparquet
RUN pip3 install biopython==1.86 pyhmmer==0.12.0

RUN pip3 install jupyterlab ipython matplotlib seaborn

RUN apt-get update && apt-get install -y hmmer

RUN mkdir /tmp/my_app
ADD . /tmp/my_app

ENTRYPOINT ["python3"]
