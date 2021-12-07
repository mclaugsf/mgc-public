FROM ubuntu:latest

# To build:
#docker build -t medicinalgenomics/megahit .

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y --no-install-recommends \
	git \
    bzip2 \
    cmake \
    gzip \
    g++ \
    libgomp1 \
    make \
    python

RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
    python3

RUN apt-get update && apt-get install -y --no-install-recommends --reinstall ca-certificates

RUN cd /opt && \
	git clone https://github.com/voutcn/megahit.git && \
	cd megahit && \
	git submodule update --init && \
	mkdir build && cd build && \
	cmake .. -DCMAKE_BUILD_TYPE=Release && \
	make -j4 && \
	make simple_test && \
	make install
