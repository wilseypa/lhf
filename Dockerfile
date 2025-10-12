# syntax=docker/dockerfile:1

FROM ubuntu:latest

# Install build dependencies
RUN apt-get update && \
    apt-get install -y \
        cmake \
        git \
        build-essential \
        libomp-dev \
        libeigen3-dev \
        openmpi-bin openmpi-common libopenmpi-dev \
        libcgal-dev \
        python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Set working directory inside container
WORKDIR /lhf

# Copy source code into container (optional)
COPY . .

# Start an interactive shell by default
CMD ["bash"]

#docker build -t lhf-container .

#docker run -it --rm --name lhf-container -v "$(pwd)":/lhf -w /lhf lhf-container bash

#docker exec -it -w /lhf lhf-container bash

#./LHF -i inputData.csv -b 1 -bm betaHighCircle -e 100 -d 2 -complexType betaPoly -m beta

