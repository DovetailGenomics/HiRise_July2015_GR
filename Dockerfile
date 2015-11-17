FROM ubuntu:14.04

# Install required packages
RUN \
  apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    liblapack-dev \
    gfortran \
    samtools \
    python3-dev \
    python3-pip \
    wget && \
  rm -rf /var/lib/apt/lists/*

# Mount project directory in Docker image
ADD . /opt/hirise

# Install python packages
RUN pip3 install -r /opt/hirise/requirements.txt

# Install Meraculous
RUN \
  wget http://downloads.sourceforge.net/project/meraculous20/release-2.0.5.tgz -P /opt && \
  tar -zxvf /opt/release-2.0.5.tgz -C /opt
WORKDIR /opt/release-2.0.5
RUN ./install.sh /usr/local/

# Install HiRise scripts
RUN ln -s /opt/hirise/scripts/* /usr/local/bin/

# Download Example data
RUN \
  wget https://s3-us-west-2.amazonaws.com/dovetail-public-data1/examples.tgz -P /opt/hirise && \
  tar -zxvf /opt/hirise/examples.tgz -C /opt/hirise

# Define default command.
CMD ["bash"]
