From debian:9
MAINTAINER Bob Clapp <bob@sep.stanford.edu>
RUN apt-get update && apt-get install -y git gcc gfortran  make libxt-dev git
RUN  apt update && \
 apt install -y build-essential cmake git gcc g++ cmake \
        libc-ares-dev libc-ares2 libcurl4-openssl-dev libssl1.0-dev make \
        pkg-config tar wget zlib1g-dev ninja-build
RUN apt update &&\
apt install -y libtbb-dev libboost-all-dev  libboost-dev python3-dev python3-numpy-dbg python3-pytest
RUN mkdir -p /tmp/cmake/build
RUN cd /tmp/cmake &&\
   wget https://github.com/Kitware/CMake/releases/download/v3.14.3/cmake-3.14.3.tar.gz &&\
   tar xzf /tmp/cmake/cmake-3.14.3.tar.gz  &&\
 cd /tmp/cmake/build &&\
   cmake -GNinja -DMAKE_CXX_FLAGS="-O3" ../cmake-3.14.3  &&\
   cd /tmp/cmake/build &&\
	ninja install &&\
  rm -rf /tmp/cmake*
RUN apt-get -y update && \
    apt-get -y install golang  python3-numpy

ENV CFLAGS "-O3 -fPIC"
ENV CXXFLAGS "-O3 -fPIC"
RUN mkdir /opt/google && \
   cd /opt/google &&\
  wget https://github.com/googleapis/google-cloud-cpp/archive/v0.14.0.tar.gz  &&\
  tar xf v0.14.0.tar.gz &&\
   cd /opt/google/google-cloud-cpp-0.14.0/  &&\
   cmake -Hsuper -Bcmake-out -GNinja &&\
   cd /opt/google/google-cloud-cpp-0.14.0/cmake-out &&\
   ninja
  

RUN mkdir -p /opt/genericIO/build
RUN git clone http://zapad.Stanford.EDU/bob/genericIO.git /opt/genericIO/src  &&\
   cd  /opt/genericIO/build &&  \
   cmake  -DCMAKE_CXX_FLAGS="-O3  -fPIC" -DBUILD_EXAMPLES=yes   -DBUILD_GCP=yes -DPYTHON_EXECUTABLE=/usr/bin/python3.5 -DCMAKE_INSTALL_PREFIX=/opt/SEP -DCMAKE_PREFIX_PATH=/opt/google/google-cloud-cpp-0.14.0/cmake-out/external/lib/cmake  -DPYTHON_PATH=/usr/lib/python3.5 ../src &&\
   cd  /opt/genericIO/build &&  \
   make  -j 12 install 
#   rm -rf /opt/genericIO

