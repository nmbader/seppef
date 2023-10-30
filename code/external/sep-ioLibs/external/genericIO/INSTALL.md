# Installing genericIO


GenericIO uses submodules for most of its dependencies. There are additional optional functionality
enabled by preinstalling other packages.

GenericIO works on any modern linuxOS and normally MacOS (there is currently big in the latest version of MacOS that
should be fixed soon). 


# Pre-steps on Debian/Ubuntu:

apt-get -y update &&\
    apt-get -y  install g++ python3-numpy git make gcc libboost-all-dev  libboost-dev &&\
    apt-get -y install  cmake python3-dev python3-pytest python3-numpy-dbg libtbb-dev&& \
git clone https://github.com/pybind/pybind11.git /opt/pybind11/src && \
  mkdir /opt/pybind11/build &&\
  cd /opt/pybind11/build && \
  cmake -DCMAKE_INSTALL_PREFIX=/usr ../src  &&\
  make -j 4 install
update-alternatives --install /usr/bin/python python /usr/bin/python3 1


# Pre-steps on CentOS



yum -y install gcc gcc-c++ epel-release make git cmake  && \
yum -y install https://centos7.iuscommunity.org/ius-release.rpm && \
yum -y install tbb-devel  boost-devel boost-filesystem python36u-devel python36u-pip  python36u-pytest cmake3  && \
yum -y clean all
update-alternatives --install /usr/bin/python python /usr/bin/python3 1
/usr/bin/python3.6 -m pip install pytest
git clone https://github.com/pybind/pybind11.git /opt/pybind11/src && \
  mkdir /opt/pybind11/build &&\
  cd /opt/pybind11/build && \
  cmake3 -DCMAKE_INSTALL_PREFIX=/usr ../src  &&\
  make -j 4 install
pip3.6 install numpy




# Optional packages

1. You can build the documentation if you preinstall doxygen.
2. You can build support for old SEPlib style packages by first installing IO-libs package from
  http://cees-gitlab.stanford.EDU/SEP-external/sep-iolibs.git
3. You can also build for support for GCP by preinstalling Google GCP cloud library from
  https://github.com/googleapis/google-cloud-cpp.git.  You must follow the recommended install
  almost exactly. 
   a. The pre-installed packages must be installed standard system paths. The library will fail to build correctly if you specify alternate paths.
   b. You need to specify the -DCMAKE_CXX_FLAGS='-fPIC' (on linux) when building the library. Otherwise a relocation error will occur when linking against the library.
   c. Currently GCP can not be succesfully built in a docker (at least on a mac). It will hang when compiling.
   d. You can not include the GCP library as a sub-module.  It does not correctly handle all its internal dependencies.


# Configuring genericIO

GenericIO uses cmake for configuration.  There are many options available in the build process but for the first time users. The standard options are recommended with the possible additions of:
  - -DBUILD_GCP=yes  Build with GCP support
  - -DBUILD_DOC=no   Turn on/off documentation building

  For example a standard configuration might look like

  cmake -DCMAKE_INSTALL_PREFIX=/build/SEP -DBUILD_GCP=yes -DBUILD_DOC=no -DBUILD_TESTS=yes ../src

  where src is the location of the genericIO code
