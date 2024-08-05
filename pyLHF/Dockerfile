FROM quay.io/pypa/manylinux_2_28_x86_64

#Install LHF dependency
RUN yum install -y eigen3-devel gmp-devel mpfr-devel openmpi-devel && yum clean all

# Download and install Boost
RUN curl -L https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.bz2 | tar --bzip2 -xf - \
    && cd /boost_1_85_0/ \
    && ./bootstrap.sh --prefix=/usr/local && ./b2 install -d0 \
    && cd / && rm -rf /boost_1_85_0

ENV PATH="$PATH:/usr/lib64/openmpi/bin"

COPY script.sh /script.sh

ENTRYPOINT ["/bin/bash", "-c", "/script.sh && exec /bin/bash"]