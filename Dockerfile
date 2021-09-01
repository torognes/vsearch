FROM alpine:latest
WORKDIR /opt/vsearch
COPY . .
RUN apk add --no-cache \
        libstdc++ zlib-dev bzip2-dev \
        autoconf automake make g++ && \
    ./autogen.sh && \
    ./configure CFLAGS="-O3" CXXFLAGS="-O3" && \
    make clean && \
    make && \
    make install && \
    make clean && \
    apk del autoconf automake make g++ && \
    rm -rf /opt/vsearch
ENTRYPOINT ["/usr/local/bin/vsearch"]
