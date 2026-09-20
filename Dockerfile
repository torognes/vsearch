FROM alpine:latest
WORKDIR /opt/vsearch
COPY . .
RUN apk add --no-cache \
        libstdc++ zlib-dev bzip2-dev \
        make g++ && \
    make clean && \
    make && \
    make install && \
    make clean && \
    apk del make g++ && \
    rm -rf /opt/vsearch
ENTRYPOINT ["/usr/local/bin/vsearch"]
