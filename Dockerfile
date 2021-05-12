FROM ubuntu:18.04 AS crass_build

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils build-essential wget autotools-dev autoconf automake libtool libz-dev

RUN wget --no-check-certificate https://github.com/apache/xerces-c/archive/refs/tags/v3.1.1.tar.gz && tar -xzf v3.1.1.tar.gz && cd xerces-c-3.1.1 && ./reconf && ./configure && make && make install && cd ..

COPY . /

RUN ./autogen.sh && ./configure && make && make install

FROM ubuntu:18.04

COPY --from=crass_build /usr/local/lib/libxerces-c-3.1.so /usr/lib/
COPY --from=crass_build /usr/local/bin/crass /usr/bin/
COPY --from=crass_build /usr/local/bin/crisprtools /usr/bin/
