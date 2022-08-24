## Use docker-compose to spin up this job

FROM emscripten/emsdk

RUN apt-get -y -m update && DEBIAN_FRONTEND=noninteractive apt-get install -y dos2unix

COPY build.sh /build.sh
RUN dos2unix /build.sh

WORKDIR /
CMD bash /build.sh