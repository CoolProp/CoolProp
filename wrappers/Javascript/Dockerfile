## Use docker-compose to spin up this job

# Pinned to match .github/workflows/javascript_builder.yml (emcc 4.0.12).
# emcc 5.x dropped the JS-side runtime embind dispatch from -lembind; the
# documented `-sEMBIND_AOT=1` migration path fails the
# `<<< EMBIND_AOT_INVOKERS >>>` substitution against our binding shape (tsgen
# itself succeeds, but the sentinel isn't surfaced in the final JS template;
# upstream emscripten-core#24540 notes AOT method-caller coverage is
# incomplete). Bumping this image is gated on either upstream fixing AOT or
# a deeper migration; bumping it also requires bumping the matching sha in
# the CI workflow.
FROM emscripten/emsdk@sha256:744fb6a68941970951bacf9d6632041a0398260492232691ef22bbf54b0585c6

RUN apt-get -y -m update && DEBIAN_FRONTEND=noninteractive apt-get install -y dos2unix

COPY build.sh /build.sh
RUN dos2unix /build.sh

WORKDIR /
CMD bash /build.sh