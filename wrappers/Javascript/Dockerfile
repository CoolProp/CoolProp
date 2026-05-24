## Use docker-compose to spin up this job

# Pinned to match .github/workflows/javascript_builder.yml (emcc 4.0.12).
# emcc 5.x dropped the JS-side __embind_register_function / __embind_register_class
# runtime dispatch from -lembind, leaving free-function and class bindings
# inaccessible from JS (the wasm still builds and links; the resolved Module
# instance just has no F2K / PropsSI / factory / AbstractState etc.). The
# documented -sEMBIND_AOT=1 path was attempted but currently fails the link
# step ('expected to find pattern in input JS: <<< EMBIND_AOT_INVOKERS >>>')
# against our binding shape, and upstream issue emscripten-core#24540 notes
# that AOT method-caller coverage is incomplete. Lift this pin once the 5.x
# binding-exposure migration is sorted (tracked separately in beads).
FROM emscripten/emsdk@sha256:744fb6a68941970951bacf9d6632041a0398260492232691ef22bbf54b0585c6

RUN apt-get -y -m update && DEBIAN_FRONTEND=noninteractive apt-get install -y dos2unix

COPY build.sh /build.sh
RUN dos2unix /build.sh

WORKDIR /
CMD bash /build.sh