## Use docker compose up --build to spin up this job
## which will generate a .tar.xz file in the working directory

FROM ubuntu:24.04

RUN apt-get -y -m update && DEBIAN_FRONTEND=noninteractive apt-get install -y cmake g++ git zip wget xz-utils

RUN mkdir /boost && \
	wget -c --no-check-certificate https://boostorg.jfrog.io/artifactory/main/release/1.87.0/source/boost_1_87_0_rc1.tar.gz -O - | tar -xz -C /boost && \
	cd /boost/boost_1_87_0/ && \
	./bootstrap.sh && \
	./b2 tools/bcp

WORKDIR /boost/boost_1_87_0/
RUN mkdir /boost_CoolProp && \
	./bin.v2/tools/bcp/gcc-13/debug/x86_64/link-static/threading-multi/visibility-hidden/bcp predef/other/endian.h boost/fusion/sequence/intrinsic/size.hpp boost/fusion/algorithm/iteration/for_each.hpp boost/fusion/include/mpl.hpp boost/fusion/sequence/intrinsic/at.hpp boost/utility/string_ref.hpp boost/utility/string_view.hpp boost/mpl/size.hpp boost/variant.hpp boost/assert.hpp boost/preprocessor.hpp boost/fusion/support/is_sequence.hpp boost/optional.hpp boost/operators.hpp boost/version.hpp /boost_CoolProp && \
	zip -r /boost_CoolProp.zip /boost_CoolProp &&  \
	tar cJf /boost_CoolProp.tar.xz /boost_CoolProp	

CMD cp /*.tar.xz /output
