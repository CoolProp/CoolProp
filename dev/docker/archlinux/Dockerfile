## Use docker compose up --build to spin up this job

FROM archlinux:latest

RUN pacman -Sy && pacman -S --noconfirm cmake git gcc python-setuptools python-packaging cython make

RUN git clone --recursive https://github.com/CoolProp/CoolProp && \
	cd CoolProp/wrappers/Python && \
	# python setup.py build --verbose install # should work, but python on arch does not seem to pass include paths to the linker, which is needed for incbin inclusion of z-compressed fluid data
	python setup.py build --verbose install cmake=default,any