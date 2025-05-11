## Use docker compose up --build to spin up this job

FROM archlinux:latest

RUN pacman -Sy && pacman -S --noconfirm cmake git gcc python-setuptools python-packaging cython

RUN git clone --recursive https://github.com/CoolProp/CoolProp && \
	cd CoolProp/wrappers/Python && \
	python setup.py build --verbose install