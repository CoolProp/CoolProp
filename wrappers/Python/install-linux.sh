#!/bin/bash
#python setup.py clean build 
python setup.py build
sudo rm -r /usr/local/lib/python2.7/dist-packages/CoolProp* 
sudo python setup.py install 
sudo chmod -R a+r /usr/local/lib/python2.7/dist-packages/CoolProp
exit 0
