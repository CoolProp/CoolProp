#!/bin/bash
#python setup.py clean build 
python setup.py build
sudo rm -r /Library/Python/2.7/site-packages/CoolProp* 
sudo python setup.py install 
sudo chmod -R a+r /Library/Python/2.7/site-packages/CoolProp
exit 0
