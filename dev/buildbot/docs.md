Using and configuring the buildbot masters and slaves
=====================================================

Master
------

```
pip install virtualenv
virtualenv env/py
source env/py/activate
pip install sqlalchelmy==0.7.10 buildbot
buildbot create-master master
mv master/master.cfg.sample master/master.cfg
buildbot start master
```

Slaves
------

To start a slave connected to a buildbot master at IP address 10.0.0.2 (default for host for VirtualBox), with a slave named ``example-slave` and passsword ``pass``, run the command

```
buildslave create-slave slave 10.0.0.2:9989 example-slave pass
buildslave start slave
```

If the master is somewhere else, just change the IP address.  
