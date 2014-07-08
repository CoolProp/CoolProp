Buildbot masters and slaves
===========================

Master
------

From the root of the git checkout (this will use the master.cfg from CoolProp)
```
pip install buildbot
cd dev/buildbot
buildbot create-master master
buildbot start master
```

The file ``buildbot-private.py`` (which is a python module with the passwords for the slaves as well as the buildbot website), should also be placed in the master folder next to master.cfg.

If you want to completely restart the master, you can do
```
buildbot restart master
```
but usually a
```
buildbot reconfig master
```
will do the job since it will just reparse the configuration file without signing you out of the server

You can add the following lines to the end of your ``.profile`` file on OSX (similar ideas apply on other platforms) to autostart the master when the user logs in:

```
# Startup the buildbot master
buildbot start ~/master
```

Slaves
------

To start a slave connected to a buildbot master at IP address 10.0.0.2 (default for host for VirtualBox), with a slave named ``example-slave`` and passsword ``pass``, run the command

```
buildslave create-slave slave 10.0.0.2:9989 example-slave pass
buildslave start slave
```

If the master is somewhere else, just change the IP address.  

On linux, you can add the following lines to the end of your ``.xsessionrc`` file (similar ideas apply on other platforms) to autostart the slave when the user logs in:

```
# Connect to the buildbot master
buildslave start ~/slave
```
