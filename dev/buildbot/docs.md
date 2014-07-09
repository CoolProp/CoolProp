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

The file ``buildbot-private.py`` (which is a python module with the passwords for the slaves as well as the buildbot website), should also be placed in the master folder next to master.cfg.  Alternatively, you can put the ``buildbot_private.py`` in another folder on the master's computer and make a soft-link in the master folder to point to the buildbot_private.py file.

If you want to completely restart the master, you can do
```
buildbot restart master
```
but usually a
```
buildbot reconfig master
```
will do the job since it will just reparse the configuration file without signing you out of the server

To ensure that the buildbot server stays online, you can make a script with the contents
```
buildbot start /path/to/master_folder
```
and add it to a cron job

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
