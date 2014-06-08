Set up linux slave
==================

In some folder (probably the home folder), do::

    sudo apt-get buildbot buildbot-slave
    buildslave create-slave slave 10.0.2.2:9989 example-slave pass
    buildslave start slave

``slave`` is the folder that the configuration will go into, ``example-slave`` is the name of the slave, ``pass`` is the password, ``10.0.0.2`` is the IP address of the host that the VM uses (shouldn't need to change it, see OSX host instructions below)

This should start the Twisted server, yielding output like::

    ian@ian-VirtualBox:~$ buildslave start slave/
    Following twistd.log until startup finished..
    2014-06-08 20:17:50+0200 [-] Log opened.
    2014-06-08 20:17:50+0200 [-] twistd 14.0.0 (/usr/bin/python 2.7.5) starting up.
    2014-06-08 20:17:50+0200 [-] reactor class: twisted.internet.epollreactor.EPollReactor.
    2014-06-08 20:17:50+0200 [-] Starting BuildSlave -- version: 0.8.8
    2014-06-08 20:17:50+0200 [-] recording hostname in twistd.hostname
    2014-06-08 20:17:50+0200 [-] Starting factory <buildslave.bot.BotFactory instance at 0x8a70a4c>
    2014-06-08 20:17:50+0200 [-] Connecting to 10.0.2.2:9989
    2014-06-08 20:17:50+0200 [Broker,client] message from master: attached
    The buildslave appears to have (re)started correctly.`

Configuration of OSX Virtualbox host
====================================

Thanks to http://stackoverflow.com/questions/1261975/addressing-localhost-from-a-virtualbox-virtual-machine and the comment of spsaucier you basically need to do the following, copied verbatim:

To enable this on OSX I had to do the following:

1. Shut your virtual machine down.
2. Go to ``VirtualBox Preferences -> Network -> Host-only Networks ->`` click the "+" icon. Click OK.
3.Select your box and click the "Settings" icon -> Network -> Adapter 2 -> On the "Attached to:" dropdown, select "Host-only Adapter" and your network (vboxnet0) should show up below by default. Click OK.
4. Once you start your box up again, you should be able to access localhost at http://10.0.2.2/

You can refer to it by localhost and access other localhosted sites by adding their references to the hosts file (C:\windows\system32\drivers\etc\hosts) like the following::

	10.0.2.2    localhost
	10.0.2.2    subdomain.localhost
