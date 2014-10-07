
********
Buildbot
********

Buildbot masters and slaves
===========================

Master
------

From the root of the git checkout (this will use the master.cfg from CoolProp)::

    pip install buildbot
    cd dev/buildbot
    buildbot create-master master
    buildbot start master

The file ``buildbot-private.py`` (which is a python module with the passwords for the slaves as well as the buildbot website), should also be placed in the master folder next to master.cfg.  Alternatively, you can put the ``buildbot_private.py`` in another folder on the master's computer and make a soft-link in the master folder to point to the buildbot_private.py file.

If you want to completely restart the master, you can do::

    buildbot restart master

but usually a::

    buildbot reconfig master

will do the job since it will just reparse the configuration file without signing you out of the server

To ensure that the buildbot server stays online, you can make a script with the contents::

    buildbot start /path/to/master_folder

and add it to a cron job

Slaves
------

To start a slave connected to a buildbot master at IP address 10.0.0.2 (default for 
host for VirtualBox), with a slave named ``example-slave`` and passsword ``pass``, 
run the command::

    buildslave create-slave slave 10.0.0.2:9989 example-slave pass
    buildslave start slave


If the master is somewhere else, just change the IP address.  As of Sept, 2014, the 
master was at www.coolprop.dreamhosters.com.  The buildbot_private.py on the master 
holds the required passwords.

On linux, you can add the following lines to the end of your ``~/.profile`` file (similar 
ideas apply on other platforms) to autostart the slave when the user logs in::

    # Connect to the buildbot master
    buildslave start ~/slave
    
... or even better, you install a service that gets started and shutdown together with 
your computer. For Debian/Ubuntu, we recommend as script like::


    #! /bin/sh
    # /etc/init.d/buildbotslave
    #
    
    # Some things that run always
    touch /var/lock/buildbotslave
    
    # Carry out specific functions when asked to by the system
    case "$1" in
      start)
        echo "Starting script buildbotslave "
        sudo -u username buildslave start "/home/username/slave/"
        ;;
      stop)
        echo "Stopping script buildbotslave"
        sudo -u username buildslave stop "/home/username/slave/"
        ;;
      restart)
        echo "Restarting script buildbotslave"
        sudo -u username buildslave stop "/home/username/slave/"
        sudo -u username buildslave start "/home/username/slave/"
        ;;
      *)
        echo "Usage: /etc/init.d/buildbotslave {start|stop|restart}"
        exit 1
        ;;
    esac
    
    exit 0
    
Which the can be added to the scheduler with ``update-rc.d buildbotslave defaults``. 
This should gracefully terminate the bot at shutdown and restart it again after reboot. 
To disable the service, run ``update-rc.d -f buildbotslave remove``.


    


Setting MIME type handler
=========================

To change the MIME types on the server so that unknown file types will map properly to ``application/octet-stream``, modify the ``buildbot.tac`` file to add the following block::

  from twisted.web.static import File

  webdir = File("public_html")
  webdir.contentTypes['.mexw32'] = 'application/octet-stream'
  webdir.contentTypes['.mexw64'] = 'application/octet-stream'
  webdir.contentTypes['.mexmaci64'] = 'application/octet-stream'
  webdir.contentTypes['.jnilib'] = 'application/octet-stream'
  webdir.contentTypes['.mexa64'] = 'application/octet-stream'
  webdir.contentTypes['.oct'] = 'application/octet-stream'
  webdir.contentTypes['.whl'] = 'application/octet-stream'
  webdir.contentTypes['.dylib'] = 'application/octet-stream'
  ...

and then do a ``buildbot restart master``


Nightly Documentation Builds
============================

Some parts of the documentation are quite involved. That is why we decided not
to rebuild the whole documentation after every commit. There is a special buildbot
slave that runs once a day and performs the most expensive jobs. This covers the
generation of validation figures for all fluids and the fitting reports for the
incompressible fluids.

If you have some tasks that take a long time, make sure to add them to that
special machine. This helps us to keep the continuous integration servers running
with an acceptable latency with regard to the commit to the git repository.
