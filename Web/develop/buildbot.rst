
********
Buildbot
********

Buildbot masters and slaves
===========================

Master
------

From the root of the git checkout (this will use the master.cfg from CoolProp)::

    virtualenv buildbot-sandbox
    source buildbot-sandbox/bin/activate
    pip install sqlalchemy==0.7.10 buildbot
    cd dev/buildbot
    buildbot create-master master
    buildbot start master

The file ``buildbot-private.py`` (which is a python module with the passwords for the slaves as well as
the buildbot website), should also be placed in the master folder next to master.cfg.  Alternatively,
you can put the ``buildbot_private.py`` in another folder on the master's computer and make a soft-link
in the master folder to point to the buildbot_private.py file.

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
host for VirtualBox), with a slave named ``a-slave`` and passsword ``pass``,
run the command::

    virtualenv a-slave-sandbox
    source a-slave-sandbox/bin/activate
    pip install sqlalchemy==0.7.10 buildbot-slave
    buildslave create-slave a-slave coolprop.dreamhosters.com:port a-slave pass
    buildslave start a-slave

If the master is somewhere else, just change the IP address.  As of Sept, 2014, the
master was at www.coolprop.dreamhosters.com.  The buildbot_private.py on the master
holds the required passwords.


Python Slaves
-------------

Based on the miniconda Python ecosystem, you can create your own virtual
environments for building the Python wheels. This requires the following
steps on a Windows machine::

  conda create -n CoolProp27 python=2
  conda create -n CoolProp34 python=3
  conda install -n CoolProp27 cython pip pywin32
  conda install -n CoolProp34 cython pip pywin32

  activate CoolProp27
  pip install wheel
  deactivate
  activate CoolProp34
  pip install wheel
  deactivate

Please repeat the steps above for both 32bit and 64bit Python environments. 

On a Linux system, things only change a little bit::

  conda create -n CoolProp27 python=2
  conda create -n CoolProp34 python=3
  conda install -n CoolProp27 cython pip
  conda install -n CoolProp34 cython pip

  source activate CoolProp27
  pip install wheel
  deactivate
  source activate CoolProp34
  pip install wheel
  deactivate

Please make sure that the standard shell ``/bin/sh`` used by the builbot is
bash or zsh. We make use of the ``source`` command, which is not part of the
POSIX specification.

At the moment, it is not possible to use several slaves for the same build job.
We have to find a new way to generate the configuration.

Information on building the single wrappers can be found on 
:ref:`this dedicated page<wrapper_common_prereqs>`.


Buildbot as a service (Windows)
-------------------------------

On Windows, you create a batch script that activates your virtual environment
and starts the buildslave::

  @echo off
  call "C:\Program Files (x86)\Miniconda32_27\Scripts\activate.bat" Buildbot
  buildslave start "C:\CoolProp-slave"

This script can then be added to the system services via::

  sc create <serviceName> binpath= <pathToBatFile> DisplayName= "CoolProp Buildbot" start= auto

You might want to run ``services.msc`` to edit the user that runs the service. If
you are tired of the error messages from the non-returning script, you could
also use a service wrapper like `NSSM <http://nssm.cc/>`_ to start the script.


Buildbot as a daemon (Linux)
----------------------------

On linux, you can add the following lines to the end of your ``~/.profile`` file (similar
ideas apply on other platforms) to autostart the slave when the user logs in::

    # Connect to the buildbot master
    buildslave start ~/slave

... or even better, you install a service that gets started and shutdown together with
your computer. For Debian/Ubuntu, we recommend as script like::

    #! /bin/sh
    ### BEGIN INIT INFO
    # Provides:          buildslave
    # Required-Start:    $remote_fs $syslog
    # Required-Stop:     $remote_fs $syslog
    # Default-Start:     2 3 4 5
    # Default-Stop:      0 1 6
    # Short-Description: A script to start the buildbot slave at boot time
    # Description:       This file activates the virtual environment and starts
    #                    the buildbot slaves. It also shuts them down if the
    #                    system is halted. Place it in /etc/init.d.
    ### END INIT INFO

    # Author: Jorrit Wronski <jowr@mek.dtu.dk>
    #
    # Please remove the "Author" lines above and replace them
    # with your own name if you copy and modify this script.

    EXECUSER=username
    NAME="a-slave"
    CTRLSCRI="/home/username/$NAME.bsh"

    # Load the VERBOSE setting and other rcS variables
    . /lib/init/vars.sh

    # Define LSB log_* functions.
    # Depend on lsb-base (>= 3.2-14) to ensure that this file is present
    # and status_of_proc is working.
    . /lib/lsb/init-functions

    #
    # Function that starts the daemon/service
    #
    do_start(){
      sudo -u $EXECUSER $CTRLSCRI start
      #start-stop-daemon --start --user $EXECUSER --chuid $EXECUSER --startas $CTRLSCRI -- start
      RETVAL="$?"
      return "$RETVAL"
    }

    #
    # Function that stops the daemon/service
    #

    # Function that stops the daemon/service
    #
    do_stop() {
      #start-stop-daemon --stop --user $EXECUSER --startas
      sudo -u $EXECUSER $CTRLSCRI stop
      RETVAL="$?"
      return "$RETVAL"
    }

    case "$1" in
    start)
        log_action_msg "Starting $NAME"
        do_start
        ;;
    stop)
        log_action_msg "Stopping $NAME"
        do_stop
        ;;
    restart)
        log_action_msg "Restarting $NAME"
        do_stop
        do_start
        ;;
    *)
        log_action_msg "Usage: $0 {start|stop|restart}"
        exit 2
        ;;
    esac
    exit 0

Which the can be added to the scheduler with ``update-rc.d buildslave defaults``.
This should gracefully terminate the bot at shutdown and restart it again after reboot.
To disable the service, run ``update-rc.d -f buildslave remove``. You can enable and
disable the daemon by runnning ``update-rc.d buildslave enable|disable``. Note that the
example above call a user-script that activates the virtual environment and starts
the buildslave. Such a script could look like this::

    #! /bin/bash
    #
    # Description:       This file activates the virtual environment and starts
    #                    the buildbot slaves. It also shuts them down if the
    #                    system is halted.
    #
    # Author: Jorrit Wronski <jowr@mek.dtu.dk>
    #
    # Please remove the "Author" lines above and replace them
    # with your own name if you copy and modify this script.
    #
    VIRTENV=/home/username/a-slave-sandbox
    SLAVEDIR=/home/username/a-slave
    #
    # Carry out specific functions when asked to by the system
    case "$1" in
      start)
        echo "Starting script buildbotslave "
        source $VIRTENV/bin/activate
        $VIRTENV/bin/buildslave start $SLAVEDIR
        ;;
      stop)
        echo "Stopping script buildbotslave"
        $VIRTENV/bin/buildslave stop $SLAVEDIR
        ;;
      restart)
        echo "Restarting script buildbotslave"
        source $VIRTENV/bin/activate
        $VIRTENV/bin/buildslave stop $SLAVEDIR
        $VIRTENV/bin/buildslave start $SLAVEDIR
        ;;
      *)
        echo "Usage: $0 {start|stop|restart}"
        exit 1
        ;;
    esac
    exit 0




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




Documentation Builds
====================

Some parts of the documentation are quite involved. That is why we decided not
to rebuild the whole documentation after every commit. There is a special python
script that runs a day and performs the most expensive jobs during
documentation rebuild. This covers the generation of validation figures for all
fluids and the fitting reports for the incompressible fluids.

If you have some tasks that take a long time, make sure to add them to that
special script in ``Web/scripts/__init__.py``. This helps us to keep the continuous
integration servers running with an acceptable latency with regard to the commits
to the git repository. However, if you are unlucky and your commit coincides with
figure generation, you will experience a long
delay between your commit and the appearance of the freshly generated documentation
on the website. You can follow the progress in the logfiles on the buildbot master though.
