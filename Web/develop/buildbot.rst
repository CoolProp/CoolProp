
********
Buildbot
********

Buildbot masters and workers
============================

Master
------

From the root of the git checkout (this will use the master.cfg from CoolProp)::

    virtualenv buildbot-sandbox
    source buildbot-sandbox/bin/activate
    pip install sqlalchemy==0.7.10 buildbot
    cd dev/buildbot
    buildbot create-master master
    buildbot start master

The file ``buildbot-private.py`` (which is a python module with the passwords for the workers as well as
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

and add it to a cron job.


The work with the ``master.cfg`` has proven to be a little more complicated. Storing changes in the repository 
and pulling them to the repository on the server is a little too cumbersome, especially when many iterations 
are needed to fix issues with the configurations. We already discussed this a couple of `times <https://github.com/CoolProp/CoolProp/issues/1052>`, 
but here the latest version of the preferred work flow for editing the ``master.cfg`` file: 

1. SSH to server and play with configuration file.
2. Run `./scripts/buildbot.sh reconfig` on server to activate changes, **ignore** git warnings.
3. Repeat 1. and 2. until you are happy with you configuration 
4. Copy the server file to your local repository.
5. Commit and push local repository.
6. Run `pushd buildbot/CoolProp.git/ && git reset --hard origin/master && popd` on server.
7. Run `./scripts/buildbot.sh reconfig` on server, make sure there are **no** warnings from git.



Workers
-------

To start a worker connected to a buildbot master at IP address 10.0.0.2 (default for host for VirtualBox), with a worker named ``a-worker`` and passsword ``pass`` in the directory ``workdir``,
run the command::

    virtualenv a-worker-sandbox
    source a-worker-sandbox/bin/activate
    pip install sqlalchemy==0.7.10 buildbot-worker
    buildslave create-slave workdir coolprop.dreamhosters.com:port a-worker pass
    buildslave restart workdir

If the master is somewhere else, just change the IP address.  As of Sept, 2014, the
master was at www.coolprop.dreamhosters.com, and the port is specified in the master.cfg ``file``.  The file ``buildbot_private.py`` on the master holds the required passwords.

OSX Virtualbox host
-------------------

Thanks to http://stackoverflow.com/questions/1261975/addressing-localhost-from-a-virtualbox-virtual-machine and the comment of spsaucier you basically need to do the following, copied verbatim:

To enable this on OSX I had to do the following:

1. Shut your virtual machine down.
2. Go to ``VirtualBox Preferences -> Network -> Host-only Networks ->`` click the "+" icon. Click OK.
3. Select your box and click the "Settings" icon -> Network -> Adapter 2 -> On the "Attached to:" dropdown, select "Host-only Adapter" and your network (vboxnet0) should show up below by default. Click OK.
4. Once you start your box up again, you should be able to access localhost at http://10.0.2.2/

You can refer to it by localhost and access other localhosted sites by adding their references to the hosts file (C:\windows\system32\drivers\etc\hosts) like the following::

	10.0.2.2    localhost
	10.0.2.2    subdomain.localhost
    

Python workers
--------------

Based on the miniconda Python ecosystem, you can create your own virtual
environments for building the Python wheels. This requires the following
steps on a Windows machine::

    conda create -y -n CoolProp27 python=2.7 cython pip pywin32 requests jinja2 pyyaml pycrypto wheel ndg-httpsclient
    conda create -y -n CoolProp36 python=3.6 cython pip pywin32 requests jinja2 pyyaml pycrypto wheel 
    conda create -y -n CoolProp37 python=3.7 cython pip pywin32 requests jinja2 pyyaml pycrypto wheel 
    conda create -y -n CoolProp38 python=3.8 cython pip pywin32 requests jinja2 pyyaml pycrypto wheel 
    #conda create -y -n CoolPropWorker python=2.7 pip && conda activate CoolPropWorker && pip install buildbot-worker && conda deactivate
    conda create -y -n CoolPropWorker python=2.7 pip && conda activate CoolPropWorker && pip install buildbot-slave==0.8.14 && conda deactivate

Please repeat the steps above for **both 32bit and 64bit** Python environments.

On a Linux system, things only change a little bit::

    conda create -n CoolProp27 python=2.7 cython pip requests jinja2 pyyaml pycrypto wheel
    conda create -n CoolProp36 python=3.6 cython pip requests jinja2 pyyaml pycrypto wheel
    conda create -n CoolProp37 python=3.7 cython pip requests jinja2 pyyaml pycrypto wheel
    conda create -n CoolProp38 python=3.8 cython pip requests jinja2 pyyaml pycrypto wheel
    conda create -n CoolPropWorker pip && source activate CoolPropWorker && pip install buildbot-worker && source deactivate

Please make sure that the standard shell ``/bin/sh`` used by the builbot is
bash or zsh. We make use of the ``source`` command, which is not part of the
POSIX specification. In Debian, ``dpkg-reconfigure dash`` can be used.

At the moment, it is not possible to use several workers for the same build job.
We have to find a new way to generate the configuration.

Information on building the single wrappers can be found on
:ref:`this dedicated page<wrapper_common_prereqs>`.

For uploading generated binary python files to PYPI, you should create a file ``~\.pypirc`` with the contents::

	[distutils]
	index-servers=
	    pypi
	    test

	[test]
	repository = https://testpypi.python.org/pypi
	username = user
	password = XXXXXXXXXXXXXXXX

	[pypi]
	repository = https://pypi.python.org/pypi
	username = user
	password = XXXXXXXXXXXXXXXX

Buildbot as a service (Windows)
-------------------------------

On Windows, you create a batch script that activates your virtual environment
and starts the buildbot worker::

    @echo off
    call "C:\Program Files (x86)\Miniconda32_27\Scripts\activate.bat" Buildbot
    buildbot-worker start "C:\CoolProp-worker"

This script can then be added to the system services via::

    sc create <serviceName> binpath= <pathToBatFile> DisplayName= "CoolProp Buildbot" start= auto

You might want to run ``services.msc`` to edit the user that runs the service. If
you are tired of the error messages from the non-returning script, you could
also use a service wrapper like `NSSM <http://nssm.cc/>`_ to start the script.

Buildbot and launchd (Mac OS)
-----------------------------
As written in the `Buildbot Wiki <http://trac.buildbot.net/wiki/UsingLaunchd>`_,
you can start your workers automatically with a so called ``plist`` or property list.
Place the example content below in a file called ``/Library/LaunchDaemons/org.coolprop.a-worker.plist``
and make sure it is owned by the user ``root`` and the group ``wheel``::

    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
    <plist version="1.0">
    <dict>
        <key>StandardOutPath</key>
        <string>org.coolprop.a-worker.log</string>
        <key>StandardErrorPath</key>
        <string>org.coolprop.a-worker-err.log</string>
        <key>Label</key>
        <string>org.coolprop.a-worker</string>
        <key>Program</key>
        <string>/Users/buildbot/bin/a-worker.command</string>
        <key>RunAtLoad</key>
        <true/>
        <key>KeepAlive</key>
        <dict>
            <key>SuccessfulExit</key>
            <false/>
        </dict>
        <key>GroupName</key>
        <string>staff</string>
        <key>UserName</key>
        <string>buildbot</string>
        <key>WorkingDirectory</key>
        <string>/Users/buildbot/worker/logs</string>
        <key>SessionCreate</key>
        <true/>
    </dict>
    </plist>

Please change the file above according to your needs and pay special attention
to username and path definitions. The script ``a-worker.command`` that is called
by ``launchd`` could look like this one::

    #!/bin/bash
    #
    # Description: This file call the control script to start and
    #              stop the buildbot worker. It stays open when being
    #              called and waits for a signal to terminate running
    #              and endless while-loop. After catching a signal
    #              to terminate, it shuts down the build worker and
    #              returns. It is a wrapper for another Bash script
    #              allowing us to use launchd on MacOS.
    #
    # Author: Jorrit Wronski <jowr@mek.dtu.dk>
    #
    # Please remove the "Author" lines above and replace them
    # with your own name if you copy and modify this script.
    #
    # If you experience any problems with the PATH variable on OSX,
    # this setting might be for you:
    if [ -x /usr/libexec/path_helper ]; then
      eval `/usr/libexec/path_helper -s`
    fi
    #
    CTRLSCRI="/Users/username/a-worker.bsh"
    #
    trap "$CTRLSCRI stop; exit 0; " TERM SIGINT SIGTERM
    #
    $CTRLSCRI start & wait
    # Just idle for one hour and keep the process alive
    # waiting for SIGTERM.
    while : ; do
      sleep 3600 & wait
    done
    #
    echo "The endless loop terminated, something is wrong here."
    exit 1

Note that this script calls another Bash script that does the actual work. We hope
to simplify maintenance by using a common control script for Linux and MacOS as
shown in :ref:`workerscript`.

Or alternatively, you can just launch buildbot worker directly if you do not use conda environment::

    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
    <plist version="1.0">
    <dict>
        <key>KeepAlive</key>
        <true/>
        <key>Label</key>
        <string>com.start.buildbot</string>
        <key>ProgramArguments</key>
        <array>
            <string>/Users/Ian/anaconda/bin/buildworker</string>
            <string>restart</string>
            <string>worker</string>
        </array>
        <key>RunAtLoad</key>
        <true/>
        <key>StandardErrorPath</key>
        <string>/Users/Ian/.buildbot_stderr</string>
        <key>StandardOutPath</key>
        <string>/Users/Ian/.buildbot_stdout</string>
        <key>UserName</key>
        <string>Ian</string>
        <key>WorkingDirectory</key>
        <string>/Users/Ian</string>
    </dict>
    </plist>

Buildbot as a daemon (Linux)
----------------------------

On Linux, you can add the following lines to the end of your ``~/.profile`` file (similar
ideas apply on other platforms) to start the worker automatically at user log in::

    # Connect to the buildbot master
    buildworker start ~/worker

... or even better, you install a service that gets started and shutdown together with
your computer. For Debian/Ubuntu, we recommend a script like::

    #! /bin/sh
    ### BEGIN INIT INFO
    # Provides:          buildworker
    # Required-Start:    $remote_fs $syslog
    # Required-Stop:     $remote_fs $syslog
    # Default-Start:     2 3 4 5
    # Default-Stop:      0 1 6
    # Short-Description: A script to start the buildbot worker at boot time
    # Description:       This file activates the virtual environment and starts
    #                    the buildbot workers. It also shuts them down if the
    #                    system is halted. Place it in /etc/init.d.
    ### END INIT INFO

    # Author: Jorrit Wronski <jowr@ipu.dk>
    #
    # Please remove the "Author" lines above and replace them
    # with your own name if you copy and modify this script.

    EXECUSER=username
    NAME="a-worker"
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

Which then can be added to the scheduler with ``update-rc.d buildworker defaults``.
This should gracefully terminate the bot at shutdown and restart it again after reboot.
To disable the service, run ``update-rc.d -f buildworker remove``. You can enable and
disable the daemon by running ``update-rc.d buildworker enable|disable``.

If you run a distribution that uses systemd, like CentOS, you might find the
following unit file helpful, which can be placed in ``/etc/systemd/system/coolpropworker.service``
or in ``~/.config/systemd/user/coolpropworker.service``::

    [Unit]
    Description=CoolProp Linux buildbot
    
    [Service]
    User=buildbot
    Type=forking
    WorkingDirectory=/home/buildbot
    StandardOutput=syslog
    StandardError=syslog
    SyslogIdentifier=CoolPropBuilder
    ExecStartPre=/bin/bash --login -c 'env > /tmp/buildbot-environment-file'
    EnvironmentFile=-/tmp/buildbot-environment-file
    ExecStart=/home/buildbot/buildbot.bsh start
    ExecStop=/home/buildbot/buildbot.bsh stop
    ExecReload=/home/buildbot/buildbot.bsh restart
    Restart=on-abnormal
    
    [Install]
    WantedBy=multi-user.target

Install the service with ``sudo systemctl enable coolpropworker.service`` and
activate it using ``sudo systemctl start coolpropworker.service``.


.. _workerscript:

Buildbot worker management (Mac OS and Linux)
---------------------------------------------

Note that the two examples above call a user-script to activate the virtual
environment and start the buildbot worker. Such a script could look like this::

    #!/bin/bash
    #
    # Description: This file activates the virtual environment and starts
    #              the buildbot workers. It is also used to shut them down
    #              during system shutdown.
    #
    # Author: Jorrit Wronski <jowr@ipu.dk>
    #
    # Please remove the "Author" lines above and replace them
    # with your own name if you copy and modify this script.
    #
    VIRTENV="a-worker-sandbox"
    WORKERDIR="/home/username/a-worker"
    #
    ## For virtualenv
    #ACTICM="source $VIRTENV/bin/activate"
    ##DEACCM="source $VIRTENV/bin/deactivate"
    #
    # For miniconda
    MINICO="/home/username/miniconda/bin/activate"
    ACTICM="source $MINICO $VIRTENV"
    #DEACCM="source deactivate"
    #
    # Carry out specific functions when asked to by the system
    case "$1" in
      create)
        echo "Creating buildbot worker"
        buildbot-worker create-worker $WORKERDIR coolprop.dreamhosters.com:port a-worker pass
        #$DEACCM
      start)
        echo "Starting buildbot worker"
        $ACTICM
        buildbot-worker start $WORKERDIR
        #$DEACCM
        ;;
      stop)
        echo "Stopping buildbot worker"
        $ACTICM
        buildbot-worker stop $WORKERDIR
        #$DEACCM
        ;;
      restart)
        echo "Restarting buildbot worker"
        $ACTICM
        buildbot-worker restart $WORKERDIR
        #$DEACCM
        ;;
      *)
        echo "Usage: $0 {create|start|stop|restart}"
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


Starting VirtualBox images at boot
==================================

You can use the built-in functionality https://www.virtualbox.org/manual/ch09.html#autostart on Linux and Mac or use
your own configuration and create a daemon entry in Library/LaunchDaemons.  Make sure you use full paths to VBoxManage::

    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
    <plist version="1.0">
    <dict>
        <key>GroupName</key>
        <string>staff</string>
        <key>InitGroups</key>
        <true/>
        <key>KeepAlive</key>
        <false/>
        <key>Label</key>
        <string>com.start.windows.vm</string>
        <key>ProgramArguments</key>
        <array>
            <string>/usr/bin/Vboxmanage</string>
            <string>startvm</string>
            <string>xp</string>
        </array>
        <key>RunAtLoad</key>
        <true/>
        <key>StandardErrorPath</key>
        <string>/Users/Ian/.virtualbox_window_stderr</string>
        <key>StandardOutPath</key>
        <string>/Users/Ian/.virtualbox_windows_stdout</string>
        <key>UserName</key>
        <string>Ian</string>
    </dict>
    </plist>


Documentation Builds
====================

Some parts of the documentation are quite involved. That is why we decided not
to rebuild the whole documentation after every commit. There is a special python
script that runs once a day and performs the most expensive jobs during
documentation rebuild. This covers the generation of validation figures for all
fluids and the fitting reports for the incompressible fluids.

If you have some tasks that take a long time, make sure to add them to that
special script in ``Web/scripts/__init__.py``. This helps us to keep the continuous
integration servers running with an acceptable latency with regard to the commits
to the git repository. However, if you are unlucky and your commit coincides with
figure generation, you will experience a long
delay between your commit and the appearance of the freshly generated documentation
on the website. You can follow the progress in the logfiles on the buildbot master though.


Work in Progress - Dockerfile Generator
=======================================

To make it short, here is what you need to know if you trust us and the docker 
build system: 

* Make sure to set the correct environment variables in an additional file before 
  you run a container, call it for example ``Dockerfile.worker.env.list``::

    WORKERDIR=/home/buildbot/workerdir
    BUILDMASTER=bots.coolprop.org
    BUILDMASTER_PORT=port
    WORKERNAME=workername
    WORKERPASS=pass
    WORKER_ENVIRONMENT_BLACKLIST=notused
    BOTADMIN=Author Name
    BOTEMAIL=noreply@coolprop.org
    BOTHOST=A short description of the host computer

* You can then run the official coolprop buildbot configuration with::

    docker run -d --env-file ./Dockerfile64.worker.env.list --name=CoolProp64-worker coolprop/workerpython 
    docker run -d --env-file ./Dockerfile32.worker.env.list --name=CoolProp32-worker coolprop/workerpython32
    
  The above commands launch background processes using the docker contains for the Python buildworkers in 
  64bit and 32bit, respectively. 

* Some steps require the upload of files to different servers. In such cases, you 
  should copy your SSH configuration or other login information to the container to 
  make use of the automatic login that is required for rsync to work properly::

    docker cp ${HOME}/.ssh ${WORKERNAME}:/home/buildbot/
    docker cp ${HOME}/.pypirc ${WORKERNAME}:/home/buildbot/
    docker exec --user root ${WORKERNAME} chown -R buildbot /home/buildbot/.ssh /home/buildbot/.pypirc
	docker exec --user root ${WORKERNAME} chgrp -R buildbot /home/buildbot/.ssh /home/buildbot/.pypirc

.. note::
  If you cannot copy the SSH keys, you can change the upload function in the 
  master configuration to employ the built-in upload framework of buildbot. 

Why the containers? In 2015, some of the buildbot workers did not perform as expected. 
Especially the Python builds on the 64bit Linux machine took ages to complete and we 
could not find any obvious reason for this behaviour. 

To make sure that there are no hidden flaws in the configuration of the buildbots 
or the virtual machines. Special configuration files can be used to build 
docker containers. Storing all configuration tasks in a structured ``Dockerfile`` 
reduces the risk of data loss and allows us to move the workers between different 
machines. 

.. warning::
  Remember that **each** command in the ``Dockerfile`` leads to the creation of a 
  **new** layer of files that cannot be deleted. Be careful here and try to bundle 
  commands to save disk space and to keep garbage out of the image. See 
  http://jrruethe.github.io/blog/2015/09/20/dockerfile-generator/ and 
  https://docs.docker.com/articles/dockerfile_best-practices/ for more good
  advice on this topic.

Some more useful commands when working with docker are::

    docker stop `docker ps -aq`; docker rm `docker ps -aq`; #delete all docker containers
    docker rmi `docker images -f "dangling=true" -q`; #delete all dangling docker images

The workflow to generate the images locally could look like::

    git clone --recursive https://github.com/CoolProp/Dockerfiles.git CoolProp.Dockerfiles.git
    cd CoolProp.Dockerfiles.git
    cd workerbase/64bit      ; docker build -t coolprop/workerbase      -f Dockerfile . ; cd ..
    cd workerpython/64bit    ; docker build -t coolprop/workerpython    -f Dockerfile . ; cd ..
    cd workerlinuxopen/64bit ; docker build -t coolprop/workerlinuxopen -f Dockerfile . ; cd ..

Please also have a look at the CoolProp repository on Docker Hub to see which 
images are available for download https://hub.docker.com/r/coolprop/ and do not hesitate to 
contribute to the sources at https://github.com/CoolProp/Dockerfiles

