
Pushing files to the release channels
=====================================

Get your keys ready
-------------------

CoolProp uses SourceForge to distribute binaries and the compressed source code 
archives. If you want to be part of the release team, you have to be registered 
`SourceForge <http://ww.sourceforge.net>`_ user. 

For seamless integration, you need to set up the ssh keys to enable automated 
``rsync`` file transfers. Use ``ssh-keygen -C "myownname@shell.sf.net"`` to generate
a pair of ssh keys on the machine you would like to use during the release process.
This would typically be the buildbot master, see also :ref:`Buildbot`.

Log in to the SourceForge web site and go to ``Me -> Account Settings -> SSH Settings``.
Set your login shell to ``/bin/bash`` and paste the content of the public key 
file into the text box. The file you want to use to that is normally ``~/.ssh/id_rsa.pub``.

