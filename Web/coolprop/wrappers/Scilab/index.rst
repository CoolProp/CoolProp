.. _Scilab:

**************
Scilab wrapper
**************

There are at least three different ways of interfacing with CoolProp from scilab:

* Call the default CoolProp shared library directly from CoolProp
* Call CoolProp using the PIMS module which calls CoolProp through Python
* SWIG wrapper for CoolProp (no longer actively supported)

Calling the shared library directly from scilab is the most straightforward to configure and get working, and requires no other scilab modules, but the only two functions that can be called are the high-level functions ``PropsSI`` and ``HAPropsSI``

For more complete coverage of CoolProp functionality, you can use the PIMS interface (see below)

Shared Library from Scilab
==========================
Download a precompiled shared library appropriate to the computer you are using from :sfdownloads:`sourceforge <shared_library>` or the development version from :sfnightly:`the nightly snapshots <shared_library>`.  Or you could build it yourself following these instructions :ref:`build your own shared library <shared_library>`.
    
Download the sample file: 

* :download:`sample.sce`

Place both files in the same folder.

Run the sample file using something like this at the scilab prompt::

    exec('sample.sce',-1)
    
The ``-1`` is to ensure that all lines are not echoed (unelated question: why is echoing all lines in a script the default behavior in scilab?)

If you always use CoolProp, you might find it convenient to have the shared library load automatically when you start Scilab. To do that, just rename the 'sample.sce' file to 'scilab.ini' and put it in the scilab home directory which can be found by executing 'SCIHOME' in the scilab console window, the result will be something like that:

For Linux:

-->SCIHOME
 SCIHOME  =
 
 /home/username/.Scilab/scilab-5.5.2

For Windows:

-->SCIHOME
 SCIHOME  =
 
 C:\\Users\\username\\AppData\\Roaming\\Scilab\\scilab-5.5.2

.. note:: 

    It is possible that on linux, you might need to add the path to the directory containing ``libCoolProp.so`` to the ``LD_LIBRARY_PATH`` environmental variable.  At the terminal, or in your ~/.bash_profile file, add: ``export LD_LIBRARY_PATH=/path/to/directory``. It is also possible to put the ``libCoolProp.so`` in '/usr/lib/scilab' folder, so you don't need to specify the path in the 'scilab.ini' file.

Calling through PIMS
====================
`PIMS <https://atoms.scilab.org/toolboxes/PIMS>`_ is an interface layer that allows you to call functions/modules/packages in Python.  It can be installed with the command at a scilab prompt::

    atomsInstall('PIMS')

.. warning::

    Only python 2.7 is supported.  The default python on your system must be python 2.7 and not python 3.4.  I have no idea if you can switch which python it uses, it doesn't seem so.
    
    Also, numpy is required, but you should be using anaconda anyway, which includes numpy

Once it is installed, make sure that it can find python, and that the python it finds is a 2.7.x version of python::

    pyEvalStr('import sys; print sys.version')
    
which should yield something like (OS-dependent)::

    2.7.11 |Continuum Analytics, Inc.| (default, Dec  7 2015, 14:10:42) [MSC v.1500 64 bit (AMD64)]
    
Then, you can use CoolProp just like in Python. Here is an example of calling the low-level interface to get the normal-boiling-point temperature of water::

    pyImport CoolProp.CoolProp
    AS = CoolProp.AbstractState('HEOS','Water')
    AS.update(CoolProp.CoolProp.PQ_INPUTS, 101325, 1)
    AS.T()

Or the same thing from the high-level interface::

    pyImport CoolProp.CoolProp
    CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')
