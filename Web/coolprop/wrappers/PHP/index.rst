.. _PHP:

***********
PHP Wrapper
***********

Pre-compiled Binaries
=====================
Pre-compiled PHP binaries are no longer supplied.  The last release with PHP files on SourceForge is 6.1.0 (2016); see `its PHP folder <https://sourceforge.net/projects/coolprop/files/CoolProp/6.1.0/PHP/>`_ for the last known binaries.  Those files were presumably built for PHP 5.x and are not expected to load in current PHP, because a compiled extension is tied to the exact PHP version it was built for.  The wrapper is built from source (see below), which takes a few minutes.

User-Compiled Binaries
======================

Since most servers that will serve PHP will be Linux, only Linux instructions are shown here.  In principle very minor changes should be required to get this working on Windows.

Common Requirements
-------------------
Compilation of the php wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Additionally, SWIG and the PHP command-line program with its development files (which provide ``php-config``) are needed.  On Ubuntu 26.04 (PHP 8.5, SWIG 4.4), these are installed with::

    sudo apt-get install swig php-cli php-dev

The wrapper was built and run with this combination on Ubuntu 26.04 (WSL).  Other Linux distributions and PHP versions are not tested.  The extension is built for one specific PHP version (and thread-safety setting), so it needs to be rebuilt after PHP is upgraded.

Linux
-----

1. Check out CoolProp::

    git clone https://github.com/CoolProp/CoolProp

2. Folder creating::

    mkdir -p CoolProp/build && cd CoolProp/build

3. Build the php module::

    cmake .. -DCOOLPROP_PHP_MODULE=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON

4. Build (verbosely so we can see if there are any problems)::

    cmake --build . -j4

   This will generate the file ``CoolPropPHP.so`` in the build folder.  With SWIG 4.4 no ``CoolProp.php`` file is generated.  A single-threaded build (without ``-j4``) takes about 8 minutes.

   .. note::

     With SWIG 4.1 and newer the PHP classes and functions are registered natively and
     the separate ``CoolProp.php`` proxy file is no longer generated.  Code written for
     older SWIG versions that starts with ``require "CoolProp.php"`` needs that
     line removed, since the file no longer exists.  Builds made with SWIG older than 4.1 are not tested.

5. Check that the module loads and returns a value (the full path to the module is given, so nothing needs to be installed)::

    php -d extension=$PWD/CoolPropPHP.so -r 'echo PropsSI("T","P",101325,"Q",1.0,"Water"), "\n";'

   The expected output is ``373.12429584767`` (the normal boiling point of water, in K).

Loading the Module
------------------

PHP loads the module through the ``extension`` setting.  Any of the following works; only the last one needs ``sudo``.

* **Command line, nothing installed.**  Pass the full path to the module with ``-d`` (as in step 5 above), or from a script, run ``php -d extension=/full/path/to/CoolPropPHP.so script.php``.

* **A user-owned configuration file.**  Create a file ``coolprop.ini`` containing the single line ``extension=/full/path/to/CoolPropPHP.so`` and start PHP with ``php -n -c coolprop.ini script.php``.  Note that ``-n`` skips the system ``php.ini`` and its ``conf.d`` folder, so other extensions are not loaded unless they are listed in this file as well.

* **An additional scan directory.**  Put ``coolprop.ini`` (as above) in a folder of the user's choice and set ``PHP_INI_SCAN_DIR`` to that folder followed by the system's own scan directory, so the other extensions stay enabled::

    PHP_INI_SCAN_DIR=$HOME/php-ini:$(php -r 'echo PHP_CONFIG_FILE_SCAN_DIR;') php script.php

* **System-wide, for a web server.**  Copy the module into PHP's extension directory and enable it in the ``php.ini`` that the server uses::

    sudo cp CoolPropPHP.so `php-config --extension-dir`

  Then add the line ``extension=CoolPropPHP.so`` to the ``php.ini`` (or a file in the ``conf.d`` folder next to it) of each PHP flavor that needs it (``php -i | grep "Loaded Configuration File"`` shows the file for the command-line program; the web server's PHP, such as PHP-FPM or Apache, has its own).  The file used by a web server can be found by creating a page on the server with the contents

  .. code-block:: php

     <?php
       phpinfo();
     ?>

  and looking at ``Loaded Configuration File``.  The web server (or PHP-FPM) needs to be restarted afterwards.  This route was not tested with a web server.

The ``-c`` and ``PHP_INI_SCAN_DIR`` routes were tested with the PHP command-line program; both returned the same value as above.  They only apply to the program they are used with, which makes them suitable for scripts and cron jobs but not for a web server that is already running.

Examples
--------

* To enable some useful debugging at runtime (turn off for deployment), the following can be added at the top of a script

  .. code-block:: php

     <?php
     error_reporting  (E_ALL);
     ini_set ('display_errors', true);
     ...
     ?>

* This stub file shows how to call the module (tested with PHP 8.5 and SWIG 4.4; no ``require`` is needed)

  .. code-block:: php

     <?php
     $p = 101325;
     $Q = 1.0;
     $T = PropsSI("T","P",$p,"Q",$Q,"Water");
     print "NBP of water is $T\n";
     ?>

* And here is another example demonstrating how to call the low-level interface (same tested versions):

  .. code-block:: php

     <?php
        $water = AbstractState::factory("HEOS", "Water");

        $p = 101325;
        $Q = 1.0;
        $water->update(PQ_INPUTS, $p, $Q);
        $T = $water->T();

        print "NBP of water is $T\n";
     ?>

* Errors (for example an unknown fluid name) are thrown as PHP ``Exception`` objects, so they can be handled with ``try``/``catch``.

* If REFPROP is installed in a non-standard location, the path is set with (REFPROP itself was not available for testing here, only that the setting is accepted)

  .. code-block:: php

     <?php
        set_config_string(config_string_to_key("ALTERNATIVE_REFPROP_PATH"), "/opt/refprop/");
     ?>
