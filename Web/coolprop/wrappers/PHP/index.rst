.. _PHP:

***********
PHP Wrapper
***********

Pre-Compiled Binaries
=====================

* Download the shared library for your architecture from :sfdownloads:`PHP`, or from the development buildbot server at :sfnightly:`PHP`.  Also download the CoolProp.php platform-independent file.  Or build it yourself (see below)

* Copy the libCoolProp.so file into the extension-dir for php::

    sudo cp libCoolProp.so `php-config --extension-dir`

  This copies the shared library into a location that PHP can load.  sudo is needed to make the copy. Alternatively, in the following step you will need to use the full path to the extension, which is just a bit more annoying.

  Put the CoolProp.php with the rest of your sources.  This is the interface file between your PHP code and the CoolProp module.

* Modify the PHP.ini file that PHP will load to add::

    extension = "libCoolProp.so"

  after ``[PHP]``. If you didn't copy libCoolProp.so into the folder given by ```php-config --extension-dir``` you will need to use the absolute path

* You can determine the php.ini file that you should be modifying by creating a file on the server with the contents

  .. code-block:: php

     <?php
       phpinfo();
     ?>

  Modify the file given by ``Loaded Configuration File``.

* To enable some useful debugging at runtime (turn off for deployment), you can add

  .. code-block:: php

     <?php
     error_reporting  (E_ALL);
     ini_set ('display_errors', true);
     ...
     ?>

* This stub file shows how to call the module

  .. code-block:: php

     <?php
     require "CoolProp.php";

     $p = 101325;
     $Q = 1.0;
     $T = PropsSI("T","P",$p,"Q",$Q,"Water");
     print "NBP of water is $T\n";
     ?>
     
 * And here is another example demonstrating how to call the low-level interface:
 
   .. code-block:: php
   
     <?php
        include_once "CoolProp.php";

        // set the REFPROP path
        set_config_string(ALTERNATIVE_REFPROP_PATH, "/opt/refprop/");

        // extend the abstract class AbstractState
        class ConcreteState extends AbstractState {
            static function factory($backend, $fluid_names) {
                $r = AbstractState_factory($backend, $fluid_names);
                if (!is_resource($r)) return $r;
                return new ConcreteState($r);
            }
        }

        // instantiate the class
        $water = ConcreteState::factory("REFPROP", "Water");

        $p = 101325;
        $Q = 1.0;
        $water->update(PQ_INPUTS, $p, $Q);
        $T = $water->T();

        print "NBP of water is $T\n";
    ?>

User-Compiled Binaries
======================

Since most servers that will serve PHP will be linux, only linux instructions are shown here.  In principle very minor changes should be required to get this working on windows.

Common Requirements
-------------------
Compilation of the php wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Additionally, you need SWIG, which can be obtained on Debian-based OS with::

    sudo apt-get install swig

On Ubuntu 16.04 LTS, the default PHP is 7.0, which is not compatible with SWIG as of Aug 10, 2016.  So to get PHP 5.6 on ubuntu 16.04, you can do::

    sudo add-apt-repository ppa:ondrej/php
    sudo apt-get update
    sudo apt-get install php5.6 php5.6-dev

Linux
-----

1. Check out CoolProp::

    git clone https://github.com/CoolProp/CoolProp --recursive

2. Folder creating::

    mkdir -p CoolProp/build && cd CoolProp/build

3. Build the php module::

    cmake .. -DCOOLPROP_PHP_MODULE=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON

4. Build (verbosely so we can see if there are any problems)::

    cmake --build .

  This will generate the file libCoolProp.so and the php module CoolProp.php

5. See the above instructions in the Pre-Compiled Binaries section for installation instructions
