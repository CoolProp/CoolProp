.. _Delphi:

**************
Delphi Wrapper
**************

Delphi `(information) <http://www.embarcadero.com/products/delphi>`_ is a programming environment based on the Pascal language.  It is used for doing GUI development and other general programming.  There is an open-source equivalent of Delphi called `Lazarus <http://www.lazarus.freepascal.org/>`_ .

A very simple example of linking Lazarus and CoolProp can be found obtained by downloading the following 3 files:

* :download:`CoolPropDelphiExample.lpi`
* :download:`CoolPropDelphiExample.lpr`
* :download:`CoolPropDelphiExample.lps`

You will need to download the release shared library for your platform from :sfdownloads:`sourceforge <shared_library>` or the development version of the shared library from :sfnightly:`the nightly snapshots <shared_library>`.

Open Lazarus and then open the project file you downloaded.  Build and run. Place the downloaded shared library for CoolProp in the same folder with generated executable.


Floating Point Exceptions
-------------------------

We have noticed a run time error that occurs with certain inputs. For example, call to

    t := PropsSI('T', 'P', p, 'H', h, ref);

This will crash if the dll is called from Delphi (not confirmed for Lazarus) at the point p = 4.4e6, h=330e3 and ref=R22.

If you disable the floating point exceptions, then you get the right answer!

    Set8087CW($133f);

Don't get me wrong here, this is not a fix! Apparently, most compilers disable the floating point exceptions so all we are doing here is masking the actual fault.
