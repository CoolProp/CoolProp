***************
CoolProp Online
***************

The online version of CoolProp is stored at https://github.com/CoolProp/CoolProp-Online.  It is an application for web2py.  It is deployed to PythonAnywhere, and served from there.  To get it working:

Getting going on PythonAnywhere
-------------------------------

1. In PythonAnywhere, add a new web app, select your username
2. Make sure that you can get to https://YOURUSERNAME.pythonanywhere.com and you don't get an error (ok, good, web2py is running). Your admin page is https://ibell.pythonanywhere.com/admin/default/index
2. Open a shell in ``/home/YOURUSERNAME/web2py/applications``
3. ``git clone https://github.com/CoolProp/CoolProp-Online coolpropgit`` to check out the application
4. Reset the web application in Web tab in PythonAnywhere, go to admin webpage (https://YOURUSERNAME.pythonanywhere.com/admin/default/index), make sure you can see ``coolpropgit`` application
5. Deposit a file called ``routes.py`` in ``/home/YOURUSERNAME/web2py`` with the contents::

    routers = dict(
        BASE = dict(
            default_application='coolpropgit',
        )
    )

6. Restart the webservice
7. Page should serve properly, and redirect to the appropriate page

Running application locally
---------------------------

.. warning::
    
    Make sure you install ``mpld3``!

Same basic protocol...

1. Check out the source of web2py (http://www.web2py.com/init/default/download)
2. Open shell in unzipped web2py folder
3. ``cd applications``
4. ``git clone https://github.com/CoolProp/CoolProp-Online coolpropgit`` to check out the application
5. ``cd ..``
6. ``python web2py.py``

Useful information
------------------

* How to update web2py safely: https://www.pythonanywhere.com/forums/topic/881/#id_post_6624
* Info on deploying on PythonAnywhere: http://web2py.com/books/default/chapter/29/13#Deploying-on-PythonAnywhere
* Setting default app: http://www.ridgesolutions.ie/index.php/2013/02/20/web2py-make-your-app-the-default-web-application/
