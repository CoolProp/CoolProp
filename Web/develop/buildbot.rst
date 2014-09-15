
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