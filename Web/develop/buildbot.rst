
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
