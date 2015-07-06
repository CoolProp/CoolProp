.. _release:

******************
Release Checklist
******************

We have made a serious effort to automate the release of new binaries. Even 
though things have become easier, there are still many things to remember. 
Here is your new best friend, a checklist that helps you to keep track of all
the small things that have to be done when releasing a new version of the CoolProp 
library. 

* **Version**: Edit CMakeLists.txt and remove all qualifiers (alpha, dev, ...) from the version number.
* **Changelog**: Update the changelog and generate a list of closed GitHub issues: 

    - Move into the ``dev/scripts`` folder and do ``python milestone2rst.py vX.X.X`` where ``vX.X.X`` is the version number of the milestone on github.
    - Copy the generated text (goes to stdout) into the changelog file in ``Web/coolprop/changelog.rst``
    
* **release branch**: Merge all code from *master* into *release* branch
* **build bots**: Force all buildbots to run on the *release* branch, this will also change the upload folder from *binaries* to *release*.
* **script**: Wait for all bots to finish and run the release script by launching the ``release version`` bot with dry run disabled and the correct version number. This uploads binaries to pypi and sourceforge. Ignore the warning ``failed to set times on "/home/project-web/coolprop/htdocs/jscript/coolprop-latest.js"``, it is a symlink and will be overwritten... 
* **clean up**: If everything went well, you can proceed: 

    - Tag the release branch in your git software. It is a good idea to make an annotated tag and include the information on the closed issues here as well. 
    - Change the default download file on sourceforge to point to the new zipped sources.
    - Bump the version number in the CMake file and commit. 
    - Announce the new features if you like...
  
That's all folks.
