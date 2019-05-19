.. _release:

******************
Release Checklist
******************

We have made a serious effort to automate the release of new binaries. Even
though things have become easier, there are still many things to remember.
Here is your new best friend, a checklist that helps you to keep track of
all the small things that have to be done when releasing a new version of
the CoolProp library. 

* **Merge to master**: Close issues and merge pull requests where
  appropriate. Merge all relevant changes from your feature branches into
  the *master* branch.
* **Issues and PRs**: Make sure that all issues that were closed since the last
  release are labelled as described below or have a milestone attached to
  them. The label "invalid" means that the reported
  issue was not related to CoolProp, whereas "duplicate" means that the
  issue has been reported earlier. Finally, the label "wontfix" means that
  we do not consider the issue for future work and closed it. The latter
  could also apply to real bugs that we cannot fix. All of these issues
  should not be mentioned in the changelog. We have created a script that
  check the issues for you: 
  
    - Move into the ``dev/scripts`` folder and run
      ``python milestone2rst.py check v6.3.0``. This command finds the date
      of the latest release and looks at all issues that have been closed
      since then. It reports problems such as missing labels.
    - Take the time to fix all problems and label issues and PRs.

* **Version**: Edit CMakeLists.txt and remove all qualifiers (alpha, dev,
  ...) from the version number and update the number for the upcoming
  release.
* **Changelog**: Update the changelog and generate a list of closed GitHub
  issues: 
  
    - Move into the ``dev/scripts`` folder and do ``python milestone2rst.py
      changelog vX.X.X`` where ``vX.X.X`` is the version number of the
      milestone on GitHub.
    - Copy the generated text (goes to stdout) into the changelog file in
      ``Web/coolprop/changelog.rst``. You might also use the same text for
      the annotated tag / release on GitHub.
    
* **Merge to release**: Merge *master* into *release* branch.
* **Build Bots**: Force all buildbots to run on the *release* branch, this
  will also change the upload folder from *binaries* to *release*.
* **Release**: Wait for all bots to finish and run the release script by
  launching the ``release version`` bot with dry run disabled and the
  correct version number. This uploads binaries to pypi and sourceforge.
  Ignore the warning ``failed to set times on 
  "/home/project-web/coolprop/htdocs/jscript/coolprop-latest.js"``,
  it is a symlink and will be overwritten. If you encounter problems, log
  in via SSH and have a look at the logs. If you would like to finished the
  release manually, consider editing the case statement in bash script and
  run the commands from ``release.bsh.cmds.txt`` manually.
* **Clean and Tag**: If everything went well, you can proceed: 

    - Create a new tag and a new release on GitHub. Remember to
      make an annotated tag and include the information on the closed
      issues here as well. 
    - Change the default download file on sourceforge to point to the new
      zipped sources.
    - Bump the version number in the CMake file and commit.
    - Announce the new features if you like...

That's all folks.
