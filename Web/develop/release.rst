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
      ``python milestone2rst.py check vX.X.X``. This command finds the date
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
    - Review the generated text from ``snippet_issues.rst.txt`` and 
	  ``snippet_pulls.rst.txt`` and update the changelog file in
      ``Web/coolprop/changelog.rst``. You might also use the same text for
      the annotated tag / release on GitHub.
	  
* **Delete cache**: Delete the cached property plots and consistency plots stored
  at https://github.com/CoolProp/CoolProp/actions/caches. This forces their 
  recreation when building the documentation. Recalculating all the consistency
  plots is our last line of defence when anything goes wrong with an EoS prior
  to a scheduled release.    
* **Push to master**: Merge your changes to the *master* branch and wait for the 
  CI system to complete the work. Only proceed if all builds finish successfully.
  If anything goes wrong, you should be able to debug the workflows locally
  using act (https://github.com/nektos/act/).
* **Tag a release**: Tag the master branch for using ``vX.X.X`` and wait once more 
  for all CI actions to complete. Make sure that the Python wheels get uploaded 
  to PyPi automatically.
* **Release**: Wait for all actions to finish and manually launch the release action
  with the version number vX.X.X as input. This updates the homepage and uploads the
  binaries to SourceForge. 
* **Clean up**: If everything went well, you can proceed: 
    - Create a new release on GitHub using the vX.X.X tag. 
	- Add a note with the download link: https://sourceforge.net/projects/coolprop/files/CoolProp/X.X.X/
    - Change the default download file on SourceForge to point to the new
      zipped sources.
    - Bump the version number in the CMake file and commit.
    - Announce the new features if you like.

That's all folks.
