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
  
    - Move into the ``dev/scripts`` folder and run ``python milestone2rst.py check vX.X.X``. This command finds the date of the latest release and looks at all issues that have been closed since then. It reports problems such as missing labels.
    - Take the time to fix all problems and label issues and PRs.

* **Version**: Edit CMakeLists.txt and remove all qualifiers (alpha, dev,
  ...) from the version number and update the number for the upcoming
  release.
* **Changelog**: Update the changelog and generate a list of closed GitHub
  issues: 
  
    - Move into the ``dev/scripts`` folder and do ``python milestone2rst.py changelog vX.X.X`` where ``vX.X.X`` is the version number of the milestone on GitHub.
    - Review the generated text from ``snippet_issues.rst.txt`` and ``snippet_pulls.rst.txt`` and update the changelog file in ``Web/coolprop/changelog.rst``. You might also use the same text for the annotated tag / release on GitHub.
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
    - Remove artifacts from TestPYPI. First with ``pypi-cleanup -t https://test.pypi.org -p coolprop --query-only -r ".*\.dev\d.*"``
      to determine what will be deleted and then with your username ``-u username --do-it``
      (and without ``--query-only``) to do it. The ``.dev`` pattern matches the
      ``X.Y.Z.dev<timestamp>`` nightlies that CI publishes to TestPyPI.

That's all folks.

*****************************
Desktop GUI release (Tauri)
*****************************

The cross-platform desktop GUI (``wrappers/GUI``, a Tauri + React app) is
**released on its own ``gui-vX.X.X`` tag** by the ``gui_builder.yml``
workflow — independently of the library's ``vX.X.X`` tag and the
SourceForge/PyPI flow above. A ``gui-v*`` tag builds all three platforms,
code-signs them, and creates a single **draft** GitHub Release with the
installers plus the ``latest.json`` auto-updater manifest attached.

Signing is engaged on every ``gui-v*`` tag (and on any build whose HEAD
commit message contains the ``[gui-sign]`` marker). macOS is signed +
notarized; Windows is signed through SignPath. **Routine pushes build
unsigned** and send no signing email.

**Preconditions (verify before tagging — all fail *open* to an
unsigned-but-functional build, so check them explicitly):**

* **GUI version**: the bundle is stamped from the **committed** version in
  ``wrappers/GUI/package.json`` / ``src-tauri/tauri.conf.json`` /
  ``src-tauri/Cargo.toml`` (plus their lockfiles). A build-time sync
  (``scripts/generate-licenses.mjs``, run via ``beforeBuildCommand``) derives
  the version from ``CMakeLists.txt`` and rewrites these files — but Tauri
  reads ``tauri.conf.json`` *before* that sync runs, so the sync does **not**
  fix the in-flight bundle; the committed value is what lands on the
  ``.dmg`` / ``.msi`` / ``.deb``. So bump the four version fields (run
  ``npm run build:licenses`` in ``wrappers/GUI`` after bumping
  ``CMakeLists.txt``, or edit them directly) **and commit them** before
  tagging. Use plain semver (e.g. ``8.0.0``); a CMake ``REVISION`` qualifier
  like ``b1`` is not valid for Cargo/npm and is dropped by the sync.
* **Apple agreements current**: notarization returns HTTP 403 *"a required
  agreement is missing or has expired"* if either the **Apple Developer
  Program License Agreement** (developer.apple.com/account, accepted by the
  Account Holder) or the **App Store Connect → Business** agreement has
  lapsed. Clear any pending banner first.
* **Apple secrets present**: ``APPLE_CERTIFICATE`` /
  ``APPLE_CERTIFICATE_PASSWORD`` / ``APPLE_SIGNING_IDENTITY`` /
  ``APPLE_TEAM_ID`` and the App Store Connect API-key trio
  (``APPLE_API_ISSUER`` / ``APPLE_API_KEY`` / ``APPLE_API_KEY_BASE64``).
* **Auto-updater secrets present**: ``TAURI_SIGNING_PRIVATE_KEY`` (+
  password) — without them the release ships but in-app updates won't
  verify.
* **Windows / SignPath**: the repo variable ``SIGNPATH_POLICY_SLUG`` must be
  ``release-signing`` (it defaults to the untrusted ``test-signing`` cert,
  which does **not** clear SmartScreen), the org secret
  ``SIGNPATH_API_TOKEN`` must be set, and the ``CI builds`` submitter must
  be on the ``release-signing`` policy's submitters list.

**Cutting the release:**

#. Push the tag: ``git tag gui-vX.X.X && git push origin gui-vX.X.X``.
#. **Approve the Windows signing request.** The ``release-signing`` policy
   requires a manual approval (1 approver). When the Windows job reaches the
   SignPath step it submits the request and **blocks Pending approval** —
   approve it in the SignPath dashboard
   (``coolprop`` → ``release-signing``) within the action's
   ``wait-for-completion-timeout-in-seconds`` (currently 3600 s / 1 h) or
   the job times out and fails. (A real approval in validation took ~41 min,
   so don't trim this much.)
#. Let all three platform builds finish. The release job's ``latest.json``
   step is **all-or-nothing** — a single failed/unsigned platform aborts the
   whole release, so a macOS notarization or Windows approval failure means
   *no* release is produced.
#. The ``Overlay signed Windows installers onto bundle`` step prints each
   Windows installer's Authenticode ``Status`` — confirm ``Valid`` (the
   release cert), not ``UnknownError`` (the test cert) or ``NotSigned``.
#. Edit the **draft** release body (testers only see what's there) and
   publish. Publishing is what arms the in-app auto-updater
   (``releases/latest/download/latest.json``), so don't publish a build you
   don't want pushed to existing users.

.. note::

   The draft's "prerelease" flag is auto-set only for tags containing
   ``-alpha`` / ``-beta`` / ``-rc``. A PEP 440-style ``gui-v8.0.0b1`` is
   **not** matched — tick *Set as a pre-release* by hand before publishing a
   beta.
