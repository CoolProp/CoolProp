"""Bulk-delete old TestPyPI (or PyPI) CoolProp dev releases -- *locally*, with
zero credential handling, by reusing the session you already have open in your
browser.

Why this exists
---------------
Deletion on PyPI/TestPyPI is NOT possible with an upload API token; it needs a
real authenticated web session.  The old CI approach (``pypi-cleanup`` driven by
a stored password + a ``pyotp``-generated TOTP code) stopped working on
2025-11-14, when PyPI shipped *email verification for TOTP logins from new
devices* (https://blog.pypi.org/posts/2025-11-14-login-verification/).  A fresh
GitHub Actions runner is always an "unrecognized device", so after the TOTP code
is accepted PyPI serves an email-confirmation page instead of completing the
login -- and the deletion step then dies with a misleading
``No CSFR found in /manage/project/CoolProp/release/...`` (see pypi-cleanup
issues #42 / #48).  No amount of better scraping fixes that: an ephemeral runner
cannot click the link in the verification email.

The workaround is to stop authenticating from scratch.  Your normal browser is
already logged in *and already a trusted device*, so this script lifts the
TestPyPI session cookies straight out of it (Chrome by default) and replays the
exact delete-confirmation POST that the web UI sends.  No password, no TOTP, no
stored secret -- just "stay logged in to TestPyPI in your browser, then run me".

Usage
-----
::

    # Dry run (default): show what would be deleted, verify the session works.
    python dev/testpypi_delete.py --keep 10

    # Actually delete the old dev builds.
    python dev/testpypi_delete.py --keep 10 --do-it

    # Pull the session from a different browser, or paste a cookie explicitly
    # if browser extraction is unavailable (e.g. Keychain access denied).
    python dev/testpypi_delete.py --keep 10 --browser firefox --do-it
    PYPI_SESSION_COOKIE='session_id=...' python dev/testpypi_delete.py --do-it

If anything goes wrong, the script prints explicit, copy-pasteable recovery steps
(re-login, the one-time macOS Keychain grant, or the paste-a-cookie fallback).

The delete *set* is computed by the same read-only logic the CI prune job uses
(:mod:`testpypi_prune`), so final/tagged releases are never candidates.
"""
import os
import sys
import time
import argparse
from html.parser import HTMLParser

import requests

# Single source of truth for "which versions are prune candidates".  This file
# lives next to testpypi_prune.py in dev/.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from testpypi_prune import fetch_release_keys, select_delete_set  # noqa: E402

# How this script is invoked, for copy-pasteable advice in error messages.
SCRIPT = "dev/testpypi_delete.py"


def _host(base):
    """'https://test.pypi.org' -> 'test.pypi.org'."""
    return base.split("//", 1)[1]


def _login_help(base, browser):
    """Actionable recovery steps when the browser session isn't logged in."""
    return (
        "\nHOW TO FIX -- the browser session is not logged in to this site:\n"
        f"  1. Open {base}/manage/projects/ in {browser} and log in\n"
        "     (username + password + TOTP).  You should NOT be asked for email\n"
        "     verification if this browser is already a trusted device -- the\n"
        "     'remember_device' cookie skips it.\n"
        "  2. Confirm you can see your project(s) listed on that page.\n"
        f"  3. Re-run the same command:  python3 {SCRIPT} ...\n"
        "\nKEYCHAIN-FREE ALTERNATIVE (no browser cookie reading at all):\n"
        f"  In {browser} DevTools -> Application -> Cookies -> {base}, copy the\n"
        "  'session_id' cookie value, then run:\n"
        f"    PYPI_SESSION_COOKIE='session_id=<value>' python3 {SCRIPT} ...\n"
        "  (paste it fresh each time -- the session cookie is short-lived.)"
    )


def _keychain_help(base, browser):
    """Actionable recovery steps when the browser cookie store can't be read."""
    host = _host(base)
    return (
        f"\nHOW TO FIX -- macOS could not read {browser}'s encrypted cookies:\n"
        "  - This needs a ONE-TIME Keychain grant.  Run this in your own terminal\n"
        "    (foreground, so the dialog is visible) and click 'Always Allow' on\n"
        "    the 'Chrome Safe Storage' prompt:\n"
        f"      python3 -c \"import browser_cookie3 as b; print(len(list(b.chrome(domain_name='{host}'))), 'cookies')\"\n"
        f"  - Then re-run:  python3 {SCRIPT} ...\n"
        "  - Safari is NOT a usable fallback here: macOS blocks its cookie file\n"
        "    unless your terminal has Full Disk Access.\n"
        "\nOR SKIP KEYCHAIN ENTIRELY (paste the cookie):\n"
        f"  In {browser} DevTools -> Application -> Cookies -> {base}, copy the\n"
        "  'session_id' cookie value, then run:\n"
        f"    PYPI_SESSION_COOKIE='session_id=<value>' python3 {SCRIPT} ..."
    )


class _DeleteFormCsrf(HTMLParser):
    """Extract the CSRF token Warehouse expects for a release deletion.

    Warehouse renders the delete-release dialog as a ``<form>`` whose ``action``
    is the release-management path and which contains a hidden ``csrf_token``
    input plus a ``confirm_delete_version`` text input (see Warehouse's
    ``manage/manage_base.html`` ``confirm_modal`` macro).  The token itself is
    the *session* CSRF token -- identical for every form on the page -- so we
    prefer the one inside the matching delete form but fall back to the first
    ``csrf_token`` seen anywhere, which guards against future markup churn.
    """

    def __init__(self, action_target):
        super().__init__()
        self._target = action_target
        self._in_form = False
        self._form_has_confirm = False
        self._form_csrf = None
        self.csrf = None       # token from the matching delete form (preferred)
        self.any_csrf = None   # first csrf_token anywhere (fallback)

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "form":
            action = a.get("action") or ""
            self._in_form = action == self._target or action.startswith(self._target)
            self._form_has_confirm = False
            self._form_csrf = None
            return
        if tag == "input":
            if a.get("name") == "csrf_token":
                value = a.get("value")
                if self.any_csrf is None:
                    self.any_csrf = value
                if self._in_form:
                    self._form_csrf = value
            elif a.get("name") == "confirm_delete_version" and self._in_form:
                self._form_has_confirm = True

    def handle_endtag(self, tag):
        if tag == "form" and self._in_form:
            if self._form_has_confirm and self._form_csrf and self.csrf is None:
                self.csrf = self._form_csrf
            self._in_form = False


def load_session(browser, cookie_string, base):
    """Build a requests.Session carrying the browser's logged-in cookies.

    Cookie source precedence: an explicit ``cookie_string`` (``--session-cookie``
    or ``$PYPI_SESSION_COOKIE``) wins; otherwise we extract cookies for the host
    from the named browser via ``browser_cookie3``.  Returns the session.
    """
    session = requests.Session()
    session.headers["User-Agent"] = "coolprop-prune/1.0 (local browser-session reuse)"

    host = _host(base)

    if cookie_string:
        for part in cookie_string.split(";"):
            part = part.strip()
            if not part or "=" not in part:
                continue
            name, value = part.split("=", 1)
            session.cookies.set(name.strip(), value.strip(), domain=host)
        if "session_id" not in {c.name for c in session.cookies}:
            print("WARNING: the pasted cookie has no 'session_id' -- auth will likely fail.\n"
                  "         Copy the 'session_id' cookie value, not some other cookie.",
                  file=sys.stderr)
        return session

    try:
        import browser_cookie3
    except ImportError:
        raise SystemExit(
            "browser_cookie3 is not installed and no --session-cookie / $PYPI_SESSION_COOKIE was given."
            + _keychain_help(base, browser))

    loader = getattr(browser_cookie3, browser, None)
    if loader is None:
        raise SystemExit(
            f"Unknown --browser {browser!r}; expected one of: chrome, firefox, edge, brave, opera, safari.")
    try:
        jar = loader(domain_name=host)
    except Exception as exc:  # decryption / Keychain / locked-DB failures
        raise SystemExit(f"Could not read {host} cookies from {browser}: {exc}" + _keychain_help(base, browser))
    session.cookies.update(jar)
    if not list(session.cookies):
        raise SystemExit(f"No {host} cookies found in {browser}." + _login_help(base, browser))
    return session


def verify_authenticated(session, base):
    """Confirm the session is a logged-in, trusted web session.

    Returns True iff GET /manage/projects/ stays on the management page (not a
    login or device-verification redirect).
    """
    r = session.get(f"{base}/manage/projects/", timeout=30)
    r.raise_for_status()
    final = r.url
    ok = "/manage/projects" in final and "login" not in final and "two-factor" not in final
    return ok, final


def delete_release(session, base, package, version, do_it, browser):
    """Delete a single release via the real web confirm-delete POST.

    In dry-run mode we still GET the release page (proving the session can reach
    it and the delete form is present) but never POST.  Returns a status string.
    """
    # The release-management path uses the PEP 440-*normalized* version, which is
    # exactly what select_delete_set returns (str(version.parse(...))).  If a raw
    # JSON key were ever non-normalized the URL wouldn't match the stored path --
    # but that fails closed: the GET 404s and we report "already-gone" rather than
    # deleting some other release.
    action = f"/manage/project/{package}/release/{version}/"
    url = f"{base}{action}"

    r = session.get(url, timeout=30)
    if r.status_code == 404:
        return "already-gone"
    r.raise_for_status()
    if "login" in r.url or "two-factor" in r.url:
        raise SystemExit(
            f"Session lost authentication while fetching {action} (redirected to {r.url})."
            + _login_help(base, browser))

    parser = _DeleteFormCsrf(action)
    parser.feed(r.text)
    csrf = parser.csrf or parser.any_csrf
    if not csrf:
        raise SystemExit(
            f"No csrf_token found on {action} -- the page served was probably a login or "
            f"device-verification page, not the release-management page (final URL: {r.url})."
            + _login_help(base, browser))

    if not do_it:
        return "would-delete"

    resp = session.post(
        url,
        data={"csrf_token": csrf, "confirm_delete_version": version},
        headers={"Referer": url, "Origin": base},
        timeout=30,
    )
    resp.raise_for_status()
    # A successful delete redirects away from the (now-gone) release page.
    if resp.url.rstrip("/").endswith(f"/release/{version}"):
        return f"FAILED (still on release page: {resp.url})"
    return "deleted"


def main():
    parser = argparse.ArgumentParser(
        description="Bulk-delete old TestPyPI/PyPI dev releases by reusing your browser's logged-in session.")
    parser.add_argument("--package", default="CoolProp", help="package name (default: CoolProp)")
    parser.add_argument("--keep", type=int, default=10,
                        help="number of newest development (non-final) builds to retain (default: 10)")
    parser.add_argument("--pypi", action="store_true",
                        help="target real PyPI instead of TestPyPI (DANGEROUS)")
    parser.add_argument("--browser", default="chrome",
                        help="browser to pull the session cookie from (chrome/safari/firefox/edge/brave)")
    parser.add_argument("--session-cookie", default=os.environ.get("PYPI_SESSION_COOKIE"),
                        help="paste a 'name=value; name2=value2' cookie string instead of reading the browser "
                             "(or set $PYPI_SESSION_COOKIE)")
    parser.add_argument("--do-it", action="store_true",
                        help="actually delete (omit for a dry run)")
    parser.add_argument("--delay", type=float, default=1.0,
                        help="seconds to pause between deletions (default: 1.0)")
    args = parser.parse_args()

    base = "https://pypi.org" if args.pypi else "https://test.pypi.org"
    remote = "PyPI" if args.pypi else "TestPyPI"

    keys = fetch_release_keys(args.package, pypi=args.pypi)
    keep_keys, delete_keys = select_delete_set(keys, args.keep)
    print(f"{remote} {args.package}: {len(keys)} total releases; "
          f"keeping {len(keep_keys)} newest dev builds, {len(delete_keys)} candidates to delete.")
    if not delete_keys:
        print("Nothing to delete; already within the keep window.")
        return 0

    session = load_session(args.browser, args.session_cookie, base)
    authed, final = verify_authenticated(session, base)
    if not authed:
        raise SystemExit(
            f"Browser session is NOT authenticated to {remote} (landed on {final})."
            + _login_help(base, args.browser))
    print(f"Authenticated to {remote} via {('pasted cookie' if args.session_cookie else args.browser)} session. OK.")

    if not args.do_it:
        print(f"\nDRY RUN -- would delete {len(delete_keys)} releases (re-run with --do-it):")
    else:
        # Deleting from the real PyPI is irreversible and production-facing, so
        # require an explicit typed confirmation on top of --do-it.  (TestPyPI
        # deletes proceed with just --do-it.)
        if args.pypi:
            answer = input(f"About to DELETE {len(delete_keys)} releases from PRODUCTION {remote}. "
                           f"Type the package name {args.package!r} to confirm: ")
            if answer.strip() != args.package:
                raise SystemExit("Confirmation did not match; aborting without deleting anything.")
        print(f"\n!!! DELETING {len(delete_keys)} releases from {remote} !!!")

    failures = 0
    for i, version in enumerate(delete_keys, 1):
        try:
            status = delete_release(session, base, args.package, version, args.do_it,
                                    "pasted cookie" if args.session_cookie else args.browser)
        except requests.RequestException as exc:
            status = f"ERROR ({exc})"
        if status.startswith(("FAILED", "ERROR")):
            failures += 1
        print(f"  [{i}/{len(delete_keys)}] {version}: {status}")
        if args.do_it and args.delay and i < len(delete_keys):
            time.sleep(args.delay)

    if failures:
        print(f"\n{failures} deletion(s) failed.")
        return 1
    print("\nDone." if args.do_it else "\nDry run complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
