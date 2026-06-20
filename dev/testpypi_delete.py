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
    python dev/testpypi_delete.py --keep 10 --browser safari --do-it
    PYPI_SESSION_COOKIE='session=...; ...' python dev/testpypi_delete.py --do-it

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

    if cookie_string:
        host = base.split("//", 1)[1]
        for part in cookie_string.split(";"):
            part = part.strip()
            if not part or "=" not in part:
                continue
            name, value = part.split("=", 1)
            session.cookies.set(name.strip(), value.strip(), domain=host)
        return session

    try:
        import browser_cookie3
    except ImportError:
        raise SystemExit(
            "browser_cookie3 is not installed and no --session-cookie was given.\n"
            "  Install it:  python -m pip install browser_cookie3\n"
            "  or paste a cookie:  PYPI_SESSION_COOKIE='session=...' python dev/testpypi_delete.py ...")

    host = base.split("//", 1)[1]
    loader = getattr(browser_cookie3, browser, None)
    if loader is None:
        raise SystemExit(f"Unknown --browser {browser!r}; expected chrome/safari/firefox/edge/brave.")
    try:
        jar = loader(domain_name=host)
    except Exception as exc:  # decryption / Keychain / locked-DB failures
        raise SystemExit(
            f"Could not read {host} cookies from {browser}: {exc}\n"
            "  If macOS Keychain access was denied, re-run and click 'Always Allow',\n"
            "  or paste the cookie instead:  PYPI_SESSION_COOKIE='session=...' python dev/testpypi_delete.py ...")
    session.cookies.update(jar)
    if not list(session.cookies):
        raise SystemExit(
            f"No {host} cookies found in {browser}. Log in to {base} in {browser} first, "
            "then re-run.")
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


def delete_release(session, base, package, version, do_it):
    """Delete a single release via the real web confirm-delete POST.

    In dry-run mode we still GET the release page (proving the session can reach
    it and the delete form is present) but never POST.  Returns a status string.
    """
    action = f"/manage/project/{package}/release/{version}/"
    url = f"{base}{action}"

    r = session.get(url, timeout=30)
    if r.status_code == 404:
        return "already-gone"
    r.raise_for_status()
    if "login" in r.url or "two-factor" in r.url:
        raise SystemExit(
            f"Session lost authentication while fetching {action} (redirected to {r.url}). "
            "Re-log in to TestPyPI in your browser and retry.")

    parser = _DeleteFormCsrf(action)
    parser.feed(r.text)
    csrf = parser.csrf or parser.any_csrf
    if not csrf:
        raise SystemExit(
            f"No csrf_token found on {action} -- the page served was probably not the "
            f"release-management page (auth/verification issue). Final URL: {r.url}")

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
            f"Browser session is not authenticated to {remote} (landed on {final}).\n"
            f"  Log in to {base} in {args.browser} (complete any email verification), then re-run.")
    print(f"Authenticated to {remote} via {('pasted cookie' if args.session_cookie else args.browser)} session. OK.")

    if not args.do_it:
        print(f"\nDRY RUN -- would delete {len(delete_keys)} releases (re-run with --do-it):")
    else:
        print(f"\n!!! DELETING {len(delete_keys)} releases from {remote} !!!")

    failures = 0
    for i, version in enumerate(delete_keys, 1):
        try:
            status = delete_release(session, base, args.package, version, args.do_it)
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
