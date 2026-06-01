/* Injects a dismissable Sponsor banner at the top of every docs page.
   Dismissal is remembered in localStorage (versioned key) so it does not
   reappear on subsequent page loads. Theme-agnostic: prepends to <body>. */
(function () {
  "use strict";
  var DISMISS_KEY = "coolprop.sponsorBanner.dismissed.v1";
  var SPONSOR_URL = "https://github.com/sponsors/CoolProp";

  function dismissed() {
    try { return window.localStorage.getItem(DISMISS_KEY) === "1"; }
    catch (e) { return false; }
  }
  function remember() {
    try { window.localStorage.setItem(DISMISS_KEY, "1"); } catch (e) { /* ignore */ }
  }

  function build() {
    if (dismissed() || document.getElementById("coolprop-sponsor-banner")) return;

    var bar = document.createElement("div");
    bar.id = "coolprop-sponsor-banner";

    var text = document.createElement("span");
    text.className = "cp-sponsor-text";
    text.innerHTML =
      "💚 <strong>CoolProp is free &amp; open-source</strong>, " +
      "maintained on volunteer time. If it helps your work, please consider sponsoring →";

    var cta = document.createElement("a");
    cta.className = "cp-sponsor-cta";
    cta.href = SPONSOR_URL;
    cta.target = "_blank";
    cta.rel = "noopener";
    cta.textContent = "Sponsor";

    var x = document.createElement("button");
    x.className = "cp-sponsor-dismiss";
    x.type = "button";
    x.setAttribute("aria-label", "Dismiss sponsor banner");
    x.innerHTML = "&times;";
    x.addEventListener("click", function () {
      remember();
      if (bar.parentNode) bar.parentNode.removeChild(bar);
    });

    bar.appendChild(text);
    bar.appendChild(cta);
    bar.appendChild(x);
    document.body.insertBefore(bar, document.body.firstChild);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", build);
  } else {
    build();
  }
})();
