#!/bin/bash
#DEVBINDIR="/home/coolprop/buildbot/server-master/public_html/binaries"
#find "$DEVBINDIR" -type d ! -perm -a+rx -exec chmod -v a+rx {} \;
#find "$DEVBINDIR" -type f ! -perm -a+r -exec chmod -v a+r {} \;
DEVBINDIR="/home/coolprop/buildbot/server-master/public_html"
chmod -R a+rX,u+w,go-w "$DEVBINDIR"
exit 0