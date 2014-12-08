#!/bin/bash
# Work around for Cron:
USER=coolprop
source /home/$USER/.bash_profile
source /home/$USER/buildbot/server-master-sandbox/bin/activate
#
function start {
    buildbot start /home/$USER/buildbot/server-master/
}
function reconfig {
    buildbot reconfig /home/$USER/buildbot/server-master/
}
function git_pull {
    pushd /home/$USER/buildbot/CoolProp.git
    if git pull; then
        echo "Updated Git repository"
    else
        echo "\"git pull\" failed, aborting."
        exit 1
    fi
    popd
}
function stop {
    buildbot stop /home/$USER/buildbot/server-master/
}
#
# Check for input
CMD="$1"
if [ "$CMD" = "restart" ]; then
    stop
    git_pull
    start
elif [ "$CMD" = "reconfig" ]; then
    git_pull
    reconfig
elif [ "$CMD" = "start" ]; then
    git_pull
    start
elif [ "$CMD" = "stop" ]; then
    stop
else
    git_pull
    start
fi
#
echo "Script executed, terminating"
exit 0
