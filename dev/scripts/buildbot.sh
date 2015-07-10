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
function git_cfg {
    pushd /home/$USER/buildbot/CoolProp.git
    git reset HEAD dev/buildbot/master/master.cfg
    git checkout -- dev/buildbot/master/master.cfg
    popd
} 
function stop {
    buildbot stop /home/$USER/buildbot/server-master/
}
function clean {
    rm -f /home/$USER/buildbot/server-master/buildbot_private.pyc
    python /home/$USER/buildbot/server-master/buildbot_private.py
}
#
# Check for input
CMD="$1"
if [ "$CMD" = "restart" ]; then
    stop
    git_pull
    clean
    start
elif [ "$CMD" = "reconfig" ]; then
    git_pull
    clean
    reconfig
elif [ "$CMD" = "reconfigmaster" ]; then
    git_cfg
    git_pull
    clean
    reconfig
elif [ "$CMD" = "start" ]; then
    git_pull
    clean
    start
elif [ "$CMD" = "stop" ]; then
    stop
else
    git_pull
    clean
    start
fi
#
echo "Script executed, terminating"
exit 0
