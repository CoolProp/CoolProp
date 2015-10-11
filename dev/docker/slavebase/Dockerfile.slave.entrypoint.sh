#!/bin/bash
CTRLAPP="/usr/local/bin/buildslave"
function shutdown()
{
  $CTRLAPP stop ${SLAVEDIR}
  exit 0
}

function startup()
{
  if [ ! -d "${SLAVEDIR}" ]; then
    /usr/local/bin/buildslave create-slave ${SLAVEDIR} ${MASTERHOST} ${SLAVENAME} ${SLAVEPASSWORD}
    echo "${BOTADMIN} <${BOTEMAIL}>" > ${SLAVEDIR}/info/admin
    echo "${BOTHOST}" > ${SLAVEDIR}/info/host
  fi
  $CTRLAPP start ${SLAVEDIR}
}

trap shutdown TERM SIGTERM SIGKILL SIGINT

startup;

# Just idle for one hour and keep the process alive
# waiting for SIGTERM.
while : ; do
sleep 3600 & wait
done
#
echo "The endless loop terminated, something is wrong here."
exit 1
