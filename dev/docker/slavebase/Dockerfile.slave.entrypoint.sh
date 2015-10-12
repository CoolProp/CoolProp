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
  # Remove the security relevant variable from the environment
  unset MASTERHOST
  unset SLAVEPASSWORD
  unset BOTADMIN
  unset BOTEMAIL
  unset BOTHOST
  $CTRLAPP start ${SLAVEDIR}
  echo "Remember to update the SSH configuration:"
  echo "docker cp \${HOME}/.ssh ${SLAVENAME}:${HOME}/"
  echo "docker cp \${HOME}/.pypirc ${SLAVENAME}:${HOME}/"
  echo "docker exec ${SLAVENAME} chown -R buildbot ${HOME}/"
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
