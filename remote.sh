#!/usr/bin/env bash
# remote.sh: the only entry point for running commands in the cluster repository
# Usage:
#   ./remote.sh run    '<command>'          run in the foreground and return its exit code
#   ./remote.sh bg     <name> '<command>'   run in the background, logged to results/logs/<name>.log
#   ./remote.sh status <name>               show background job status and the end of its log
set -euo pipefail

HOST=Ai43
RDIR=/root/project/brain_atlas   # absolute path of the repository on the cluster
LOGDIR=results/logs

remote() {
  # PATH is set explicitly because a non-interactive login shell does not read ~/.bashrc
  local setup='export PATH="$HOME/.local/bin:$PATH"'
  ssh "$HOST" "bash -lc $(printf %q "$setup; cd $RDIR && $1")"
}

check_name() {
  [[ "$1" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "job name may contain only letters, digits, underscores and hyphens" >&2; exit 2; }
}

case "${1:-}" in
  run)
    remote "$2"
    ;;
  bg)
    check_name "$2"
    name=$2
    remote "mkdir -p $LOGDIR
if [ -f $LOGDIR/$name.pid ] && kill -0 \$(cat $LOGDIR/$name.pid) 2>/dev/null; then
  echo 'a job with this name is still running' >&2; exit 1
fi
echo \"start \$(date -Is) code \$(git describe --always --dirty 2>/dev/null)\" > $LOGDIR/$name.log
nohup setsid bash -c $(printf %q "$3") >> $LOGDIR/$name.log 2>&1 < /dev/null &
echo \$! > $LOGDIR/$name.pid
echo \"started pid \$!\""
    ;;
  status)
    check_name "$2"
    name=$2
    remote "pid=\$(cat $LOGDIR/$name.pid 2>/dev/null || true)
if [ -z \"\$pid\" ] || ! ps -p \$pid -o pid,etime,%cpu,rss,args; then echo 'job is not running'; fi
tail -n 40 $LOGDIR/$name.log"
    ;;
  *)
    sed -n '2,6p' "$0"; exit 2
    ;;
esac
