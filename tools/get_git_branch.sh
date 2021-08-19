#!/bin/sh
# See https://stackoverflow.com/questions/6059336/how-to-find-the-current-git-branch-in-detached-head-state
$1 for-each-ref --format='%(objectname) %(refname:short)' refs/heads | awk "/^$($1 rev-parse HEAD)/ {print \$2}"
