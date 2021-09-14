#!/bin/sh
# See https://stackoverflow.com/questions/6059336/how-to-find-the-current-git-branch-in-detached-head-state
$1 git describe --contains --all HEAD | rev | cut -d/ -f1 | rev
