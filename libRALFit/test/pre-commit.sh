#!/bin/bash
# Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
# All rights reserved.

# tests to be run before a commit
# to enable, run
#  ln -s pre-commit.sh ../../.git/hooks/pre-commit

# stash uncommitted code
git stash -q --keep-index

# run the tests
./run_tests.sh
RESULT=$?

# and restore code base
git stash pop -q

[ $RESULT -ne 0 ] && exit 1
exit 0
