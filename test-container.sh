#!/bin/bash

TEST_IMAGE="click_mergevcfs_test_image"

cd $( dirname "${BASH_SOURCE[0]}" )

if [ "$1" = "--skip-build" ]; then
    echo "skipping build..."
else
    echo "building image, to skip run with --skip-build"
    docker build -t $TEST_IMAGE .
fi

# see https://explainshell.com/explain?cmd=set+-euxo%20pipefail
set -euxo pipefail

# There was issues that docker container has a memory limit of 1.95Gb,
# while the callers requires at least 2Gb of memory.
# To raise container memory limit, see
# https://stackoverflow.com/questions/32834082/how-to-increase-docker-machine-memory-mac

# remove pytest cache
echo "testing docker image..."
find . -name '*.pyc' -exec rm {} +
find . -name '__pycache__' -exec rm -rf {} +

docker run --rm -it $TEST_IMAGE --version
docker run --rm -it --entrypoint ''  -v `pwd`:/test -w /test  \
    $TEST_IMAGE bash -c 'cp -r /test /tmp && cd /tmp/test/ && pip install pytest && pytest tests' # && cp .coverage /test'

# move container coverage paths to local, see .coveragerc [paths] and this comment:
# https://github.com/pytest-dev/pytest-cov/issues/146#issuecomment-272971136
# echo "combining container coverage..."
# command -v coverage > /dev/null 2>&1 || pip install coverage
# mv .coverage .coverage.tmp
# coverage combine --append

echo "tests finished..."
