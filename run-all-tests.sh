#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
GEN_DIR="${SCRIPT_DIR}/build-release/tools"
TARGET_DIR="${SCRIPT_DIR}/tests/data"
TESTS_DIR="${SCRIPT_DIR}/build-release/tests"

$SCRIPT_DIR/build.sh "generator AllTests"

SMALL_GRAPH="$TARGET_DIR/gnp_100k_400k.metis"
BIG_GRAPH="$TARGET_DIR/gnp_1m_4m.metis"

if [ ! -f "$SMALL_GRAPH" ]; then
    echo "generationg small gnp graph"
    $GEN_DIR/generator --output=$SMALL_GRAPH gnp -n 100000 -a 4 -d 1
else
    echo "skip generationg small gnp graph -> file exists already"
fi

if [ ! -f "$BIG_GRAPH" ]; then
    echo "generationg big gnp graph"
    $GEN_DIR/generator --output=$BIG_GRAPH gnp -n 1000000 -a 4 -d 1
else
    echo "skip generationg big gnp graph -> file exists already"
fi

echo "running tests"
cd $TESTS_DIR
./AllTests


