#!/bin/sh
mkdir tmp || true

HYPERGRAPH_FILE=$1
P=$2
N_TRIALS=$3
CUTOFF=$4

if [ -z $P ]
then
  echo "Expected P variable"
  exit 1
fi

cp $HYPERGRAPH_FILE tmp/hypergraph.json

if [ -z $HYPERGRAPH_FILE ]
then
  echo "Expected hypergraph file"
  exit 1
fi


if [ -z $N_TRIALS ]
then
  echo "Expected n trials variable"
  exit 1
fi

if [ -z $CUTOFF ]
then
  echo "Expected cutoff variable"
  exit 1
fi

PARAMS="{\"p\": $P, \"n_trials\": $N_TRIALS, \"cutoff\": $CUTOFF, \"file_out\": \"out/$HYPERGRAPH_FILE-$P-$N_TRIALS-$CUTOFF.json\"}"
echo $PARAMS > tmp/params.json
