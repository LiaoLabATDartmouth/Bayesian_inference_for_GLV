#!/bin/bash

cwd="$(pwd)"
base="BF_interact"
cmdstan="$HOME/.cmdstan/cmdstan-2.35.0"
make "$cwd/$base" -C "$cmdstan"

# run three replicates
for i in 1 2 3
do
    ./$base sample num_samples=10000 num_warmup=1000 data file=${base}.data.json output file=${base}_${i}.csv &
done

wait
${cmdstan}/bin/stansummary ${base}_*.csv > summary_123.txt

