#!/bin/bash

cwd="$(pwd)"
base=`basename "$cwd"`
cmdstan="$HOME/Desktop/cmdstan-2.33.1/"
make "$cwd/$base" -C "$cmdstan"

# run three replicates
for i in 1 2 3
do
    ./$base sample num_samples=1000 num_warmup=100 data file=${base}.data.json output file=${base}_${i}.csv &
done

