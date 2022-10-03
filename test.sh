#!/bin/bash
julia --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'

parallel --joblog /tmp/log -j8 jupyter nbconvert --to notebook --ExecutePreprocessor.timeout=600 --execute --inplace {} ::: docs/*.ipynb
cat /tmp/log

parallel -j8 jupyter nbconvert --clear-output --inplace {} ::: docs/*.ipynb
