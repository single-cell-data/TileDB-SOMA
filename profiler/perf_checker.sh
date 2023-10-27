#!/bin/bash

python -m venv perf
source perf/bin/activate
pip install gitpython
pip install psutil
pip install comacore
pip install profiler
pip install tiledbsoma
pip install cellxgene_census
python -m profiler "python ann_data.py" -t gtime

python ./top_profiler.py