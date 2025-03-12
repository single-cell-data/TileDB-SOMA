#!/usr/bin/env bash

jq -r '.mem | map_values(.peak_mem) | to_entries[] | (.key + " " + (.value | tostring))' "$1" | \
parallel -k -j+0 --env PATH --colsep ' ' 'echo {1} `numfmt --to=iec {2}`'
