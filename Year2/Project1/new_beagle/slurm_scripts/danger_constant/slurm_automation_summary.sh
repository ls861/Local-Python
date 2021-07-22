#!/bin/bash

# all 8 files
beagle_results=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/danger_constant')
names=('0' '1' '2' '3' '4' '5' '6' '7')
totals=('3043' '3070' '2974' '3354' '2866' '2641' '2935' '3106')
sub_folder=('fit_002')
empty_limit=('50')

for name in "${names[@]}"
do
  total_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' | wc -l)"
  complete_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size +0k | wc -l)"
  empty_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size 0k | wc -l)"
  limit_gz="$((${totals[${name}]}-${empty_limit}))"
  echo "  field ${name} # max/total/complete/empty/limit .gz files: ${totals[${name}]} $total_gz $complete_gz $empty_gz $empty_limit"
done
