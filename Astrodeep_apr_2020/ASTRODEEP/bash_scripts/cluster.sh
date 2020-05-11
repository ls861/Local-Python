#!/bin/bash

#names = ('0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P')
names=('0A2744C')

# CHECK IF QUEUE IS EMPTY HERE, BREAK IF NOT, PRINT TIME AND AN ERROR

for name in "${names[@]}"
do
  #echo 'test'
  #turd = $(find /Users/lester/Documents/PhD/${name}/fit_001 -name '*.gz' -type 'f' -size 0k | wc -l)
  #${turd} > /Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/bash_scripts/test_output.txt

  now=$(date)
  echo "$now"

  test1=$(find /Users/lester/Documents/PhD/${name}/fit_001 -name '*.gz' -type 'f' -size 0k | wc -l)
  echo "$test1"


done

#find /Users/lester/Documents/PhD/${name}/fit_001 -name '*.gz' -type 'f' -size 0k | wc -l > /Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/bash_scripts/test_output.txt
#find /Users/lester/Documents/PhD/${name}/fit_001 -name '*.gz' -type 'f' -size +0k | wc -l > /Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/bash_scripts/test_output.txt
