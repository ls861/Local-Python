#!/bin/bash

# all 8 files
log_file=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/new_beagle/danger_delayed_uniform_logtau/slurm_automation_log.txt')
beagle_results=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/danger_delayed_uniform_logtau')
names=('0' '1' '2' '3' '4' '5' '6' '7')
totals=('3043' '3070' '2974' '3354' '2866' '2641' '2935' '3106')
sub_folder=('fit_002')
slurm_folder=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/new_beagle/danger_delayed_uniform_logtau')
empty_limit=('50')

echo "     " >> "$log_file"
echo "SCRIPT START" >> "$log_file"

now=$(date)
echo "$now" >> "$log_file"

squeue="$(/usr/local/software/slurm/current/bin/squeue -u ls861 -r | wc -l)" # -r expands pending
# exit if there is more than 300 pending or running tasks in squeue
if [ $squeue -gt '251' ]; then
  echo "SQUEUE $squeue lines" >> "$log_file"
  exit 555
fi

for name in "${names[@]}"
do
  total_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' | wc -l)"
  complete_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size +0k | wc -l)"
  empty_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size 0k | wc -l)"
  limit_gz="$((${totals[${name}]}-${empty_limit}))"
  echo "  field ${name} # max/total/complete/empty/limit .gz files: ${totals[${name}]} $total_gz $complete_gz $empty_gz $empty_limit" >> "$log_file"

  # if there are not many empty + complete gz files, run some more nodes
  if [ $total_gz -lt $limit_gz ]; then
    echo "  Running slurm file for field $name" >> "$log_file"
    /usr/local/software/slurm/current/bin/sbatch ${slurm_folder}/slurm_field_${name}_fit_001
    exit 666
  fi

done

echo "SCRIPT COMPLETE" >> "$log_file"
