#!/bin/bash

#log_file=('/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/bash_scripts/test_output.txt')
#beagle_results=('/Users/lester/Documents/PhD')
#names=('0A2744C' '1A2744P' '2M0416C' '3M0416P' '4M0717C' '5M0717P' '6M1149C' '7M1149P')
#names=('0A2744C')
#sub_folder=('test')
#slurm_folder=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run')
#empty_limit=('50')

# first 4 files
#log_file=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/slurm_automation_log.txt')
#beagle_results=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results')
##names=('0A2744C')
#names=('0A2744C' '1A2744P' '2M0416C' '3M0416P')
#sub_folder=('fit_001')
#slurm_folder=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run')
#empty_limit=('20')

# all 8 files
# log_file=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/slurm_automation_log.txt')
# beagle_results=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results')
# names=('0A2744C' '1A2744P' '2M0416C' '3M0416P' '4M0717C' '5M0717P' '6M1149C' '7M1149P')
# sub_folder=('fit_001')
# slurm_folder=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run')
# empty_limit=('20')

# all 8 files
log_file=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/new_beagle/danger_delayed_uniform_logtau/slurm_automation_log.txt')
beagle_results=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/danger_delayed_uniform_logtau')
names=('0' '1' '2' '3' '4' '5' '6' '7')
sub_folder=('fit_001')
slurm_folder=('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/new_beagle/danger_delayed_uniform_logtau')
empty_limit=('20')

echo "     " >> "$log_file"
echo "SCRIPT START" >> "$log_file"

now=$(date)
echo "$now" >> "$log_file"

squeue="$(/usr/local/software/slurm/current/bin/squeue -u ls861 | wc -l)"
# exit if there is something in squeue
if [ $squeue -ne '1' ]; then
  echo "SQUEUE $squeue lines" >> "$log_file"
  exit 555
fi

for name in "${names[@]}"
do

  complete_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size +0k | wc -l)"
  echo "  ${name} # complete .gz files: $complete_gz" >> "$log_file"

  # if for this field, there are no complete .gz files, then run 50 nodes
  if [ $complete_gz -eq '0' ]; then
    echo "  Running slurm file for $name" >> "$log_file"
    /usr/local/software/slurm/current/bin/sbatch ${slurm_folder}/slurm_field_${name}_fit_001
    # exit 666
  fi

  # if for this field, there are some complete .gz files, then count empty .gz files and .lock files
  if [ $complete_gz -gt '0' ]; then
    empty_gz="$(find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size 0k | wc -l)"
    empty_lock="$(find ${beagle_results}/${name}/${sub_folder} -name '*.lock' -type 'f' | wc -l)"
    echo "  ${name} # empty .gz files: $empty_gz (limit: $empty_limit)" >> "$log_file"
    echo "  ${name} # empty .lock files: $empty_lock (limit: $empty_limit)" >> "$log_file"

    # if there are more than the limit of empty .gz files, then remove empty .gz and lock files, and then run 50 nodes
    if [ $empty_gz -gt $empty_limit ]; then
      echo "  Deleting lock and empty .gz files" >> "$log_file"
      find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size 0k -delete
      find ${beagle_results}/${name}/${sub_folder} -name '*.lock' -type 'f' -delete
      echo "  Re-running slurm file for $name" >> "$log_file"
      /usr/local/software/slurm/current/bin/sbatch ${slurm_folder}/slurm_field_${name}_fit_001
      # exit 777

    # elif [ $empty_lock -gt $empty_limit ]; then
    #   echo "  Deleting lock and empty .gz files" >> "$log_file"
    #   find ${beagle_results}/${name}/${sub_folder} -name '*.gz' -type 'f' -size 0k -delete
    #   find ${beagle_results}/${name}/${sub_folder} -name '*.lock' -type 'f' -delete
    #   echo "  Re-running slurm file for $name" >> "$log_file"
    #   /usr/local/software/slurm/current/bin/sbatch ${slurm_folder}/slurm_field_${name}_fit_001
    #   exit 888
    fi
  fi

done

# echo "SCRIPT COMPLETE" >> "$log_file"
