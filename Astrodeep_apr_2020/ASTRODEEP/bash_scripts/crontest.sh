#!/bin/bash

log_file=('/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/bash_scripts/cron_test_output.txt')

echo "     " >> "$log_file"
echo "THIS IS A TEST" >> "$log_file"

now=$(date)
echo "$now" >> "$log_file"

echo "SCRIPT COMPLETE" >> "$log_file"
