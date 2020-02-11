#!/bin/bash
# Script to download ASTRODEEP tar files and then extract them

wget -nc http://ia2-guest.oats.inaf.it/wp-content/uploads/2016/02/A2744cl_26012016.tar.gz
wget -nc http://ia2-guest.oats.inaf.it/wp-content/uploads/2016/02/A2744PAR_26012016.tar.gz
wget -nc http://ia2-guest.oats.inaf.it/wp-content/uploads/2016/02/M0416cl_26012016.tar.gz
wget -nc http://ia2-guest.oats.inaf.it/wp-content/uploads/2016/02/M0416PAR_26012016.tar.gz
wget -nc http://www.astrodeep.eu/wp-content/uploads/2017/02/MACS0717cl_10022017.tar.gz
wget -nc http://www.astrodeep.eu/wp-content/uploads/2017/02/MACS0717par_10022017.tar.gz
wget -nc http://www.astrodeep.eu/wp-content/uploads/2017/02/MACS1149cl_10022017.tar.gz
wget -nc http://www.astrodeep.eu/wp-content/uploads/2017/02/MACS1149par_10022017.tar.gz

for i in *.tar.gz; do echo working on $i; tar xvzf $i ; done
