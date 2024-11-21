#!/usr/bin/env bash

# Bash script to run HMM3 searches

#new to this version is the use of the Bash script ini_file_parser.sh written by WOLF Software and available via the MIT license at https://github.com/DevelopersToolbox/ini-file-parser

#check which unix platform (linux or OSX)
if [[ "$OSTYPE" == "linux"* ]]; then
    # for Linux machines
    #$Bin="$(dirname $(readlink -f -- "$0"))"
    Bin="$( readlink -f -- "${0%/*}" )"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # for OSX maachines
    # readlink does not work for OSX
    Bin="$(cd -- "$(dirname "$0")" && pwd)"
fi

# set the phage home directory
phome='/opt/Phage_Finder'

# set the base or working directory
base=`pwd`

# set peptide filename
pepfile=$1

# Load in the ini file parser
source ${phome}/bin/ini-file-parser.sh

# Load and process the ini/config file
process_ini_file "${phome}/phage_finder.ini"

# set paths for executibles from the phage_finder.ini file

hmmsearch=$(get_value 'hmmsearch' 'hmmsearch')

if [ ! -d HMM3_searches_dir ] # check if directory is present
then
    mkdir ${base}/HMM3_searches_dir
else
    rm -rf ${base}/HMM3_searches_dir/
    mkdir ${base}/HMM3_searches_dir
fi
cd ${base}/HMM3_searches_dir/
total=`cat ${phome}/hmm3.lst | wc -l | sed 's/^ *//'`
for i in `cat ${phome}/hmm3.lst`
do
  let "count = count + 1"
  Result=`echo $count $total | awk '{printf( "%3.1f\n", (($1/$2) * 100))}'`
  echo -ne "\r    $Result% complete"
  #echo -ne "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b$Result% complete"
  ${hmmsearch} ${phome}/PHAGE_HMM3s_dir/${i}.HMM ${base}/${pepfile} >> ${i}.out
done
echo
cat *.out >> ../combined.hmm3
cd ..
rm -rf HMM3_searches_dir/
