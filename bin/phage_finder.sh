#!/usr/bin/env bash

# This Bash script is the Phage_Finder pipeline that runs BLASTP searches, HMMer3 searches, tRNAscan-SE, Aragorn and cmscan

# Usage: phage_finder.sh <prefix of .pep/.faa, .ptt and .con/.fna file>

# NOTE: a phage_finder_info.txt file will be searched for before a .ptt file
# .pep/.faa is the multifasta protein sequence file
# .ptt is a GenBank .ptt file that has the coordinates and ORF names with annotation. 
# .con/.fna is a file that contains the complete nucleotide sequence of the genome being searched
# If .ptt is not avaliable, Phage_finder_info.txt can be used. It can be generated from the BVBRC annotation files using the phage_finder_info-generator-BVBRC_annot.py script. 

#very cool Perl FindBin functional equivalent for BASH by Nicholas Dronen of crashingdaily.wordpress.com

#new to this version is the use of the Bash script ini_file_parser.sh written by WOLF Software and available via the MIT license at https://github.com/DevelopersToolbox/ini-file-parser
#also added the ability to perform blastp searches using diamond
NO_ARGS=0
OPTERROR=65

usage() { echo "Usage: $0 [-s <blastp or diamond>] <prefix of .pep/.faa, .ptt and .con/.fna file>" 1>&2; exit 1; }

if [ $# -eq "$NO_ARGS" ]; then
    usage
fi

while getopts ":hs:" Option
do
  case $Option in
    s) s=${OPTARG};;
    h) usage;;
  esac
done
shift $(($OPTIND - 1))  # Decrements the argument pointer so it points to next argument.

#if no value for option s or not equal to blastp or diamond, exit with usage
if ! [[ "$s" =~ ^(blastp|diamond)$ ]]; then
    usage
fi

#check which unix platform (linux or OSX)
if [[ "$OSTYPE" == "linux"* ]]; then
    # for Linux machines
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

# Load in the ini file parser
source ${phome}/bin/ini-file-parser.sh

# Load and process the ini/config file
process_ini_file "${phome}/phage_finder.ini"

# set paths for executibles from the phage_finder.ini file

blastp=$(get_value 'pipeline' 'blastp')
diamond=$(get_value 'pipeline' 'diamond')
tRNAscan=$(get_value 'pipeline' 'tRNAscan')
aragorn=$(get_value 'pipeline' 'aragorn')
seqstat=$(get_value 'pipeline' 'seqstat')
cmscan=$(get_value 'pipeline' 'cmscan')
fig2dev=$(get_value 'pipeline' 'fig2dev')

prefix=$1

  if [ -s ${base}/${prefix}.pep ] # check if .pep file is present
  then
      pepfile="${prefix}.pep"
  elif [ -s ${base}/${prefix}.faa ]
  then
      pepfile="${prefix}.faa"
  else
     echo "Could not find a ${prefix}.pep or ${prefix}.faa file.  Please check to make sure the file is present and contains data"
     exit 1
  fi 
    if [ -s ${base}/phage_finder_info.txt ] # check for phage_finder info file and if it has contents
    then
        infofile="phage_finder_info.txt"
    elif [ -s ${base}/${prefix}.ptt ]
    then
          infofile="${prefix}.ptt"
    else
      echo "Could not find a phage_finder_info.txt file or $prefix.ptt file.  Please make sure one of these files is present and contains data."
      exit 1
    fi
    if [ ! -s ${base}/combined.hmm3 ] # if GLOCAL HMM results not present, search
    then
        ## conduct GLOCAL HMMER3 searches
        echo "  HMMER3 searches ..."
        ${phome}/bin/HMM3_searches.sh ${pepfile}
    fi
    if [ "${s}" = "blastp" ] # if blastp searches desired
    then
        echo "BLASTP searches requested"
        if [ ! -s ${base}/ncbi.out ] # if BLAST results not present, search
        then
            ## do NCBI BLASTP searches
            echo "  BLASTing ${pepfile} against the Phage DB ..."
            ${blastp} -db ${phome}/DB/phage_03_25_19.db -outfmt 6 -evalue 0.001 -query ${pepfile} -out ncbi.out -max_target_seqs 5 -num_threads 6
        fi
    elif [ "${s}" = "diamond" ] # if diamond searches desired
    then
        echo "Diamond searches requested"
        if [ ! -s ${base}/diamond.out ] # if DIAMOND results not present, then search
        then
            ## do DIAMOND BLASTP searches
            echo "  Searching ${pepfile} against the Phage DB using Diamond..."
            ${diamond} blastp --db ${phome}/DB/phage_03_25_19_diamond.db.dmnd --outfmt 6 --evalue 0.001 --query ${pepfile} --out diamond.out
        fi
    fi
    if [ -s ${base}/${prefix}.con ]
    then
        contigfile="${prefix}.con"
    elif [ -s ${base}/${prefix}.fna ]
    then
        contigfile="${prefix}.fna"
    else
        echo "Could not find a phage_finder_info.txt file or ${prefix}.ptt file.  Please make sure one of these files is present and contains data.  In the meantime, I will go ahead and run phage_finder.pl without this information, but beware... NO att sites will be found!"
        contigfile=""
    fi
    if [ ! -s ${base}/tRNAscan.out ] && [ -s ${base}/${contigfile} ] # if tRNAscan.out file not present, and contig file present, then search
    then
        ## find tRNAs
        echo "  find tRNA sequences ..."
        ${tRNAscan} -B -o tRNAscan.out ${base}/${contigfile} > /dev/null 2>&1
    fi

    if [ ! -s ${base}/tmRNA_aragorn.out ] && [ -s ${base}/${contigfile} ] # if tRNAscan.out file not present, and contig file present, then search
    then
        ## find tmRNAs
        echo "  find tmRNA sequences ..."
        ${aragorn} -m -o tmRNA_aragorn.out ${base}/${contigfile}
    fi

    if [ ! -s ${base}/ncRNA_cmscan.out ] && [ -s ${base}/${contigfile} ] # if ncRNA_cmscan.out file not present, and contig file present, then search
    then
        ## find select ncRNAs
        echo "  find ncRNA sequences ..."
        Z=`${seqstat} ${base}/${contigfile}| perl -ne 'chomp; if (/^Total # residues:\s+(\d+)/) {$n = $1; $Z=($n*2)/1000000; print "$Z\n";}'`
        #echo "Z = $Z"        
        ${cmscan} -Z $Z --cut_ga --rfam --nohmmonly --tblout ${base}/ncRNA_cmscan.out --fmt 2 ${phome}/RfamDB/Rfam_PhageFinder.cm ${base}/${contigfile} > ${prefix}.cmscan
    fi

    ## find the prophage
    echo "  searching for Prophage regions ..."
    if [ "${s}" = "blastp" ] # if blastp searches used
    then
        ${phome}/bin/Phage_Finder_v2.6.pl -t ncbi.out -i ${infofile} -r tRNAscan.out -n tmRNA_aragorn.out -c ncRNA_cmscan.out -A ${contigfile} -S
    elif [ "${s}" = "diamond" ] # if diamond searches used
    then
        ${phome}/bin/Phage_Finder_v2.6.pl -t diamond.out -i ${infofile} -r tRNAscan.out -n tmRNA_aragorn.out -c ncRNA_cmscan.out -A ${contigfile} -S
    fi

    ## generate linear figures of each prophage region

    echo "  generating linear illustrations..."
    ${phome}/bin/LinearDisplay_detailed.pl -A strict_dir/PFPR.frag -F strict_dir/PFPR.att -L -nt > strict_dir/PFPR_linear.fig
    ${fig2dev} -L pdf strict_dir/PFPR_linear.fig strict_dir/PFPR_linear.pdf
