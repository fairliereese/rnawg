#!/bin/bash
#SBATCH --job-name=talon
#SBATCH -n 16
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=32G
#SBATCH --mail-user=freese@uci.edu

# usage
# sbatch sbatch_talon.sh -d <database> -n (if nanopore) <oprefix> <sample name>

ont=false

# get -n (nanopore) and -d (database) args
while getopts nd: flag; do
  case "$flag" in
    n) ont=true ;;
    d) db=${OPTARG} ;;
  esac
done

# shift flag args away and get opref / sample
shift $((OPTIND - 1))
opref=$1
sample=$2

ipref=$(dirname $(dirname $opref))/processing/$(basename $opref)
sam=${ipref}_merged_primers.sam
gtf=~/mortazavi_lab/ref/gencode.vM21/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf


# deal with different annotation names for the different tissues
if [ "$sample" == "adrenal" ]
  then
    build=mm10.fasta
  else
    build=mm10
fi

if "$ont"
  then
    printf "${sample},MinION,${sam}" > ${opref}_config.csv
  else
    printf "${sample},SequelII,${sam}" > ${opref}_config.csv
fi

if [ -z "$db" ]
  then
    echo "No database given. Will make new database."

    talon_initialize_database \
        --f ${gtf} \
        --g ${build} \
        --a gencode_vM21 \
        --l 0 \
        --idprefix ENCODEM \
        --5p 500 \
        --3p 300 \
        --o ${opref}

      db=${opref}.db

  else
    echo "Adding reads to database ${db}"
fi

if "$ont"
  then
    # use lower coverage requirements
    echo "Using ONT settings..."
      talon \
          --f ${opref}_config.csv \
          --cb \
          --db ${db} \
          --build $build \
          -t 32 \
          -c 0.8 \
          --o ${opref}
  else
    echo "Using PacBio settings..."
    talon \
          --f ${opref}_config.csv \
          --cb \
          --db ${db} \
          --build $build \
          -t 32 \
          --o ${opref}
fi
