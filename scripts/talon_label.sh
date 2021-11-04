opref=$1
files=$2

# please deal with the damn line endings
dos2unix $files

n=`wc -l $files | cut -d' ' -f1`

d="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/"
sbatch --array=1-${n} ${d}sbatch_talon_label.sh $opref $files
# sbatch --array=1 ${d}sbatch_talon_label.sh $opref $files
