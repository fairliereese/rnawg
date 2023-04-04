### Combining CAGE/RAMPAGE test and train peaks, before defining the new test/train sets for this study
### And filtering out peaks based on their chromosomal location

data_dir="data_dir/"
input_dir="$data_dir"/TSS_test_train_files
output_dir="$data_dir"/all_chr

for train_file in "$input_dir"/train/*txt; do
    #echo $train_file
    experiment=$(basename $train_file)
    experiment=${experiment##Training.}
    experiment=${experiment%%.txt}
    # echo $experiment

    test_file="$input_dir"/test/Test."$experiment".txt
    # echo $test_file

    tmp_output_file="$output_dir"/"$experiment".1.tmp
    awk 'BEGIN {FS=OFS="\t"}{NF--; print}' $train_file > $tmp_output_file
    cat $test_file >> $tmp_output_file

    output_file="$output_dir"/Allchr."$experiment".txt
    sort -k1,1 -k2,2n $tmp_output_file | grep "chr[[:digit:]]\+\|chrX\|chrY" > $output_file
	awk 'BEGIN{FS=OFS="\t"}{str=sprintf("%s_%s_%s", $1, $2, $3); print $1, $2, $3, str, $11}' $output_file > $tmp_output_file
	cp $tmp_output_file $output_file
    rm $tmp_output_file

    ls $output_file
    echo "---------------------------"
done


