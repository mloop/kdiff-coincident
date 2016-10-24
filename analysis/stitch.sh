# Purpose: stitch results together into one file

# Author: Matthew Shane Loop

output_file_format="results-*-*.txt"

# $file_count variable value must match 54*number of iterations
file_count=$(ls | wc -l)
if (( $file_count != $(expr 54001) )); then
    echo "Number of files not expected."
    exit 0
fi

# Checking that number of rows in each output file from simulate.R is correct.
# The only possible value should be 5
for f in $output_file_format; do
    lc=`wc -l < $f`
    {
    if (( $lc != 5 )); then
    	echo "Error in $f: line count = $lc"
		exit 0
	fi
	}
done

# Merge files
echo "merging the outputs ...."
touch results.txt
grep -v [1] < results-1-1.txt | cat > results.txt
for f in $output_file_format; do
    grep -v ["c"] < $f | cat >> results.txt
	rm $f
done
echo "Done"
date
