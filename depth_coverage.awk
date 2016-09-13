# set the chromosome to 1
BEGIN{chr=1}
{
    # save the total sum of depths and total counts of depth > 1 (ie coverage)
    totsum += $3;
    if($3>0){totcount++}

    # get the sum and count for the current chromosome and count the number of
    # positions (rows)
    if($1 == chr) {
      sum+=$3;
      row++;
      if($3>0){count++}
    }
    # if we have reached the end of the chromosome, print the statistics for the chromosome
    # then reset chr, sum, row and count.
    else {
      print "Chr" chr ":AvgDepth = " sum/row ", Coverage = " count/row
      chr = $1;
      sum = $3;
      row = 1;
      count = 0;
      if($3>0){count++}
    }
}
END {
    if (row > 0) {
        print "Chr" chr ":AvgDepth = " sum/row ", Coverage = " count/row;
    }
    if (NR > 0) {
        print "Total: AvgDepth=" totsum/NR ", Coverage = " totcount/NR;
    }
}