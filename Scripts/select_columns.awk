# select_columns.awk
BEGIN {
    OFS = FS = "\t";
    want[1] = "rsid";
    want[2] = "n_cases";
    want[3] = "n_controls";
    want[3] = "n_total";
}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/^[ \r\n\t]+|[ \r\n\t]+$/, "", $i);  # trim whitespace
        hdr[$i] = i;
    }
    for (i in want) {
        if (!(want[i] in hdr)) {
            print "Missing column:", want[i] > "/dev/stderr";
            exit 1;
        }
    }
    print want[1], want[2], want[3];
}
NR > 1 {
    print $(hdr[want[1]]), $(hdr[want[2]]), $(hdr[want[3]]);
}




