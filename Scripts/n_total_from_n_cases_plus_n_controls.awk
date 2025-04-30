BEGIN {
    FS = OFS = "\t";
}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", $i);  # Clean header fields
        if ($i == "n_cases") c1 = i;
        if ($i == "n_controls") c2 = i;
    }
    if (!c1 || !c2) {
        print "Error: Column not found" > "/dev/stderr"; exit 1
    }
    print $0, "n_total";
}
NR > 1 {
    # Strip carriage return from $0, then split into fields
    gsub(/\r$/, "", $0);
    print $0, $(c1) + $(c2);
}
