#!/usr/bin/env gawk -f
# Usage:
#   gawk -v key1=2 -v key2=1 -v valcol=3 -f join_genes_configurable.awk gProfiler.tsv zscore.tsv > output.tsv

BEGIN {
    FS = OFS = "\t";

    # Default columns if not passed
    if (!key1) key1 = 1;      # key in gProfiler
    if (!valcol) valcol = 3; # value to append from gProfiler
    if (!key2) key2 = 1;      # key in z-score file
}

# First file (gProfiler)
FNR == NR {
    if (NR == 1) next;  # skip header
    key = $key1;
    val = $valcol;
    map[key] = val;
    next;
}

# Second file (z-score)
FNR == 1 {
    print $0, "Gene_Name";  # header
    next;
}

{
    id = $key2;
    if (id in map) {
        print $0, map[id];
    }
}
