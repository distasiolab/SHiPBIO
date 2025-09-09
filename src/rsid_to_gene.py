
#!/usr/bin/env python3
import argparse
import cyvcf2
import pandas as pd
import os
from tqdm import tqdm

CACHE_FILE = "rsid_to_gene.tsv"


def build_rsid_to_gene(vcf_file, map_file, update_cache=True):
    """
    Build an rsID → gene symbol mapping from dbSNP VCF.
    Saves TSV with columns: rsid, gene.
    Shows tqdm progress bar with percentage.
    """
    rsid_list = []
    gene_list = []

    vcf = cyvcf2.VCF(vcf_file)

    # Try to get total number of records
    total = getattr(vcf, "num_records", None)
    if total is None:
        print("Counting total variants in VCF for progress bar...")
        total = sum(1 for _ in cyvcf2.VCF(vcf_file))
        print(f"Total variants: {total}")

    bar_format = "{l_bar}{bar}| {n_fmt}/{total_fmt} variants [{elapsed}<{remaining}, {rate_fmt}] {percentage:3.0f}%"

    for var in tqdm(vcf, total=total, unit="variants", desc="Processing VCF", bar_format=bar_format):
        rsid = var.ID
        geneinfo = var.INFO.get("GENEINFO")
        if rsid and geneinfo:
            gene_symbol = geneinfo.split(":")[0]
            rsid_list.append(rsid)
            gene_list.append(gene_symbol)

    df = pd.DataFrame({"rsid": rsid_list, "gene": gene_list})
    df.to_csv(map_file, sep="\t", index=False)
    print(f"\n✅ Saved mapping: {map_file} ({len(df)} SNPs with gene info)")

    # Save/update cache
    if update_cache and map_file != CACHE_FILE:
        df.to_csv(CACHE_FILE, sep="\t", index=False)
        print(f"💾 Cached mapping at {CACHE_FILE}")

    return df


def load_mapping(map_file=None):
    """
    Load a prebuilt rsID → gene mapping TSV.
    Defaults to local cache file if no map_file provided.
    """
    path = map_file if map_file else CACHE_FILE

    if not os.path.exists(path):
        raise FileNotFoundError(
            f"No mapping file found. Provide --vcf to build one or place cache at {CACHE_FILE}"
        )

    df = pd.read_csv(path, sep="\t")
    if not {"rsid", "gene"}.issubset(df.columns):
        raise ValueError("Mapping file must have columns: rsid, gene")

    print(f"✅ Loaded mapping from {path} ({len(df)} entries)")
    return df


def annotate_snp_file(snp_file, mapping_df, output_file):
    """
    Annotate an SNP result TSV with gene symbols using the mapping DataFrame.
    Assumes SNP column is named 'SNP'.
    """
    df = pd.read_csv(snp_file, sep="\t")

    if "SNP" not in df.columns:
        raise ValueError(f"Input SNP file must have a 'SNP' column. Found: {df.columns}")

    rs2gene = dict(zip(mapping_df["rsid"], mapping_df["gene"]))
    df["GENE"] = df["SNP"].map(rs2gene).fillna("NA")

    df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Annotated SNP file saved: {output_file}")


def merge_chromosome_maps(chrom_dir, output_file):
    """
    Merge all per-chromosome mapping TSVs in a directory into a single TSV.
    Assumes each file has columns: rsid, gene
    """
    import glob

    tsv_files = sorted(glob.glob(os.path.join(chrom_dir, "*.tsv")))
    if not tsv_files:
        raise FileNotFoundError(f"No TSV files found in {chrom_dir}")

    merged_df = pd.concat(
        [pd.read_csv(f, sep="\t") for f in tsv_files],
        ignore_index=True
    )

    merged_df.drop_duplicates(subset="rsid", inplace=True)
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Merged {len(tsv_files)} chromosome maps into {output_file} ({len(merged_df)} SNPs)")
    return merged_df

    

def main():
    parser = argparse.ArgumentParser(
        description="Build/load rsID → gene mapping and annotate SNP result files."
    )
    parser.add_argument("--vcf", help="Path to dbSNP VCF file (.gz) to build mapping")
    parser.add_argument("--map", help="Mapping TSV file (rsid,gene). If --vcf is given, this will be created/overwritten. Defaults to cache (rsid_to_gene.tsv).")
    parser.add_argument("--annotate", nargs="+", help="One or more SNP TSV files to annotate (must contain 'SNP' column)")
    parser.add_argument("--annotated-out", help="Output file for single annotated SNP TSV (ignored if multiple input files)")
    parser.add_argument("--force-cache", action="store_true", help="Force overwrite cached mapping")
    parser.add_argument("--merge-dir", help="Directory with per-chromosome mapping TSVs to merge into a single map")
    args = parser.parse_args()


    # Step 0: Merge chromosome maps if requested
    if args.merge_dir:
        map_file = args.map if args.map else CACHE_FILE
        mapping_df = merge_chromosome_maps(args.merge_dir, map_file)
    else:
        # Step 1: Build or load mapping
        if args.vcf:
            map_file = args.map if args.map else CACHE_FILE
            mapping_df = build_rsid_to_gene(args.vcf, map_file, update_cache=True)
        else:
            mapping_df = load_mapping(args.map)

    # Step 2: Annotate SNP file(s) if requested
    if args.annotate:
        if len(args.annotate) == 1:
            out_file = args.annotated_out if args.annotated_out else os.path.splitext(args.annotate[0])[0] + "_with_genes.tsv"
            annotate_snp_file(args.annotate[0], mapping_df, out_file)
        else:
            if args.annotated_out:
                print("⚠️  Ignoring --annotated-out because multiple files were provided.")
            for infile in args.annotate:
                base, _ = os.path.splitext(infile)
                out_file = f"{base}_with_genes.tsv"
                if not os.path.exists(out_file):
                    annotate_snp_file(infile, mapping_df, out_file)
                else:
                    print(f"Skipping {infile}, {out_file} already exists.")

if __name__ == "__main__":
    main()
