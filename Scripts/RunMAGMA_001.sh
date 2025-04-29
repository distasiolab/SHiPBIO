magma_dir=${HOME}/Programs/magma
trait=AMD_Fritsche-26691988

mkdir -p out/step1
${magma_dir}/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}/g1000_eur.bim \
    --gene-loc ${magma_dir}/NCBI37.3.gene.loc \
    --out out/step1

mkdir -p out/step2
${magma_dir}/magma \
    --bfile ${magma_dir}/g1000_eur \
    --pval ${trait}.pval use='SNP,GC_Pvalue' ncol='N' \
    --gene-annot out/step1.genes.annot \
    --out out/step2/${trait}
