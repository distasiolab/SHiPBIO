magma_dir=${HOME}/Programs/magma
trait=Amyloidosis_-_Zhou_30104761

mkdir -p out/step2
${magma_dir}/magma \
    --bfile ${magma_dir}/g1000_eur \
    --pval ${trait}.pval use='rsid,pval_SAIGE_NoSPA' ncol='n_total' \
    --gene-annot out/step1.genes.annot \
    --out out/step2/${trait}
