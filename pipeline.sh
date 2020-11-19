inputvcffile=$1
inputnormalbamfile=$2
inputtumorbamfile=$3
echo $inputnormalbamfile
echo $inputtumorbamfile
path=$inputvcffile
path2=$inputnormalbamfile
path3=$inputtumorbamfile
directory=${path##*/}
newdirectory=${directory%.vcf.gz}
normbamdirectory=${path2##*/}
tumorbamdirectory=${path3##*/}
newnormbamdirectory=${normbamdirectory%.bam}
cd /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/
mkdir $directory
cp /home/biodata/TCGA_Tumor_Het/Lung_TCGA_mutect2/$newdirectory/$directory /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory
cd /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/
suffix=".vcf"
gunzip -c /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/$directory>$newdirectory$suffix
source activate bcftools
bcftools query /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/$newdirectory$suffix -f '%CHROM %POS [\t%FILTER] [\t%AD] [\t%REF] [\t%ALT]\n'>/home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/output_file.txt
source deactivate bcftools
cd /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory
awk '$1!="chrX" && $1!="chrY" && $1!="chrM"' output_file.txt > temp_autosomal_allelic_counts.txt
awk '{print $1,$2}' temp_autosomal_allelic_counts.txt > autosomal_allelic_counts_1,2.txt
cd /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory
suffix2=".sorted.bam"
ln -s $inputnormalbamfile
samtools index $normbamdirectory
echo comehere
ln -s $inputtumorbamfile
samtools index $tumorbamdirectory
echo comethere
cd /home/biodata/miniconda2/alleleCount/example/filtered
source activate allelecount
alleleCounter -l /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/autosomal_allelic_counts_1,2.txt -b /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/$normbamdirectory -o /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/newoutput.txt
echo thisline
alleleCounter -l /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/autosomal_allelic_counts_1,2.txt -b /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/$tumorbamdirectory -o /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/042fftumoroutput_new.txt
echo thatline
source deactivate alleleCount
cd /home/biodata/miniconda2/ascat/filtered
mkdir $directory
cd $directory
touch tumorBAF.txt
touch normalBAF.txt
touch tumorLogR.txt
touch normalLogR.txt
cd /home/p2010-217-gpfs/xinling/test
echo temp there 5
python Cal_ASCAT_input.py $directory
echo temp there 6
cd /home/biodata/miniconda2/ascat/filtered/$directory
source activate ascat
Rscript /home/p2010-217-gpfs/xinling/test/ascat.R
source deactivate ascat
cp /home/biodata/miniconda2/ascat/filtered/$directory/major_minor_copy_number.tsv /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory
cp /home/biodata/miniconda2/ascat/filtered/$directory/tumor_purity.tsv /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory
cd /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory
awk '$2!="S2"' major_minor_copy_number.tsv > S1_major_minor_copy_number.tsv
grep -E "PASS" output_file.txt > filteredvcf_output_file.txt
awk '$1!="chrX" && $1!="chrY" && $1!="chrM"' filteredvcf_output_file.txt > new_temp_autosomal_allelic_counts.txt
awk '{print $1,$2,$6}' new_temp_autosomal_allelic_counts.txt > new_temp_autosomal_allelic_counts_2.txt
awk '{print $1,$2,$8,$10}' new_temp_autosomal_allelic_counts.txt > mutation_detail.txt
tr ','  ' ' <new_temp_autosomal_allelic_counts_2.txt > final_autosomal_allelic_counts.txt
awk '{print $2}' tumor_purity.tsv > new_tumor_purity.tsv
sed -n 2p new_tumor_purity.tsv > S1_tumor_purity.tsv
touch PyClone_input.tsv
cd /home/p2010-217-gpfs/xinling/test
python combine_copynumber_and_allelecounts.py $directory
x=$(wc -l < "/home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/S1_major_minor_copy_number.tsv")
numline=${x%% *}
if [ "$numline" == "1" ]; then
    echo "ascat can't produce optimal solution because most SNPs are not germline heterozygous"
    exit 3
fi
cd /home/biodata/miniconda2
source activate PyClone
tumor_purity="$(awk 'NR==1{print $1}' /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/S1_tumor_purity.tsv)"
PyClone setup_analysis --in_files /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/PyClone_input.tsv --working_dir /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory --init_method 'connected' --density pyclone_binomial --num_iters 100000 --prior total_copy_number --tumour_contents $tumor_purity
PyClone run_analysis --config_file /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/config.yaml
PyClone build_table --config_file /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/config.yaml --out_file /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/loci.tsv --table_type loci --burnin 10000
PyClone build_table --config_file /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/config.yaml --out_file /home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/$directory/cluster.tsv --table_type cluster --burnin 10000
source deactivate PyClone

