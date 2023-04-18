# 2.2
# prepare directories and move into the right place
. ~/.bashrc
mamba activate ngs
rm -rf ~/ngs_pipeline
mkdir -p ~/ngs_pipeline
cd ~/ngs_pipeline
mkdir -p data meta results logs
mkdir -p data/trimmed_fastq/trimmed_data
mkdir -p data/untrimmed_fastq
cp ../ngs_downloads/*fastq.gz data/untrimmed_fastq/
cp ../ngs_downloads/chr22.genes.hg19.bed data/
cp ../ngs_downloads/hg19.fa.gz data/


# In order to ensure that the environment is called correctly without issues for assessment
# all programs installed above with mamba are prepended with "mamba run -n ngs"
# Where there is redirected output I also need to add --no-capture-output to avoid interfering with the output
# This should normally not be necessary in a fully controlled environment, however to ensure smooth running
# during assessment where I cannot fully control the environment I chose to add this extra code for safety

# run fastqc and move it to the fastqc_untrimmed_reads directory
cd ~/ngs_pipeline/data/untrimmed_fastq
mamba run -n ngs fastqc -t 4 *.fastq.gz
mkdir -p ~/ngs_pipeline/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_pipeline/results/fastqc_untrimmed_reads/
echo "please review html files in ~/ngs_pipeline/results/fastqc_untrimmed_reads/"
# identify any issues with the summary by doing grep except PASS
cd ~/ngs_pipeline/results/fastqc_untrimmed_reads/
for zip in *.zip;do mamba run -n ngs unzip -qo $zip;done
cat */summary.txt > ~/ngs_pipeline/logs/fastqc_untrimmed_summaries.txt
grep -v PASS ~/ngs_pipeline/logs/fastqc_untrimmed_summaries.txt

# run trimmomatic
mkdir -p ~/ngs_pipeline/data/trimmed_fastq/trimmed_data
cd ~/ngs_pipeline/data/untrimmed_fastq
# to avoid version issues use realpath to find the most recent trimmomatic adapter version in the miniconda packages
mamba run -n ngs trimmomatic PE -threads 4 -phred33 WES01_chr22m_R1.fastq.gz WES01_chr22m_R2.fastq.gz -baseout ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R ILLUMINACLIP:$(realpath ~/miniconda3/pkgs/trimmomatic-*/share/trimmomatic-*/adapters/Nextera* | head -n 1):2:30:10 TRAILING:25 MINLEN:50

# redo fastqc on trimmed data
echo "fastqc"
mamba run -n ngs fastqc -t 4 ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R_1P ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R_2P
mkdir -p ~/ngs_pipeline/results/fastqc_trimmed_reads
mv  ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R_*_fastqc* ~/ngs_pipeline/results/fastqc_trimmed_reads/
echo "please review html files in ~/ngs_pipeline/results/fastqc_trimmed_reads/"
cd  ~/ngs_pipeline/results/fastqc_trimmed_reads/
for zip in *.zip;do mamba run -n ngs unzip -qo $zip;done
cat */summary.txt > ~/ngs_pipeline/logs/fastqc_trimmed_summaries.txt
grep -v PASS ~/ngs_pipeline/logs/fastqc_trimmed_summaries.txt

#2.3
# First index the hg19.fa.gz file, then run bwa mem
mkdir ~/ngs_pipeline/data/reference -p
cd ~/ngs_pipeline/data/reference
#cp ~/ngs_downloads/hg19.fa.gz .
#mamba run -n ngs bwa index hg19.fa.gz
mkdir -p ~/ngs_pipeline/data/aligned_data
#mamba run -n ngs --no-capture-output  bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 ~/ngs_pipeline/data/reference/hg19.fa.gz ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R_1P ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R_2P > ~/ngs_pipeline/data/aligned_data/WES01_chr22m.sam

# Instead run bowtie2 
cp ~/ngs_downloads/hg19.zip .
mamba run -n ngs unzip hg19.zip
mamba run -n ngs --no-capture-output bowtie2 -x hg19 -1 ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R_1P -2 ~/ngs_pipeline/data/trimmed_fastq/trimmed_data/WES01_chr22m_trimmed_R_2P > ~/ngs_pipeline/data/aligned_data/WES01_chr22m.sam

# turn the sam file into a sorted bam file for next steps
echo "samtools"
cd ~/ngs_pipeline/data/aligned_data
mamba run -n ngs --no-capture-output samtools view -h -S -b WES01_chr22m.sam > WES01_chr22m.bam
mamba run -n ngs --no-capture-output samtools sort WES01_chr22m.bam > WES01_chr22m_sorted.bam
mamba run -n ngs samtools index WES01_chr22m_sorted.bam

# use picard to mark duplicates and filter
mamba run -n ngs picard MarkDuplicates I=WES01_chr22m_sorted.bam O=WES01_chr22m_sorted_marked.bam M=marked_dup_metrics.txt

mamba run -n ngs samtools index WES01_chr22m_sorted_marked.bam
mamba run -n ngs samtools view -F 1796 -q 20 -o WES01_chr22m_sorted_filtered.bam WES01_chr22m_sorted_marked.bam
mamba run -n ngs samtools index WES01_chr22m_sorted_filtered.bam

echo "flagstats pre-filtering"
mamba run -n ngs samtools flagstat WES01_chr22m.bam
echo "flagstats post-filtering"
mamba run -n ngs samtools flagstat WES01_chr22m_sorted_filtered.bam

mamba run -n ngs --no-capture-output samtools view WES01_chr22m_sorted_filtered.bam > WES01_chr22m_sorted_filtered.sam

echo "idxstats post-filtering"
mamba run -n ngs --no-capture-output samtools idxstats WES01_chr22m_sorted_filtered.bam | column -t | sort -k 4 -r | column -t

# to avoid version issues use realpath to find the most recent picard version in the miniconda packages
mamba run -n ngs java -jar $(readlink -f ~/miniconda3/pkgs/picard-*/share/picard-*/picard.jar | head -n1) CollectInsertSizeMetrics I=WES01_chr22m_sorted_filtered.bam O=insert_size_metrics.txt   H=insert_size_histogram.pdf M=0.5

mamba run -n ngs --no-capture-output bedtools coverage -a WES01_chr22m_sorted_filtered.bam -b ~/ngs_pipeline/data/chr22.genes.hg19.bed > coverageBed.txt

# 2.4

# index the hg19 fasta file
zcat ~/ngs_pipeline/data/hg19.fa.gz > ~/ngs_pipeline/data/reference/hg19.fa
mamba run -n ngs samtools faidx ~/ngs_pipeline/data/reference/hg19.fa

# filter by bed file, then index the filtered file
mamba run -n ngs --no-capture-output samtools view -b -h -L ~/ngs_pipeline/data/chr22.genes.hg19.bed  ~/ngs_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered.bam  > ~/ngs_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered_bed_filtered.bam
mamba run -n ngs samtools index ~/ngs_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered_bed_filtered.bam
# run the Freebayes variant calling
mamba run -n ngs freebayes --bam ~/ngs_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered_bed_filtered.bam --fasta-reference ~/ngs_pipeline/data/reference/hg19.fa --vcf ~/ngs_pipeline/results/WES01_chr22m.vcf
mamba run -n ngs bgzip ~/ngs_pipeline/results/WES01_chr22m.vcf
mamba run -n ngs tabix -p vcf ~/ngs_pipeline/results/WES01_chr22m.vcf.gz
# quality filter results, using my choice of filters
mamba run -n ngs --no-capture-output vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_pipeline/results/WES01_chr22m.vcf.gz > ~/ngs_pipeline/results/WES01_chr22m_filtered.vcf
# intersect results with the bed file, just to make sure. Testing shows no change as it's been filtered above already
mamba run -n ngs --no-capture-output bedtools intersect -header -wa -a ~/ngs_pipeline/results/WES01_chr22m_filtered.vcf -b ~/ngs_pipeline/data/chr22.genes.hg19.bed > ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.vcf
mamba run -n ngs bgzip ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.vcf
mamba run -n ngs tabix -p vcf ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz


#2.5
# run annotate using annovar
cd ~/annovar
./convert2annovar.pl -format vcf4 ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz > ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.avinput
./table_annovar.pl ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.avinput humandb/ -buildver hg19 -out ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

# prioritization to exonic not seen in dbSNP, first using annovar
cd ~/ngs_pipeline/results/
# first write the header line, then use awk to filter for values in column 6 and then column 11 and write to basic_prioritization_annovar.csv
head -n1 WES01_chr22m_filtered_chr22.hg19_multianno.csv > basic_prioritization.csv
awk -F',' '$6 ~ "exonic" { print $0 } ' WES01_chr22m_filtered_chr22.hg19_multianno.csv | awk -F',' '$11 ~ "." { print $0 }' >> basic_prioritization_annovar.csv 

# run annotate using snpEff
cd ~/snpEff
mamba run -n ngs --no-capture-output java -Xmx8g -jar snpEff.jar hg19 ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz >  ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.hg19_snpEff.vcf

# prioritization to exonic not seen in dbSNP, next using snpEff
mamba run -n ngs --no-capture-output java -jar SnpSift.jar filter -f ~/ngs_pipeline/results/WES01_chr22m_filtered_chr22.hg19_snpEff.vcf "! exists ID" > ~/ngs_pipeline/results//WES01_chr22m_filtered_chr22.hg19_snpEff_not_in_dbSnp.vcf
grep ^# ~/ngs_pipeline/results//WES01_chr22m_filtered_chr22.hg19_snpEff_not_in_dbSnp.vcf > ~/ngs_pipeline/results/basic_prioritization_snpEff.vcf
grep -i exon ~/ngs_pipeline/results//WES01_chr22m_filtered_chr22.hg19_snpEff_not_in_dbSnp.vcf >> ~/ngs_pipeline/results/basic_prioritization_snpEff.vcf
