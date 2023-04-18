# 2.1
cd ~
# clear everything and start afresh
rm miniconda3/ annovar/ snpEff/ -rf
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh -nc
# make the downloaded file executable and install
chmod u+x Miniconda3-py310_23.1.0-1-Linux-x86_64.sh
./Miniconda3-py310_23.1.0-1-Linux-x86_64.sh -b
# required to start up miniconda for the first time
eval "$(~/miniconda3/bin/conda shell.bash hook)"

# refresh the terminal with conda activate the environment and install mamba
. ~/.bashrc
conda install -n base -c conda-forge mamba -y
sleep 5
# adds info into bashrc, then source bashrc again
mamba init
sleep 2
. ~/.bashrc
# make everything executable in conda (avoids future problems)
chmod -R 777 ~/miniconda3/
# install bioconda etc. sometimes fails so in a retry loop
for i in {1..5}; do mamba create -n ngs -y -c conda-forge -c bioconda -c defaults -c astrobiomike bit && break || sleep 15; done
sleep 2
# activate environment
. ~/.bashrc
mamba activate ngs
sleep 2
# This sometimes crashes, re-trying up to 5 times to get this done. Please be patient.
for i in {1..5}; do mamba run -n ngs mamba install -c bioconda bowtie2 samtools bwa freebayes picard bedtools trimmomatic fastqc vcflib prodigal unzip -y && break || sleep 15; done
mamba run -n ngs mamba install -c anaconda fontconfig -y
mamba run -n ngs mamba install -c conda-forge r-essentials -y


# centralise all downloads into downloads folder for future use (save bandwidth on repeated runs with -nc option)
mkdir -p ~/ngs_downloads
cd ~/ngs_downloads
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -nc
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/WES01_chr22m_R1.fastq.gz -nc
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/WES01_chr22m_R2.fastq.gz -nc
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/chr22.genes.hg19.bed -nc
wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip -nc
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz -nc
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -nc

cd ~
cp ~/ngs_downloads/annovar.latest.tar.gz .
tar -zxf annovar.latest.tar.gz
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

cd ~
cp ~/ngs_downloads/snpEff_latest_core.zip .

mamba run -n ngs unzip -qo snpEff_latest_core.zip
cd snpEff
mamba run -n ngs java -jar snpEff.jar download -v hg19
mamba run -n ngs java -Xmx8g -jar snpEff.jar dump -v hg19 > /dev/null

