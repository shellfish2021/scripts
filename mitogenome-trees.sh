##此脚本最初由孔老师编写
##2021-12-9 魏靖乘进行了修改与补充并进行了第一次上传
##2022-12-6 时隔几乎一年，魏靖乘进行了第二次上传，增加了从转录组获取线粒体基因组的步骤并增加了两个一键完成式脚本：AUTOMITONC和AUTOMITOAA
##该脚本用于线粒体基因组系统发育树的构建，希望读者看完后能对线粒体基因组的系统发育分析流程有一定了解，并在一定程度上能合理使用这些方法来解决科学问题。--魏靖乘 2021.12.9

###检查下载序列的质量

md5sum -c MD5.txt

#blast
for species in contigs
do
	makeblastdb -in ${species}.fasta -dbtype nucl -out ${species}.db
	tblastn -query Circularized_assembly_1_H8.fasta -db ${species}.db -num_alignments 20 -out ${species}.blastn -evalue 1e-5
	perl /home/weijc/scripts/blast2table2.pl -top -format 10 -expect 0.00001 -index H ${species}.blastn > ${species}.blast2table
	cat ${species}.blast2table | awk -F "\t" '{print $12}' | sort | uniq > hit.txt 
	perl /home/weijc/scripts/select_contigs.pl -n hit.txt ${species}.fasta ${species}_potential_mito_genome.fasta
done

##MitoZ拼线粒体基因组（无需seed）

source activate
conda activate mitozEnv

nohup python3 MitoZ.py all --genetic_code 5 --clade Arthropoda --outprefix H4 \
--thread_number 12 \
--fastq1 /home/weijc/genome/H4_1.clean.fq.gz \
--fastq2 /home/weijc/genome/H4_2.clean.fq.gz \
--fastq_read_length 150 \
--insert_size 250  \
--run_mode 2 \
--filter_taxa_method 1 \
--requiring_taxa 'Arthropoda' &


##SPADES拼线粒体基因组

spades.py -m 32 -t 12 -1 /home/weijc/genome/H4_1.clean.fq.gz -2 /home/weijc/genome/H4_2.clean.fq.gz -o H4

##用得到的contigs建库

 makeblastdb -in /home/weijc/SPADES/H9/H9/contigs.fasta -dbtype nucl -out H9

 ##使用近缘物种的线粒体基因组比对

blastn -query Circularized_assembly_1_H23.fasta -db H9 -out dendromitogenom.fasta -task dc-megablast

##得到比对的序列，生成可能的线粒体基因组的fasta文件

 blast2table2.pl -format 9 -verbose dendromitogenom.fasta > blastout.pharsed
 cat blastout.pharsed | awk -F "\t" '{print $11}' | sort | uniq > lure.txt
select_contigs.pl -n lure.txt contigs.fasta dendroselect.fasta

##从转录组中获取线粒体基因组、

#kingfisher下载数据

conda create -c conda-forge -c bioconda -n kingfisher pigz python extern curl sra-tools pandas requests aria2  ##配置环境变量
conda activate kingfisher
nohup kingfisher get -r SRR3726700 -m aws-http &    

##Trinity拼接转录组

nohup Trinity --seqType fq --left SRR3726700_1.fastq.gz --right SRR3726700_2.fastq.gz --CPU 6  --max_memory 20G  --full_cleanup --bflyGCThreads 3 --bflyHeapSpaceMax 20G --output ./Noumeaella_rubrofasciata_trinity_out_dir &

nohup Trinity --seqType fq \
--left  SRR8573928_2.fastq.gz \
--right  SRR8573928_2.fastq.gz \
--CPU 20 --max_memory 20G --full_cleanup --bflyGCThreads 3 \
--bflyHeapSpaceMax 20G \
--trimmomatic &

##tblastx比对线粒体基因组

for file in *.Trinity.fasta
do
makeblastdb -in $file -out $file -dbtype nucl 
tblastx -query *.seed.fasta -db $file -num_alignments 20 -out $file.blastx -evalue 1e-5
perl /home/weijc/scripts/blast2table2.pl -top -format 10 -expect 0.00001 -index H -header $file.blastx > $file.blast2table 
cat $file.blast2table | awk -F "\t" '{print $12}' | sort | uniq > hit.txt 
perl /home/weijc/scripts/select_contigs.pl -n hit.txt $file potential_mito_genome.fasta
done

for file in *.Trinity.fasta
do
makeblastdb -in $file -out $file -dbtype nucl -parse_seqids
blastn -query *.seed.fasta -db $file -out $file.blastx -num_threads 8 -evalue 1e-5
perl /home/weijc/scripts/blast2table2.pl -format 9 -verbose $file.blastx > $file.blast2table 
cat $file.blast2table | awk -F "\t" '{print $11}' | sort | uniq > hit.txt 
perl /home/weijc/scripts/select_contigs.pl -n hit.txt $file potential_mito_genome.fasta
done

###用bwa回贴reads

bwa index -p huitie Noumeaella_rubrofasciata.nent.fasta                                           ##用拼出来的contigs建立数据库
bwa mem -t 12 huitie ../SRR3726700_1.fastq.gz ../SRR3726700_2.fastq.gz -o huitie.sam              ##把原始reads往数据库里的contigs回贴
samtools view -@ 10 huitie.sam -b > huitie.bam                                                     
samtools sort -@ 10 huitie.bam -o sorted_huitie.bam
rm huitie.sam huitie.bam
samtools index -@ 10 sorted_huitie.bam sorted_huitie.bam.bai
bedtools bamtobed -i ./sorted_huitie.bam > sorted_huitie.bed
bamdst -p sorted_huitie.bed -o ./ sorted_huitie.bam                                               ##获得回贴深度文件（检验是否是核基因）


##bowtie2获得线粒体基因组

mkdir bowtie_base    #创建bowtie_base文件夹 
mkdir bowtie         #创建bowtie文件夹
bowtie2-build Littorina_saxatilis.fasta Mito-Mollusca                                    #在bowtie_base文件夹创建线粒体数据库 
bowtie2 -k 1 -t -p 4 --al-conc ./matches_mito.fq --un-conc ./$hancock_Cleaned.fq \
  -x /home/weijc/mitofromtrans/Hancockia/Mito-Mollusca -1 ./*rRNAcleaned.1.fq -2 ./*rRNAcleaned.2.fq  #过滤掉线粒体基因组







rename 's/.txt/.fas/g' *.txt
#rename .txt .fas *.txt

Artemis格式
替换 (>.*?) .*$ 为 “\1”


替换>lcl|DQ238598.1_prot_ABB17317.1_1 [gene=COX1] [protein=cytochrome c oxidase subunit I] [protein_id=ABB17317.1] [location=1..1536] [gbkey=CDS]为COX1
替换 (>.*?) (.*?) .*$ 为 “>\2”
[gene=]

##为每个序列开头加上文件名

for FILENAME in *.fas
do
	ORTHOLOGY_GROUP=`echo $FILENAME | cut -d . -f 1`
	sed -i "s/>/>$ORTHOLOGY_GROUP"_"/g" $FILENAME
done

##将间断的序列整合在一起
#需要安装HaMStR-master
for FILENAME in *.fas
do
perl /home/weijc/scripts/nentferner.pl -in=$FILENAME -out=$FILENAME.nent 
done

mkdir ML_tree
mv *.fas.nent ML_tree/
cd ML_tree

rename 's/.fas.nent/.fas/g' *.fas.nent        ##这种语法适用于Perl语言环境 
rename .fas.nent .fas *.fas.nent              ##这种语法适用于C语言环境（实验室服务器上用这个）

##统一基因名称

sed -i 's/COX3/cox3/g' *.fas
sed -i 's/ND3/nad3/g' *.fas
sed -i 's/ND2/nad2/g' *.fas
sed -i 's/COX1/cox1/g' *.fas
sed -i 's/COX2/cox2/g' *.fas
sed -i 's/ATP8/atp8/g' *.fas
sed -i 's/ATP6/atp6/g' *.fas
sed -i 's/ND1/nad1/g' *.fas
sed -i 's/ND6/nad6/g' *.fas
sed -i 's/CYTB/cob/g' *.fas
sed -i 's/ND4L/nad4l/g' *.fas
sed -i 's/ND4/nad4/g' *.fas
sed -i 's/ND5/nad5/g' *.fas

sed -i 's/NAD1/nad1/g' *.fas
sed -i 's/NAD2/nad2/g' *.fas
sed -i 's/NAD3/nad3/g' *.fas
sed -i 's/NAD4L/nad4l/g' *.fas
sed -i 's/NAD4/nad4/g' *.fas
sed -i 's/NAD5/nad5/g' *.fas
sed -i 's/NAD6/nad6/g' *.fas
sed -i 's/cytB/cob/g' *.fas
sed -i 's/COXIII/cox3/g' *.fas
sed -i 's/COXII/cox2/g' *.fas
sed -i 's/COXI/cox1/g' *.fas
sed -i 's/nad4L/nad4l/g' *.fas
sed -i 's/nad4l/nad7/g' *.fas

cat *.fas > All.fa

##提取每个基因的序列至单独的文件

grep -A 1 cox1 All.fa > cox1.fasta
grep -A 1 cox2 All.fa > cox2.fasta
grep -A 1 cox3 All.fa > cox3.fasta
grep -A 1 atp6 All.fa > atp6.fasta
grep -A 1 atp8 All.fa > atp8.fasta
grep -A 1 nad2 All.fa > nad2.fasta
grep -A 1 nad5 All.fa > nad5.fasta
grep -A 1 nad4 All.fa > nad4.fasta      
grep -A 1 nad7 All.fa > nad4l.fasta
grep -A 1 cob All.fa > cob.fasta
grep -A 1 nad6 All.fa > nad6.fasta
grep -A 1 nad1 All.fa > nad1.fasta
grep -A 1 nad3 All.fa > nad3.fasta

sed -i 's/nad7/nad4l/g' nad4l.fasta

sed -i '/--/d' *.fasta
sed -i '/^$/d' *.fasta


##将核酸序列转换成蛋白序列

for FILENAME in *.fasta
do
transeq -sequence $FILENAME -outseq $FILENAME.pro.fasta -table 5  
done

##用mafft对序列进行比对

for FILENAME in *pro.fasta
do
mafft --thread 10 --auto --localpair --maxiterate 1000 $FILENAME > $FILENAME.aln 
done

rm -f *.fasta
rename 's/.fasta.aln/.fasta/g' *.fasta.aln
rename .fasta.aln .fasta *.fasta.aln       ##两种语法的区别同上，此后不再赘述

##统一物种的基因名

sed -i 's/_atp6//g' *.fasta
sed -i 's/_cob//g' *.fasta
sed -i 's/_cox1//g' *.fasta
sed -i 's/_cox2//g' *.fasta
sed -i 's/_cox3//g' *.fasta
sed -i 's/_nad1//g' *.fasta
sed -i 's/_nad2//g' *.fasta
sed -i 's/_nad3//g' *.fasta
sed -i 's/_nad4l//g' *.fasta
sed -i 's/_nad4//g' *.fasta
sed -i 's/_nad5//g' *.fasta
sed -i 's/_nad6//g' *.fasta
sed -i 's/_atp8//g' *.fasta
sed -i 's/\_\[gene\=atp6\]//g' *.fasta
sed -i 's/\_\[gene\=atp8\]//g' *.fasta
sed -i 's/\_\[gene\=cob\]//g' *.fasta
sed -i 's/\_\[gene\=cox1\]//g' *.fasta
sed -i 's/\_\[gene\=cox2\]//g' *.fasta
sed -i 's/\_\[gene\=cox3\]//g' *.fasta
sed -i 's/\_\[gene\=nad1\]//g' *.fasta
sed -i 's/\_\[gene\=nad2\]//g' *.fasta
sed -i 's/\_\[gene\=nad3\]//g' *.fasta
sed -i 's/\_\[gene\=nad4l\]//g' *.fasta
sed -i 's/\_\[gene\=nad4\]//g' *.fasta
sed -i 's/\_\[gene\=nad5\]//g' *.fasta
sed -i 's/\_\[gene\=nad6\]//g' *.fasta


##用Gblocks对比对过的序列过滤
#Gblocks
echo "Trimming alignments in Gblocks..."
for FILENAME in *.fasta
do
sequences=`grep -c \> $FILENAME`
half_of_sequences=`expr $sequences / 2`
b1=`expr $half_of_sequences + 1` #Sets b1 as the floating point calculation (# of sequences)/2 + 1. If (# of sequences)/2 is not a whole number, the output is truncated so that is why 1 is added to each value.
b2=$b1 #Sets b2 = b1
b3=8 #default value
b4=2 #default is 10
b5=a #default is none
Gblocks $FILENAME -t=p -p=n -b1=$b1 -b2=$b2 -b3=$b3 -b4=$b4 -b5=$b5 # Executes Gblocks
done
echo Done
echo

#MAFFT 
for FILENAME in *.fasta-gb
do
mafft --thread 10 --auto --localpair --maxiterate 1000 $FILENAME > $FILENAME.aln 
done

mkdir tree
mv *.fasta-gb.aln ./tree
cd tree/
rename 's/.fasta-gb.aln/.fas/g' *.fasta-gb.aln
rename .fasta-gb.aln .fas *.fasta-gb.aln

##用FASconCAT生成数据矩阵
#FASconCAT_v1.0.pl
perl /home/weijc/scripts/FASconCAT_v1.0.pl    ##选项依次输入p、p、s

##使用PartitionFinder寻找最适模型

#创建python=2.7的虚拟环境
conda create -n ptFinder python=2.7
source activate
conda activate ptFinder

#下载并解压PartitionFinder
wget https://github.com/brettc/partitionfinder/archive/refs/tags/v2.1.1.tar.gz
tar -zxvf partitionfinder-2.1.1.tar.gz
cd partitionfinder-2.1.1

#安装需要的python包（可以先试运行看缺少的dependence是哪一些）
conda install -c conda-forge numpy 
conda install -c conda-forge pandas 
conda install -c conda-forge pytables 
conda install -c conda-forge pyparsing 
conda install -c conda-forge scipy 
conda install -c conda-forge sklearn-contrib-lightning 

#将数据矩阵和配置文件放在同一目录下运行脚本
python PartitionFinder.py /home/weijc/software/partitionfinder-2.1.1/nudibranchia

#使用partitionfinder的宽松聚类（rcluster）模型
#首先策略选择rlucster：[schemes]search = rcluster;

python PartitionFinder.py /home/weijc/software/partitionfinder-2.1.1/nudibranchia -r


##使用iqtree建树（-p和-spp区别在于分区是否使用edge-link模型，-b 重复自举 -B超快自举）

iqtree -s FcC_smatrix.phy -spp partition.nex -alrt 1000 -b 1000 -T AUTO

##使用ModelFinderPro自动搜索模型（这个partition.nex只包含分区信息没有模型信息）

iqtree -s FcC_smatrix.phy -p partition.nex -m MFP+MERGE -B 1000 -T AUTO

#iqtree的GHOST模型

iqtree -s FcC_smatrix.phy -m GTR+FO*H4 -B 1000 -T AUTO

##Phylobayes建树（使用CAT+GTR+G4模型，模型可自行选择，每10棵树取一次，运行1000次）

mpirun -np 10 pb_mpi -d FcC_smatrix.phy -cat -gtr -dgam 4 -x 10 10000 chain1

mpirun -np 10 pb_mpi -d FcC_smatrix.phy -cat -gtr -dgam 4 -x 10 10000 chain2


#拟合两条链的树，丢弃前100棵树，此后每10棵树取1棵

bpcomp -x 100 10 chain1 chain2
##bmge过滤


bmge -i nad6.fasta.pro.fasta -t AA -of nad6.bmge.fasta
bmge -i nad5.fasta.pro.fasta -t AA -of nad5.bmge.fasta
bmge -i nad4l.fasta.pro.fasta -t AA -of nad4l.bmge.fasta
bmge -i nad4.fasta.pro.fasta -t AA -of nad4.bmge.fasta
bmge -i nad3.fasta.pro.fasta -t AA -of nad3.bmge.fasta
bmge -i nad2.fasta.pro.fasta -t AA -of nad2.bmge.fasta
bmge -i nad1.fasta.pro.fasta -t AA -of nad1.bmge.fasta
bmge -i cox3.fasta.pro.fasta -t AA -of cox3.bmge.fasta
bmge -i cox2.fasta.pro.fasta -t AA -of cox2.bmge.fasta
bmge -i cox1.fasta.pro.fasta -t AA -of cox1.bmge.fasta
bmge -i cob.fasta.pro.fasta -t AA -of cob.bmge.fasta
bmge -i atp8.fasta.pro.fasta -t AA -of atp8.bmge.fasta
bmge -i atp6.fasta.pro.fasta -t AA -of atp6.bmge.fasta

mkdir tree
mv *.bmge.fasta ./tree
cd tree/
rename .bmge.fasta .fas *.bmge.fasta


##调用曙光服务器的公共samtools

singularity exec -B /public1:/public1 /public1/db/sif/trinityrnaseq.v2.14.0.simg samtools

##翻转序列

seqkit seq -t dna -p -r potential_mito.genome.fasta > 1.fasta


cp /home/weijc/zhushi/work/finish/H24*.cds .
cat *.cds > H24.fasta
rm *.cds

cp /home/weijc/zhushi/work/finish/H24*.pep .
cat *.pep > H24.pro.fasta
rm *.pep


