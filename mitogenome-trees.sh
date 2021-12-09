##此脚本最初由孔老师编写
##2021-12-9 魏靖乘进行了修改与补充并进行了第一次上传
##该脚本用于线粒体基因组系统发育树的构建，希望读者看完后能对线粒体基因组的系统发育分析流程有一定了解，并在一定程度上能合理使用这些方法来解决科学问题。--魏靖乘 2021.12.9

##处理Artrmis导出的fasta格式文件

Artemis格式
替换 (>.*?) .*$ 为 “\1”

##处理NCBI上下载的文件（有时不同的作者论文中上传的文件格式可能不太一样，这一步灵活处理）

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
perl /home/weijc/scripte/nentferner.pl -in=$FILENAME -out=$FILENAME.nent 
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

cat *.fas > All.fa

##提取每个基因的序列至单独的文件

grep -A 1 cox1 All.fa > cox1.fasta
grep -A 1 cox2 All.fa > cox2.fasta
grep -A 1 cox3 All.fa > cox3.fasta
grep -A 1 atp6 All.fa > atp6.fasta
grep -A 1 atp8 All.fa > atp8.fasta
grep -A 1 nad2 All.fa > nad2.fasta
grep -A 1 nad5 All.fa > nad5.fasta
grep -A 1 -w nad4 All.fa > nad4.fasta | grep -A 1 nad4$ All.fa >> nad4.fasta      ##只抓取带有nad4的行而过滤掉nad4l（同样带有nad4）的行
grep -A 1 nad4l All.fa > nad4l.fasta
grep -A 1 cob All.fa > cob.fasta
grep -A 1 nad6 All.fa > nad6.fasta
grep -A 1 nad1 All.fa > nad1.fasta
grep -A 1 nad3 All.fa > nad3.fasta

sed -i '/--/d' *.fasta
sed -i '/^$/d' *.fasta

##用mafft对序列进行比对

for FILENAME in *.fasta
do
mafft --thread 10 --auto --localpair --maxiterate 1000 $FILENAME > $FILENAME.aln 
done

rm -f *.fasta
rename 's/.fasta.aln/.fasta/g' *.fasta.aln
rename .fasta.aln .fasta *.fasta.aln       ##两种语法的区别同上，此后不再赘述

sed -i 's/_atp6//g' *.fasta
sed -i 's/_cob//g' *.fasta
sed -i 's/_cox1//g' *.fasta
sed -i 's/_cox2//g' *.fasta
sed -i 's/_cox3//g' *.fasta
sed -i 's/_nad1//g' *.fasta
sed -i 's/_nad2//g' *.fasta
sed -i 's/_nad3//g' *.fasta
sed -i 's/_nad4$//g' *.fasta
sed -i 's/_nad4l//g' *.fasta
sed -i 's/_nad5//g' *.fasta
sed -i 's/_nad6//g' *.fasta
sed -i 's/_atp8//g' *.fasta

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
Perl /home/weijc/scripts/FASconCAT_v1.pl    ##选项依次输入p、p、s

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


##使用iqtree建树

iqtree -s FcC_smatrix.phy -spp partition.nex -alrt 1000 -bb 1000 -nt AUTO

iqtree -s FcC_smatrix.phy -p partition.nex -m MFP+MERGE -B 1000 -T AUTO

##Phylobayes建树

mpirun -np 10 pb_mpi -d Fcc_smatrix.phy -cat -gtr -x 10 1000 chain1

##使用IQtree的GHOST模型（更注重进化异质性的模型）

 #最泛用版本（所有参数全部独立运算）

 iqtree -s FcC_smatrix.fas -m GTR+FO*H4