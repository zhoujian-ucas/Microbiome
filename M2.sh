#!/bin/bash
### 使用 vsearch+usearch11+QIIME 处理双端测序数据
### 目标 OTU表和注释信息
 
## 目录
# temp/: 临时文件目录
# res/: 结果文件目录
# data/:双端测序数据集
# database/:参考数据库
# err/: 错误信息
# log/: 操作日志
 
## 文件
HOME=$HOME/Lab_2		### 操作空间
temp=${HOME}/temp2	   	### 临时文件存放目录
res=${HOME}/res2		### 结果文件存放目录
data=${HOME}/seq		### 测序数据存放位置
refdb=${HOME}/database		### 参考测序数据库
err=${HOME}/err2		### 错误信息
log=${HOME}/log2		### 日志信息
 
 
########################################################################################### 
# 实验准备                                                                                #
# 实验数据为已合并为fasta格式的双端测序序列                                               #
# 安装usearch11、vsearch和qiime                                                           #
# 下载Greengene最新数据库，320MB                                                          #
# wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz     #
# 解压数据包后大小3.4G                                                                    #
# tar xvzf gg_13_8_otus.tar.gz                                                            #
###########################################################################################
 
echo "实验开始"
cd ${HOME}
 
echo "双端测序数据合并"
#echo "Step_1.1：decompress paired reads sequences"
#gunzip ${data}/*/*.fastq.gz		### 每个样品存放于一个目录中
 
prefix='ERR7765'			### 样本序号前缀
### merged.fastq 合并后的FASTQ格式reads
### merged.fasta 合并后的FASTA格式reads
### 55 样本编号末尾起始值
### 69 样本编号末尾结束值
 
echo "Step_1.2:Merge paired reads"
for i in `seq 55 69`;do
./usearch11 -fastq_mergepairs ${data}/${prefix}${i}/${prefix}${i}_1.fastq -reverse ${data}/${prefix}${i}/${prefix}${i}_2.fastq -fastqout ${temp}/${i}.merged.fastq -fastaout ${temp}/${i}.merged.fasta
done
 
echo "Step_2:Quality filtering"
for i in `seq 55 69`;do
./usearch11 -fastq_filter ${temp}/${i}.merged.fastq -fastq_maxee 1.0 -fastaout ${temp}/${i}.filtered.fasta -threads 2
done
cat ${temp}/*.filtered.fasta > ${res}/merged.fasta
#cat ${temp}/*.filtered.fastq > ${res}/merged.fastq
echo ""
ls -lh ${res}/merged.fasta
echo ""
 
echo "Step_3:Dereplication"
## 方式1
./usearch11 -fastx_uniques ${res}/merged.fasta -fastaout ${res}/unique.fasta --minuniquesize 2 -threads 2 -sizeout
## 方式2
# vsearch --derep_fulllength ${res}/filtered_v.fasta --output ${res}/unique_v.fa --minuniquesize 2 --threads 10 --sizeout
 
echo "Step_4:Discard singletons"
./usearch11 -sortbysize ${res}/unique.fasta -fastaout ${res}/discard.fasta -minsize 2
 
### 聚类OTU ###
echo "Step_5:UPARSE-OTU"
## 采用步骤3中的方式1的执行结果
./usearch11 -cluster_otus ${res}/unique.fasta -otus ${res}/otus_u51.fasta -uparseout ${res}/uparse.1.txt -relabel Otu -threads 10
## 采用步骤3中的执行结果
./usearch11 -cluster_otus ${res}/discard.fasta -otus ${res}/otus_u52.fasta -uparseout ${res}/uparse.2.txt -relabel Otu -threads 10
## Denoise 降噪得到的OTU
./usearch11 -unoise3 ${res}/unique.fasta -zotus ${res}/zotus_51.fasta -threads 30
 
## 采用QIIME的方式
echo "Picking OTUs"
pick_otus.py -i ${res}/unique.fasta -o ${res} -m uclust -r ${refdb}/gg_13_8_otus/rep_set/97_otus.fasta --threads 30
echo "Picking a representative sequence set, one sequence from each OTU"
pick_rep_set.py -i ${res}/unique_otus.txt -f ${res}/unique.fasta -o ${res}/rep_set1.fasta
 
 
 
### 统计fasta文件中序列的数量 ###
### 查看OTU数量  (分别为三种方式聚类得到的OTU数量以及一个降噪后得到的OTU数量) ###
echo "Number of OTU:"
grep '>' -c ${res}/otus_u51.fasta
grep '>' -c ${res}/otus_u52.fasta
grep '>' -c ${res}/zotus_51.fasta
 
echo "Step_6:Construct OTU table"
### Input should be reads before quality filtering and before discarding low-abundance unique sequences
### Reads must have sample indenfiers in the labels.
### The search database is the OTU (or ZOTU) sequences, which must have valid OTU identifiers in the labels.
### suggest making OTU tables for both 97% OTUs and denoised sequences (ZOTUs).
./usearch11 -otutab ${res}/merged.fasta -otus ${res}/otus_u51.fasta -otutabout ${res}/otutab_u61.txt -mapout ${res}/map_61.txt -threads 15
./usearch11 -otutab ${res}/merged.fasta -zotus ${res}/zotus_51.fasta -otutabout ${res}/zotutab_61.txt -mapout ${res}/zmap_61.txt -threads 15
 
# 基于数据库去嵌合体
echo "Step_7:Chimera detection"
otusfile='otus_u52.fasta'	### 选择上文中产生的OTU table
reffile='gold.fa'		### 选择一个参考数据库
printf "OTU表名：%s ,参考数据库名：%s \n"${otusfile} ${reffile}
### 基于所选择的参考数据库比对去除已知序列的嵌合体 ###
./usearch11 -uchime2_ref ${res}/${otusfile} -db ${refdb}/${reffile} -chimeras ${res}/otus_chimeras.fasta -notmatched ${res}/otus_rdp.fasta -uchimeout ${res}/otus_rdp.uchime -strand plus -mode sensitive -threads 30
 
### 获取嵌合体的序列ID ###
echo "Get sequenceID of chimeric:"
grep '>' ${res}/otus_chimeras.fasta|sed 's/>//g' > ${res}/otus_chimeras.id
### 剔除嵌合体的序列 ###
filter_fasta.py -f ${res}/${otusfile} -o ${res}/otus_non_chimera.fasta -s ${res}/otus_chimeras.id -n
 
### 检查是否为预期的序列数量
echo "Check whether it is the expected number of sequences:"
grep '>' -c ${res}/otus_non_chimera.fasta
 
# 去除非细菌序列
echo "Step_8:Removal of non-bacterial sequences"
### 下载Greengene最新数据库，320MB
### wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
### 解压数据包后大小3.4G
### tar xvzf gg_13_8_otus.tar.gz
### 将OTU与97%相似聚类的代表性序列多序列比对，大约8min
time align_seqs.py -i ${res}/otus_non_chimera.fasta -t ${refdb}/gg_13_8_otus/rep_set_aligned/97_otus.fasta -o ${res}/aligned/
 
# 无法比对细菌的数量
echo "Number of incomparable bacteria:"
grep -c '>' ${res}/aligned/otus_non_chimera_failures.fasta 
 
# 获得不像细菌的OTU ID
echo "Obtaining OTU IDs that are not like bacteria:"
grep '>' ${res}/aligned/otus_non_chimera_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > ${res}/aligned/otus_non_chimera_failures.id
 
# 过滤非细菌序列
filter_fasta.py -f ${res}/otus_non_chimera.fasta -o ${res}/otus_rdp_align.fasta -s ${res}/aligned/otus_non_chimera_failures.id -n
 
# 查看我们现在还有多少OTU
echo "How much is left of OTUs?:"
grep '>' -c ${res}/otus_rdp_align.fasta
 
# 产生代表性序列和OTU表
echo "Step_9:Generate representative sequences and OTU tables"
# 重命名OTU，这就是最终版的代表性序列，即Reference(可选)
echo "Step_9.1:Rename OTU"
awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' ${res}/otus_rdp_align.fasta > ${res}/rep_seqs.fasta
 
# 生成OTU表
echo "Step_9.2:Generate OTU tables"
./usearch11 -usearch_global ${res}/merged.fasta -db ${res}/rep_seqs.fasta -otutabout ${res}/otu_table.txt -biomout ${res}/otu_table.biom -strand plus -id 0.97 -threads 30
#less ${res}/otu_table.txt
 
# OTU物种注释，表操作
### 物种注释
### 必须添加rdp_classifier.jar的路径到PATH中，否则可能报错
### 如果没有rdp_classifier,可从如下地址下载,下载后进行解压
### https://jaist.dl.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.2.zip
echo "Step_10:Species annotation"
#RDP_Classifier_Path=$HOME/miniconda3/share/rdp_classifier-2.2-1/rdp_classifier.jar
#echo "export RDP_JAR_PATH=${RDP_Classifier_Path}" >> $HOME/.bashrc
echo "export RDP_JAR_PATH=$HOME/database/rdp_classifier_2.2/rdp_classifier-2.2.jar" >> $HOME/.bashrc
source $HOME/.bashrc
assign_taxonomy.py -i ${res}/rep_seqs.fasta -r ${refdb}/gg_13_8_otus/rep_set/97_otus.fasta -t ${refdb}/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -m rdp -o ${res}
#echo "Align the representative sequence set"
#align_seqs.py  -i ${res}/rep_seqs.fasta -t ${refdb}/gg_13_8_otus/rep_set_aligned/core_set_aligned.fasta.imputed -o ${res}/pynast_aligned/
#echo "Filtering the alignment prior to tre building and removing positions which are all gaps, or not useful for phylogenetic inference"
#filter_alignment.py -i ${res}/pynast_aligned/rep_seqs_aligned.fasta -o ${res}/filtered_alignment/ --remove_outliers
#echo "Building a phylogenetic tree"
#make_phylogeny.py -i ${res}/filtered_alignment/rep_seqs_aligned_pfiltered.fasta -o ${res}/rep_phylo.tre
#echo "Building an OTU table"
#make_otu_table.py -i ${res}/unique_otus.txt -t ${refdb}/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -o ${res}/otu_table_non_chimeric.biom
 
### OTU表统计、格式转换、添加信息
echo "Step_11:OTU Table Statistics, Format Conversion, Adding Information"
### 文本OTU表转换为BIOM：方便操作
biom convert -i ${res}/otu_table.txt -o ${res}/otu_table.biom --table-type="OTU table" --to-json
### 添加物种信息至OTU表最后一列，命名为taxonomy
biom add-metadata -i ${res}/otu_table.biom --observation-metadata-fp ${res}/rep_seqs_tax_assignments.txt -o ${res}/otu_table_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
### 转换biom为txt格式，带有物种注释：人类可读
biom convert -i ${res}/otu_table_tax.biom -o ${res}/otu_table_tax.txt --to-tsv --header-key taxonomy
### 查看OTU表的基本信息：样品，OUT数量统计
biom summarize-table -i ${res}/otu_table_tax.biom -o ${res}/otu_table_tax.sum
#less ${res}/otu_table_tax.sum
 
# OTU表筛选
echo "Step_12:Screen of OTU table"
# 按样品数据量过滤：选择counts>3000的样品
filter_samples_from_otu_table.py -i ${res}/otu_table_tax.biom -o ${res}/otu_table2.biom -n 3000
# 查看过滤后结果：只有15个样品，1129个OTU
echo "The result of screen step1:"
biom summarize-table -i ${res}/otu_table2.biom
# 按OTU丰度过滤：选择相对丰度均值大于万分之一的OTU
filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i ${res}/otu_table2.biom -o ${res}/otu_table3.biom
# 查看过滤后结果：只有15个样品，444个OTU
echo "The result of screen step2:"
biom summarize-table -i ${res}/otu_table3.biom
# 转换最终biom格式OTU表为文本OTU表格
biom convert -i ${res}/otu_table3.biom -o ${res}/otu_table3.txt --table-type="OTU table" --to-tsv
# OTU表格式调整方便R读取
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' ${res}/otu_table3.txt
# 筛选最终OTU表中对应的OTU序列
filter_fasta.py -f ${res}/rep_seqs.fasta -b ${res}/otu_table3.biom -o ${res}/rep_seqs4.fasta
 
printf "最终的OTU表:%s\t最终的OTU序列:%s\n"${res}/otu_table_tax.txt ${res}/rep_seqs.fasta
 
# End
echo "接下来:进化树,Alpha,Beta多样性分析以及后续操作"
echo "End