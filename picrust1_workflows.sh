#zj@zlab.ac.cn
#2020-07
#诺和致源单端测序结果产出Greengenes结果
#合并所有fna文件
cat *.fna > seq.fna
#conda安装独立安装qiime1和picrust1
conda create -n picrust1 -c bioconda -c conda-forge picrust
conda activate picrust1
conda create -n qiime1 python=2.7
conda install qiime
#qiime1环境
conda activate qiime1
#将GG138数据库移动到公共数据库db文件夹
#基于fna文件生成GG数据
echo "pick_otus:enable_rev_strand_match True"  >> $PWD/otu_picking_params_97.txt
echo "pick_otus:similarity 0.97" >> $PWD/otu_picking_params_97.txt
pick_closed_reference_otus.py -i $PWD/seq.fna -o $PWD/New/ -p $PWD/otu_picking_params_97.txt -r ~/db/gg_13_8_otus/rep_set/97_otus.fasta -t ~/db/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
#基于默认下载数据库
echo "pick_otus:enable_rev_strand_match True"  >> $PWD/otu_picking_params_97.txt
echo "pick_otus:similarity 0.97" >> $PWD/otu_picking_params_97.txt
pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/ucrC97/ -p $PWD/otu_picking_params_97.txt -r $PWD/gg_13_5_otus/rep_set/97_otus.fasta -t $PWD/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt

#制作符合picrust格式的表格
filter_otus_from_otu_table.py -i otu_table.biom -o closed_otu_table.biom --negate_ids_to_exclude -e ~/db/gg_13_8_otus/rep_set/97_otus.fasta
#去除0reads样本（可选）
filter_samples_from_otu_table.py -i closed_otu_table.biom -o closed_otu_table_filt.biom -n 1
#Picrust1环境
conda activate picrust1
#生成picrust 所需格式，选了去除必须用closed_otu_table_filt.biom
normalize_by_copy_number.py -i closed_otu_table_filt.biom -o normalized_otus.biom
#预测功能并生成txt格式文件
#predict_metagenomes.py -i normalized_otus.biom -o predicted_metagenomes.biom
predict_metagenomes.py -f -i normalized_otus.biom -o predicted_metagenomes.txt
#biom convert -i otu_table.biom -o new.txt --to-tsv
#按KEGG L3层级注释出结果
categorize_by_function.py -f -i predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o predicted_metagenomes.L3.txt
#Level 2层级
categorize_by_function.py -f -i predicted_metagenomes.biom -c KEGG_Pathways -l 2 -o predicted_metagenomes.L2.txt
#Level 1层级
categorize_by_function.py -f -i predicted_metagenomes.biom -c KEGG_Pathways -l 1 -o predicted_metagenomes.L1.txt
#注释出哪个OTU贡献代谢功能
metagenome_contributions.py -i normalized_otus.biom -l K00001,K00002,K00004 -o ko_metagenome_contributions.tab
#打包结果
tar czvf picrust.tar seq