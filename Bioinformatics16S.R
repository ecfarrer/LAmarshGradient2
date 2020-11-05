#Microbial data

##### Bacteria #####

#In terminal

#updat to newest conda
conda update conda
conda install wget

#install newest qiime
wget https://data.qiime2.org/distro/core/qiime2-2020.2-py36-osx-conda.yml
conda env create -n qiime2-2020.2 --file qiime2-2020.2-py36-osx-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2020.2-py36-osx-conda.yml


source activate qiime2-2020.2

#start 4:44pm end 5:00pm
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path Farrer_5749_19062302PhragSoil16S2017 \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path demux-paired-end.qza

#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#forward primer: GTGYCAGCMGCCGCGGTAA
#reverse primer: GGACTACNVGGGTWTCTAAT
#reverseComplement(DNAString("GTGYCAGCMGCCGCGGTAA"))

#start 5:00pm, end 5:40pm, ran fine, no error
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux-paired-end.qza \
--p-cores 4 \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--p-adapter-f ATTAGAWACCCBNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose

#reverse compement of forward primer might be at the end of a few of the R2 reads
# gunzip -k 1_S1_L001_R1_001.fastq.gz
# grep "GTGYCAGCMGCCGCGGTAA" 1_S1_L001_R1_001.fastq
# grep "TTACCGCGGCKGCTGRCAC" 1_S1_L001_R1_001.fastq
# gunzip -k 168_S165_L001_R2_001.fastq.gz
# grep "ATTAGA.ACCC" 168_S165_L001_R2_001.fastq
# grep "GGACTAC" 168_S165_L001_R2_001.fastq

#summarize, this takes a few minutes
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv

#export
qiime tools export \
--input-path trimmed-seqs.qza \
--output-path exported-trimmed-seqs
#CHECK  R1 almost all got trimmed to 272 or 281 bp (odd b/c CO dada wre trimmed to 253 bp). for the 272 ones, the forward primer did not get removed. I thnk it was only allowing one adapter per read to be removed and it checked for the -adapter-f first, then the front-f. there is some way to do this at the same time but I was having trouble figuring out how. the -times is probably the solution but -n is used in R DADA2 which is not an option in QIIME2. they give an example on https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage  for using -times to look for linked adapters (but this requires all reads to be bound by a 5' and 3' adapter, which is not exactly what I want). So I thnk I will just do two rounds which should be fine.



#trying cut adapt twice in succession to get all primers off
#start 10:20pm, end 10:49pm, ran fine, no error
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux-paired-end.qza \
--p-cores 4 \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--o-trimmed-sequences trimmed-seqsfront.qza \
--verbose

#this appears to have worked. all cases I looked at it said nearly 100% of the sequences had primers trimmed.

#didn't export
qiime tools export \
--input-path trimmed-seqsfront.qza \
--output-path exported-trimmed-seqsfront

#start 10:50, end 11:15?
qiime cutadapt trim-paired \
--i-demultiplexed-sequences trimmed-seqsfront.qza \
--p-cores 4 \
--p-adapter-f ATTAGAWACCCBNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \
--o-trimmed-sequences trimmed-seqsfrontadapt.qza \
--verbose

qiime tools export \
--input-path trimmed-seqsfrontadapt.qza \
--output-path exported-trimmed-seqsfrontadapt

qiime demux summarize \
--i-data trimmed-seqsfrontadapt.qza \
--o-visualization trimmed-seqsfrontadapt.qzv

qiime tools view trimmed-seqsfrontadapt.qzv


#start 8:20pm
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-seqsfrontadapt.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--o-denoising-stats denoising-stats \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 232 \
--p-trunc-len-r 200


#didn't run
qiime tools export \
--input-path rep-seqsN.qza \
--output-path exported-rep-seqsN

#to print the line with the longest number of characters:
cat dna-sequences.fasta|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 460bp, and blasts to a fungus, odd

#26277 are 253 bp
awk '{print length}' dna-sequences.fasta | sort | uniq -c


#export table to get otu table
qiime tools export \
--input-path table.qza \
--output-path exported-table

#convert biom to otu table text file!
biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#delete the space and # from '#OTU ID' on line 2 
sed '2s/[ ]//' otu_table.txt | sed '2s/.//' > otu_table2.txt

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file mappingfile.tsv

qiime tools view table.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv


##### Create a phylogenetic tree #####

#do the alignment, mask (filter) positions that are highly variable and will add noise to a phylogenetic tree, make tree from masked alignment, make an unrooted and rooted tree
#start 11pm, 

#crashed when I used all threads and when I used 3
#--p-n-threads auto \

qiime phylogeny align-to-tree-mafft-fasttree \
--p-n-threads 3 \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#export
qiime tools export \
--input-path rooted-tree.qza \
--output-path exported-rooted-tree


##### Assign taxonomy #####

#using greengenes

##### Greengenes database ######
#I redid this, need to if you update qiime
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path gg_13_8_otus/rep_set/99_otus.fasta \
--output-path gg_13_8_otus_99_otus.qza

#Import taxonomy
#Need to duplicate and save the txt file as tsv
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.tsv \
--output-path gg_13_8_otus_99_taxonomy.qza

#start 11:16 pm, end 11:42 pm
#note: if you have paired reads, do not set --p-trunc-len. I added the min and max - my reads are between 232-420, 100 to 400 was used by the moving pictures tutorial which use the same primers
qiime feature-classifier extract-reads \
--i-sequences gg_13_8_otus_99_otus.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--p-min-length 100 \
--p-max-length 420 \
--o-reads ref-seqs_all_99_gg_13_8.qza

#export to see what it trimmed
#qiime tools export \
#ref-seqs_all_99_gg_13_8.qza \
#--output-dir exported-ref-seqs_all_99_gg_13_8

#cd exported-ref-seqs_all_99_gg_13_8
#awk '{print length}' dna-sequences.fasta | sort | uniq -c
#Most lengths are 253, there are some longer and shorter from about 251-257

##Train classifier
#start 11:42pm, end 11:49pm
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_gg_13_8.qza \
--i-reference-taxonomy gg_13_8_otus_99_taxonomy.qza \
--o-classifier all_EMB_gg_13_8_classifier.qza


#Assign taxonomy
#start 2:45pm, end 3:03pm
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_gg_13_8_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -1 \
--o-classification taxonomy_gg.qza

qiime tools export \
--input-path taxonomy_gg.qza \
--output-path exported-taxonomy_gg


#navigate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#didn't do this
#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy_gg.qza \
--o-visualization taxonomy_gg.qzv

qiime tools view taxonomy_gg.qzv

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_gg.qza \
--m-metadata-file mappingfile.tsv \
--o-visualization taxa-bar-plots_gg.qzv

qiime tools view taxa-bar-plots_gg.qzv



####Summary: 
#I redid the bioinformatics for evrything except the tree just to check that McKenzie and I were getting the same results. As far as I can tell, the restuls are identical. even though i used a different updated qiime version and I trained my own classifier (not sure where she got hers). Even the numbers describing the features were the same.
#I will use her files just on the off chance that something was different and I want to make sure the tree file links up ok.

### Reading in my files as a check ####
otufileBacE<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/exported-table/otu_table2.txt",header=T)
head(otufileBacE)
cbind(otufileBacE$X1,otufileBac$X1)



#### Reading in McKenzie's files ####

otufileBac<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/FinalfilesfromMcKenzie/exported-feature-table-2.csv",header=T)
head(otufileBac)
taxfileBac<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/FinalfilesfromMcKenzie/taxonomy-2.csv",header=T)
head(taxfileBac)
colnames(taxfileBac)[1]<-"OTUID"

otufileBac2<-otufileBac%>%
  full_join(taxfileBac)

otufileBac2$taxonomy<-otufileBac2$Taxon
otufileBac2$Taxon<-NULL
otufileBac2$Confidence<-NULL
head(otufileBac2)

write.table(otufileBac2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/FinalfilesfromMcKenzie/otu_table.txt",sep="\t",row.names = F)

#open otu_table in excel and add '#OTU ID' as first cell name
biom convert -i otu_table.txt -o otu_table.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

otuBac <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/FinalfilesfromMcKenzie/otu_table.biom",parseFunction = parse_taxonomy_default)
otu_table(otuBac)[1:5,1:5]
tax_table(otuBac)
#remove Xs from beginning of sample names
colnames(otu_table(otuBac)) <- gsub("X","s",colnames(otu_table(otuBac)))

#need to go to LoadingPlantSoilData.R and run it
tempa<-data.frame(SampleID=colnames(otu_table(otuBac)),Plot=as.numeric(gsub("s","",colnames(otu_table(otuBac)))))
tempd<-dat17[,c(1:4,52:77)]
mapBac<-left_join(tempa,tempd)

write.table(mapBac,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/FinalfilesfromMcKenzie/mapBac.txt",sep="\t",row.names = F)

mapBac2<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/FinalfilesfromMcKenzie/mapBac.txt")

treeBac<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Bacteria/FinalfilesfromMcKenzie/tree-2.nwk")

datBac<-merge_phyloseq(otuBac,mapBac2,treeBac)




##### Bacteria cleaning ####

datBac2<-datBac%>%
  subset_taxa(Rank1=="k__Archaea"|Rank1=="k__Bacteria")%>%
  subset_taxa(is.na(Rank3)==T|Rank3!="c__Chloroplast")%>%
  subset_taxa(is.na(Rank5)==T|Rank5!="f__mitochondria")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

sort(sample_sums(datBac2))

datBac3 <-datBac2 %>%
  subset_samples(sample_sums(datBac2)>9000) %>%
  filter_taxa(function(x) sum(x) > (0), prune=T)


#Rarefying with rel abun, to 9316
# this takes out samples 12, 10, 112, 151, 165, 11, 29 
datBac4<-datBac3%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBac3)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )

#Rarefying without rel abun
datBac4c<-datBac3%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBac3)),rngseed=10,replace=F)
# removed 26275 OTUs 

sample_data(datBac4c)$Transect<-factor(sample_data(datBac4c)$Transect,levels=c("Native","Transition","Phragmites"))
sample_data(datBac4c)$Site<-factor(sample_data(datBac4c)$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
sample_data(datBac4)$Transect<-factor(sample_data(datBac4)$Transect,levels=c("Native","Transition","Phragmites"))
sample_data(datBac4)$Site<-factor(sample_data(datBac4)$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))



#Make OTU tables, this takes a while, start 10:05pm, end 1:30am
datBac4otu<-cbind(sample_data(datBac4),t(otu_table(datBac4)))
datBac4cotu<-cbind(sample_data(datBac4c),t(otu_table(datBac4c)))



#Richness - bacteria

richbac<-estimate_richness(datBac4c, split = TRUE, measures = c("Observed", "Shannon", "Chao1","se.chao1","Simpson","InvSimpson"))
colnames(richbac)<-c("RichnessBac","Chao1Bac","se.chao1Bac","ShannonBac","SimpsonBac","InvSimpsonBac")
richbac$SampleID<-rownames(richbac)
richbac2<-separate(richbac,SampleID,c(NA,"Plot"),"s")
richbac2$Plot<-as.numeric(richbac2$Plot)
#richITS2<-richITS%>%
#  left_join(data.frame(sample_data(datITS2rcsoil)))
#sample_data(datITSS5c)[,]

dat17b<-dat17%>%
  full_join(richbac2)
dat17<-dat17b

#Taking out bacterial sample with huge number of reads

dat17a<-dat17
ind<-which(dat17a$Chao1Bac>9000)
dat17a$Chao1Bac[ind]<-NA
dat17a$RichnessBac[ind]<-NA
dat17a$ShannonBac[ind]<-NA
dat17a$SimpsonBac[ind]<-NA
dat17a$InvSimpsonBac[ind]<-NA



