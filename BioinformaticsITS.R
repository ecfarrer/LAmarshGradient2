#2017 ROOT and Soil fungal samples for Phragmites plots.

#First filter and infer sequence variants for each set: Root.Lane1,Root.Lane2, Soil
#For each set, run through fist part of script, CHANGE DIRECTORY folder

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/Roots_SoilFungi2017EF2.RData")

save.image(file = "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/Roots_SoilFungi2017EF2.RData")

#BiocManager::install(version = '3.10')

#BiocManager::install(("decontam"))

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.10")

packageVersion("dada2") #I'm using an old version b/c i didn't want to update r and everything

library(dada2)
library(ShortRead)
library(Biostrings)
library("phyloseq")
library("plyr")
library("ggplot2")
library(dplyr)
library(tidyr)
library(nlme)
library(vegan)
library(reshape)
library(stringr)
library(plotrix)
library(data.table)
library(decontam)
library(vegan)

#in unix add s to begining of file names
for i in *; do mv "$i" s"$i"; done;ls -l

#change path for each sequence set
path<- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/Soil"
path<- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/Roots1" #310 files
path<- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/Roots2" #148 files

list.files(path)

#separate forward and reverse
fnFs <- sort(list.files(path,pattern= "_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path,pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))


#Identify Primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA" 
REV <- "GCTGCGTTCTTCATCGATGC" 

#verify that we have the rich primers and correct orientation 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


#The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to "pre-filter" the sequences just to remove those with Ns, but perform no other filtering. I trimmed the ends of forward and reverse reads so that the quality score was roughly 30.
fnFs.filtN <- file.path(path,"filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path,"filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, trimRight = c(50,90))#soil is trimRight = c(50,90); roots1 is trimRight = c(50,95), roots2 is trimRight = c(50,90)

#count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. Identifying and counting the primers on one set of paired end FASTQ files is sufficient, assuming all the files were created using the same library preparation, so we'll just process the first sample.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#remove primers
cutadapt <- "/Users/farrer/miniconda3/envs/qiime2-2018.2/bin/cutadapt" 
list.files(cutadapt)
system2(cutadapt, args = "--version") # Run shell commands from R

#create output filenames for the cutadapt-ed files, and define the parameters we are going to give the cutadapt command
#critical parameters are the primers, and they need to be in the right orientation,
# i.e. the FWD primer should have been matching the forward-reads in its forward orientation, and the REV primer should have been matching the reverse-reads in its forward orientation. 
path.cut <- file.path(path, "cutadapt")
list.files(path.cut)
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt, 30 min
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-m", 20,# -n 2 required to remove FWD and REV from reads, -m means it deletes reads of less than 20 length
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#sanity check:count presence of primers in cut-adapeted files
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


plotQualityProfile(cutFs[1:8])+ 
  geom_vline(xintercept=250)+
  geom_hline(yintercept=30)

plotQualityProfile(cutRs[1:8])+
  geom_vline(xintercept=210)+
  geom_hline(yintercept=30)


#filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))



#Standard filtering parameters maxN=0 (already filtered), truncQ=2, rm.phix=TRUE, maxEE=TRUE.  soil is out2, roots1 is out3, roots2 is out4
#enforce min length of 50 bp
#start 10pm, end 10:06
out4 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out4)
#write.csv(out2, "readsout.csv")

#inspect read quality profiles 1,82, 92
plotQualityProfile(filtFs[70])
plotQualityProfile(filtRs[70])

#learn Error rates, start 10:20, end 10:22
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference, 4 min
#At this step the core sample inference algorithm is applied to the dereplicated data
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#merge paired reads, 2 min
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#construct sequence table
#seqtab.Root.Lane1, seqtab.Root.Lane2, seqtab.Soil
seqtab.Soil <- makeSequenceTable(mergers)

#Track reads through pipeline
getN <- function(x) sum(getUniques(x))
tracksoil<-cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) 
trackroots1<-cbind(out3, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
trackroots2<-cbind(out4, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))


#MERGE 3 sets################################################################################################
seqtab.all<-mergeSequenceTables(seqtab.Soil,seqtab.Roots1,seqtab.Roots2)

#remove chimeras, start 11:30, end 11:31
seqtab.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)

write.csv(seqtab.nochim, "seqtab.nochim.csv")
rownames(seqtab.nochim)

#inspect distribution of sequence length
plot(table(nchar(getSequences(seqtab.nochim)))) #229 is the most common length

#normally do track here but not if you merged multiple lanes
#track <- cbind(out6, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))




#assign taxonomy###############################################################################################

#Assign taxonomy, start 10:20am, end 1:50pm
unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAMarsh/Survey/Stats/Gradient/QIIME2/sh_general_release_dynamic_s_04.02.2020.fasta" #

taxa<-assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, minBoot=70, tryRC = TRUE,outputBootstraps=T)
taxaonly<-taxa$tax

write.csv(taxa,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/taxa.csv")

#create tax table for Phyloseq####################################################################################
sequences<-as.data.frame(rownames(taxaonly))
rownames(sequences)<-paste0("OTU",1:nrow(sequences))
sequences$OTU<-rownames(sequences)
colnames(sequences)<-c("sequence","OTUID")
OTUID<-as.data.frame(sequences$OTUID)

taxa1<-cbind(OTUID, taxaonly)
write.csv(taxa1,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/taxa_OTUID.csv")

OTUID<-sequences$OTUID
rownames(taxa1)<-OTUID

#use this for phyloseq tax table
taxa1<-taxa1[,2:8]

phy.tax<-tax_table(as.matrix(taxa1))
phy.tax

#Create sample data table###########################################################################################
sam<-as.data.frame(rownames(seqtab.nochim))
write.csv(sam, "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/sam.csv")
#clean in excel and reupload
sam1<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/sam.csv")
sam1<-sam1[,2:4]

# plot info
#soilchem<-read.csv("Soilchem.csv")

#join sam1 and rpi
#SD<-left_join(sam1,soilchem, by="Plot")
#colnames(SD)<-c("SampleID","Plot","SampleID.y","Site","SiteID"     ,"Transect" ,  "Calcium"   , "Copper" ,    "Magnesium" , "Phosphorus" ,"Potassium"  ,"Sodium"  ,   "Sulfur" ,   "Zinc" ,      "Carbon"  ,   "Nitrogen"  , "pH"    ,     "plant.rich" ,"sal" )

#root sample species
#species<-read.csv("Species_Samples.csv");colnames(species)<-c("SampleID","Species")
#SD<-left_join(SD,species,by="SampleID")

#add column in excel fpr decntam and clean
#write.csv(SD,"SampleData.csv")
#SD<-read.csv("SampleData.csv")
SD<-sam1
rownames(SD)<-SD$SampleID
SD$Sample_or_Control<-ifelse(SD$SampleType=="control","Control","Sample")
SD$is.neg<-ifelse(SD$Sample_or_Control=="Control",TRUE,FALSE)

#load in dat17, no i'll add it below on just the soil data
phy.SD<-sample_data(SD)

#Create official ASV Table####################################################################################
OTU.Table<-as.data.frame(seqtab.nochim)
colnames(OTU.Table)<-OTUID

#sam<-read.csv("sam.csv")
SampleID<-sam1$SampleID

#check
#cbind(rownames(OTU.Table),as.character(SampleID))
rownames(OTU.Table)<-SampleID
#OTU.Table$SampleID_Original<-rownames(OTU.Table)
#ASV.Table<-left_join(OTU.Table,sam,by="SampleID_Original")
#rownames(ASV.Table)<-ASV.Table$SampleID
#ASV.Table<-ASV.Table[,1:16456]

phy.ASV<-otu_table(OTU.Table,taxa_are_rows = F)

#create phyloseq object
phy<-merge_phyloseq(phy.ASV,phy.SD,phy.tax)




#### Decontam ####

datITS<-phy

print(as.data.frame(unique(tax_table(datITS)[,"Kingdom"])),row.names=F)

####Filter and remove contaminant sequences with decontam()####
#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement

#Filter samples with really low numbers of reads (<1000), not sure if this is necessary for the prevelance method but it is kind of weird that some of my actual samples had lower number of reads than the negative controls. This could be because the pcr failed and susannah put them in anyway "just in case" but if the pcr failed and all we got was contamination, then that skews the ability to see that they are contaminated.
sort(sample_sums(datITS))

datITSS<-datITS%>%
  subset_samples(sample_sums(datITS)>1000) 
datITSS

sample_data(datITSS)

dfITS <- as.data.frame(sample_data(datITSS)) # Put sample_data into a ggplot-friendly data.frame
dfITS$LibrarySize <- sample_sums(datITSS)
dfITS <- dfITS[order(dfITS$LibrarySize),]
dfITS$Index <- seq(nrow(dfITS))
ggplot(data=dfITS, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

contamdf.prevITS <- isContaminant(datITSS, method="prevalence", neg="is.neg")
table(contamdf.prevITS$contaminant)
contamdf.prevITS5 <- isContaminant(datITSS, method="prevalence", neg="is.neg",threshold = .2)
table(contamdf.prevITS5$contaminant)

#threshold .1: 51/53 contaminants, threshold .2: 65/ , threshold .5: 87/88 (taking out low # read samples/not taking them out) - this is an old comment from using carolyn's phy object

datITSS.pa <- transform_sample_counts(datITSS, function(abund) 1*(abund>0))
datITSS.pa.neg <- prune_samples(sample_data(datITSS.pa)$is.neg == TRUE, datITSS.pa)
datITSS.pa.pos <- prune_samples(sample_data(datITSS.pa)$is.neg == FALSE, datITSS.pa)

# Make data.frame of prevalence in positive and negative samples
#threshold .1
datITSSdf.pa <- data.frame(pa.pos=taxa_sums(datITSS.pa.pos), pa.neg=taxa_sums(datITSS.pa.neg),contaminant=contamdf.prevITS$contaminant)
ggplot(data=datITSSdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#threshold .2
datITSSdf.pa <- data.frame(pa.pos=taxa_sums(datITSS.pa.pos), pa.neg=taxa_sums(datITSS.pa.neg),contaminant=contamdf.prevITS5$contaminant)
ggplot(data=datITSSdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#I think I will use a threshold = 0.2. It captures the "contaminants" that were in all three negative controls but doesn't overdo it.

#take out contaminants and filter negative controls, there were 78 contaminants
datITSS2 <- prune_taxa(!contamdf.prevITS5$contaminant, datITSS)
datITSS3 <-datITSS2 %>%
  subset_samples(is.neg==FALSE)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T) #odd 84 taxa additionally are removed after removing the 78 contaminants. these must be taxa only found in the negative controls (which for some reason are not called contaminants) (??)

min(sample_sums(datITSS4))
sort(sample_sums(datITSS4))

#Filter out root samples and samples with low sampling depth
datITSS4 <-datITSS3 %>%
  subset_samples(SampleType=="soil")%>%
  subset_samples(sample_sums(datITSS3)>5000) %>%
  filter_taxa(function(x) sum(x) > (0), prune=T)


####### Rarefy to 5984 ######
#This takes out samples: 155, 118, 134, 135, 88, 99
datITSS4
datITSS5<-datITSS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS4)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
datITSS5c<-datITSS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS4)),rngseed=10,replace=F)
#2941 OTUs were removed because they are no longer present in any sample after random subsampling




#Adding to the sample data with dat from dat17
#need to open and run LoadingPlantSoilData.R script

tempf<-sample_data(datITSS5c)
tempf$Plot<-as.numeric(as.character(tempf$Plot))
tempd<-dat17[,c(1:4,52:77)]

tempb<-left_join(tempf,tempd,by="Plot")#warning is ok
tempb2<-cbind(Plottest=tempb$Plot,tempb[,6:34])

sample_data(datITSS5c)[,6:35]<-tempb2
sample_data(datITSS5c)[1:5,1:35]

#checks
cbind(as.character(sample_data(datITSS5c)$SampleID),as.numeric(as.character(sample_data(datITSS5c)$Plot)),sample_data(datITSS5c)$Plottest)
sample_data(datITSS5c)$Plottest<-NULL

sample_data(datITSS5)
sample_data(datITSS5)[,6:35]<-tempb2
cbind(as.numeric(as.character(sample_data(datITSS5)$Plot)),sample_data(datITSS5)$Plottest)
sample_data(datITSS5)$Plottest<-NULL


#Make OTU tables, this takes a while, start 10:47-10:52
datITSS5otu<-cbind(sample_data(datITSS5),otu_table(datITSS5))
datITSS5cotu<-cbind(sample_data(datITSS5c),otu_table(datITSS5c))




#Richness

richITS<-estimate_richness(datITSS5c, split = TRUE, measures = c("Observed", "Shannon","Chao1","Simpson","InvSimpson"))
colnames(richITS)<-c("RichnessITS","Chao1ITS","se.chao1ITS","ShannonITS","SimpsonITS","InvSimpsonITS")
richITS$SampleID<-rownames(richITS)
richITS2<-separate(richITS,SampleID,c(NA,"Plot"),"s")
richITS2$Plot<-as.numeric(richITS2$Plot)
#richITS2<-richITS%>%
#  left_join(data.frame(sample_data(datITS2rcsoil)))
#sample_data(datITSS5c)[,]

dat17b<-dat17%>%
  full_join(richITS2)
dat17<-dat17b




###### FUNguild ######
#https://github.com/UMNFuN/FUNGuild

temp<-data.frame(cbind(tax_table(datITSS5c),1))
temp$V8<-NULL
temp2<-temp%>%
  unite("Kingdom_Species",sep=";",remove=T)
head(temp2)
colnames(temp2)<-"taxonomy"

cbind(row.names(temp2),row.names(d)

datITSS5cfunguild<-cbind(t(otu_table(datITSS5c)),temp2)
datITSS5cfunguild[1:5,161:163]
write.csv(datITSS5cfunguild,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/datITSS5cfunguild.csv")
#then open the file and write "OTU ID" in the first cell

#then open terminal and navigate to folder and run
python Guilds_v1.1.py -otu datITSS5cfunguild.csv -m -u

#then open datITSS5cfunguild.guilds.txt and delete OTU ID in the first cell

funguild<-read.delim2("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/QIIME2/Fungi/analysisSept2020/datITSS5cfunguild.guilds.txt",header=T,row.names = 1)
head(funguild)
funguild[1:5,170:172]

#rarefied to 5984
colSums(funguild[,2:5])

funguild2<-funguild%>%
  filter(Trophic.Mode%in%c("Symbiotroph","Pathotroph","Saprotroph"))%>%
  arrange(Trophic.Mode)%>%
  dplyr::select(s1:s98,Trophic.Mode)%>%
  group_by(Trophic.Mode)%>%
  summarise_all(list(sum=sum))
funguild3<-data.frame(funguild2)
row.names(funguild3)<-funguild3$Trophic.Mode;funguild3$Trophic.Mode<-NULL
funguild4<-data.frame(t(funguild3))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1:3]<-funguild5[,1:3]/5984*100

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

#To get "Plant Pathogen" or "Endophyte" or "Arbuscular Mycorrhizal"
funguild2<-funguild%>%
  filter(Guild%in%c("Plant Pathogen","Endophyte","Arbuscular Mycorrhizal"))%>%
  arrange(Guild)%>%
  dplyr::select(s1:s98,Guild)%>%
  group_by(Guild)%>%
  summarise_all(list(sum=sum))
funguild3<-data.frame(funguild2)
row.names(funguild3)<-funguild3$Guild;funguild3$Guild<-NULL
funguild4<-data.frame(t(funguild3))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1:3]<-funguild5[,1:3]/5984*100

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

#Note: any taxon classified as one thing (i.e. plant pathogen with nothing else) always gets at least a probable, no possibles.
#Trying taking out the "probables" for plant pathogens so it only includes highly probable. if i do this, there are only 7 non zero plots for plant pathogens, so that's not reasonble 
funguild2<-funguild%>% 
  filter(Guild%in%c("Plant Pathogen","Endophyte","Arbuscular Mycorrhizal"))%>%
  filter(Confidence.Ranking=="Highly Probable")%>%
  arrange(Guild)%>%
  dplyr::select(s1:s98,Guild)%>%
  group_by(Guild)%>%
  summarise_all(list(sum=sum))
funguild3<-data.frame(funguild2)
row.names(funguild3)<-funguild3$Guild;funguild3$Guild<-NULL
funguild4<-data.frame(t(funguild3))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1:3]<-funguild5[,1:3]/5984*100
colnames(funguild5)[1:3]<-c("Arbuscular.Mycorrhizalhighlyprobable","Endophytehighlyprobable","Plant.Pathogenhighlyprobable")

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

#Getting data for ectomycorrhizal
funguild2<-funguild%>% 
  filter(Guild%in%c("Ectomycorrhizal"))%>%
  arrange(Guild)%>%
  dplyr::select(s1:s98,Guild)%>%
  group_by(Guild)%>%
  summarise_all(list(sum=sum))
funguild3<-data.frame(funguild2)
row.names(funguild3)<-funguild3$Guild;funguild3$Guild<-NULL
funguild4<-data.frame(t(funguild3))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1]<-funguild5[,1]/5984*100

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

dat17$PlantMutualist<-dat17$Ectomycorrhizal+dat17$Arbuscular.Mycorrhizal+dat17$Endophyte


#Try any description that includes "plant pathogen"
ind<-grep("Plant Pathogen",funguild$Guild)
funguild2<-funguild[ind,]
funguild3<-colSums(funguild2[,1:162])
funguild4<-data.frame(t(data.frame(t(funguild3))))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
colnames(funguild4)[1]<-"Plant.Pathogenanywhere"
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1]<-funguild5[,1]/5984*100

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b
