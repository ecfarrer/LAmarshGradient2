#Genbank submission info and code


#### Info for bioproject #####
Bioproject Title:
  Marsh soil microbial communities 2017 Louisiana

We sampled soil microbial communities (ITS, 16S) in 8 marsh sites across a salinity gradient in SE Louisiana that were invaded by Phragmites australis. Sites ranged in salinity from freshwater marsh to saltmarsh. At each site, 21 plots were set up: 7 plots in a dense Phragmites stand, 7 plots in a transition zone ~50:50 Phragmites:native plants, and 7 plots in native-only areas. We tested whether soil microbial community compostition and diversity differed in invaded (Phragmites) vs. native areas. 

For Specimen type choose:
  Genome, metagenome or marker sequences (MIxS compliant)
Survey-related Marker Sequences (MIMARKS)
soil

LEQSF(2017–20)-RD-A-14	Louisiana Board of Regents	LA BOR	Global Change and Plant Invasions: Will Microbes Enhance Invasions in a Changing World?


##### Info for BioSamples mimarks files #####
#Biosamples are for example, one soil sample. Each biosample can be used in multiple sequencing data files (for this project we have 168 biosamples, and each biosample was used for 16S and ITS sequencing)

head(dat17)

latlon<-paste(format(dat17$Lat,digits=8,nsmall=8),"N",format(-dat17$Long,digits=8,nsmall=8), "W")
geoloc<-paste("USA:Louisiana",dat17$Site,sep=":")

ITSmimarks<-data.frame("sample_name"=paste(dat17$Plot,"_LAmarsh_Survey_2017_Soil",sep=""),organism="soil metagenome",collection_date=0,depth="10cm",elev="not collected",env_broad_scale="marsh",env_local_scale1=dat17$MarshClassV,env_local_scale2=dat17$Transect,env_medium="wetland soil",geo_loc_name=geoloc,lat_lon=latlon)

write.csv(ITSmimarks,"~/Dropbox/EmilyComputerBackup/Documents/Sequencing/GenBankSubmissions/ITSmimarks.csv",row.names = F)    



###### ITS #####

###### Changing file names in terminal #####

#previous submission I labeled samples: S.99.2015.ITS.R1.fastq.gz

#to count how many files in a folder in terminal
ls -1|wc -l

#sample names are:
#10_S10_L001_R1_001.fastq.gz

#I want to name it:
#9_LAmarsh_Survey_2017_Soil_ITS_R1.fastq.g

#I'm not doing the controls, so take them out of the folder before renaming files

# delete everything between _S and L001
for f in *; do mv "$f" "$(echo "$f" | sed 's/_S.*L001//')"; done;ls -l

#delete 001 at end
for f in *; do mv -n "$f" "${f/_001/}";done;ls -l

#add Soil_2017_ITS_ after the _
for f in *; do mv -n "$f" "${f/_/_Soil_2017_ITS_}";done;ls -l

###making some changes to naming
#delete everything after samplenumber at end
for f in *; do mv -n "$f" "${f/Soil_2017_ITS/}";done;ls -l

#add LAmarsh_Survey_2017_Soil_ITS after the _
for f in *; do mv -n "$f" "${f/_/_LAmarsh_Survey_2017_Soil_ITS}";done;ls -l


#making a list of file names for sra upload file

ITSSRAinfo<-data.frame(library=paste(dat17$Plot,"_LAmarsh_Survey_2017_Soil_ITS",sep=""),col1=paste(dat17$Plot,"_LAmarsh_Survey_2017_Soil_ITS_R1.fastq.gz",sep=""),col2=paste(dat17$Plot,"_LAmarsh_Survey_2017_Soil_ITS_R2.fastq.gz",sep=""))

write.csv(ITSSRAinfo,"~/Dropbox/EmilyComputerBackup/Documents/Sequencing/GenBankSubmissions/ITSSRAinfo.csv",row.names = F)    


#for uploading files with aspera
#after uploading the SRA file, click yes on the popup question about aspera
#choose FTP or Aspera Command Line file preload
#choose aspera command line uplod instructions
#click on "get the key file" and put it in the folder right outside the foler with the fastaq files
#follow the instructions for finding the path to the key file and the folder of your sequences and enter them below, then run in terminal

#commandline aspera
/Users/farrer/Applications/Aspera\ Connect.app/Contents/Resources/ascp -i /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Sequencing/GenBankSubmissions/aspera.openssh -QT -l100m -k1 -d /Users/farrer//Dropbox/EmilyComputerBackup/Documents/Sequencing/GenBankSubmissions/GenBankSRASubmissionITS subasp@upload.ncbi.nlm.nih.gov:uploads/efarrer_tulane.edu_chwkHMsu

#then click on the blue button that says "Select Preload Folder" and select your folder. make sure it has the correct number of files. If it has one more file than it should it might be fine b/c it uploaded a hidden file, on the next page that will be deleted.
#note: if the file upload doesn't finish for exmple due ot the internet cutting to or something, just run that same command again in terminal, it will skip files that are already uploaded.



###### 16S #####

# delete everything between _S and L001
for f in *; do mv "$f" "$(echo "$f" | sed 's/_S.*L001//')"; done;ls -l

#delete 001 at end
for f in *; do mv -n "$f" "${f/_001/}";done;ls -l

#add LAmarsh_Survey_2017_Soil_16S after the _
for f in *; do mv -n "$f" "${f/_/_LAmarsh_Survey_2017_Soil_16S_}";done;ls -l


#making a list of file names for sra upload file

SRAinfo16S<-data.frame(library=paste(dat17$Plot,"_LAmarsh_Survey_2017_Soil_16S",sep=""),col1=paste(dat17$Plot,"_LAmarsh_Survey_2017_Soil_16S_R1.fastq.gz",sep=""),col2=paste(dat17$Plot,"_LAmarsh_Survey_2017_Soil_16S_R2.fastq.gz",sep=""))

write.csv(SRAinfo16S,"~/Dropbox/EmilyComputerBackup/Documents/Sequencing/GenBankSubmissions/GenBankSRASubmissionSoil16S_5749_19062302/SRAinfo16S.csv",row.names = F)    

#commandline aspera
/Users/farrer/Applications/Aspera\ Connect.app/Contents/Resources/ascp -i /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Sequencing/GenBankSubmissions/aspera.openssh -QT -l100m -k1 -d /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Sequencing/GenBankSubmissions/GenBankSRASubmissionSoil16S_5749_19062302/GenBankSRASubmission16S subasp@upload.ncbi.nlm.nih.gov:uploads/efarrer_tulane.edu_chwkHMsu


