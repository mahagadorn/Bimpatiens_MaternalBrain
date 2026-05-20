#File for Merging individual volumetric data files into one table
#This is in lieu of copy and pasting data between files
examplesetup <- read.csv("../../Bimp_MaternalBrainVolumetric/TEST/K18B040_volume")

#Function for creating table
create_empty_table <- function(num_rows, num_cols, type_vec) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  for(i in 1:ncol(frame)) {
    print(type_vec[i])
    if(type_vec[i] == 'numeric') {frame[,i] <- as.numeric(frame[,i])}
    if(type_vec[i] == 'character') {frame[,i] <- as.character(frame[,i])}
    if(type_vec[i] == 'logical') {frame[,i] <- as.logical(frame[,i])}
    if(type_vec[i] == 'factor') {frame[,i] <- as.factor(frame[,i])}
  }
  return(frame)
}

#create empty table
#testfile <- create_empty_table(62,18, c("character", rep("numeric",17)))
volumetricData <- create_empty_table(65,18, c("character", rep("numeric",17)))

#Assign column names as objects from example file
colnames(volumetricData) <- c("sample_ID", examplesetup$Object)

#Get a list of all the files
individualfiles <- list.files(path = "../../Bimp_MaternalBrainVolumetric/StructureVolume/", pattern = "K18*", full.names = TRUE, recursive=FALSE)
#individual file names
library(qdapRegex)
names <- ex_between(individualfiles, "StructureVolume/", "_")


#read in files
LSdata<- lapply(individualfiles, function(i){read.csv2(i,header=TRUE, sep=",")})
 
#set index names in file
LSdata <- setNames(LSdata, names)
#Add row names
rownames(volumetricData) <- names #Works




volumetricData <- t(volumetricData)
volumetricData <- volumetricData[-1,]




#Adding files 

#K18B040
volumetricData[,"K18B040"] <- LSdata$K18B040$Volume[match(rownames(volumetricData), LSdata$K18B040$Object)]

#K18B061
volumetricData[,"K18B061"] <- LSdata$K18B061$Volume[match(rownames(volumetricData), LSdata$K18B061$Object)]

#K18B136
volumetricData[,"K18B136"] <- LSdata$K18B136$Volume[match(rownames(volumetricData), LSdata$K18B136$Object)]

#K18B139
volumetricData[,"K18B139"] <- LSdata$K18B139$Volume[match(rownames(volumetricData), LSdata$K18B139$Object)]

#K18B164 
volumetricData[,"K18B164"] <- LSdata$K18B164$Volume[match(rownames(volumetricData), LSdata$K18B164$Object)]

#K18B190
volumetricData[,"K18B190"] <- LSdata$K18B190$Volume[match(rownames(volumetricData), LSdata$K18B190$Object)]

#K18B191 
volumetricData[,"K18B191"] <- LSdata$K18B191$Volume[match(rownames(volumetricData), LSdata$K18B191$Object)]

#K18B236
volumetricData[,"K18B236"] <- LSdata$K18B236$Volume[match(rownames(volumetricData), LSdata$K18B236$Object)]

#K18C012 
volumetricData[,"K18C012"] <- LSdata$K18C012$Volume[match(rownames(volumetricData), LSdata$K18C012$Object)]

#K18C026 
volumetricData[,"K18C026"] <- LSdata$K18C026$Volume[match(rownames(volumetricData), LSdata$K18C026$Object)]

#K18C068
volumetricData[,"K18C068"] <- LSdata$K18C068$Volume[match(rownames(volumetricData), LSdata$K18C068$Object)]

#K18C119
volumetricData[,"K18C119"] <- LSdata$K18C119$Volume[match(rownames(volumetricData), LSdata$K18C119$Object)]

#K18C120
volumetricData[,"K18C120"] <- LSdata$K18C120$Volume[match(rownames(volumetricData), LSdata$K18C120$Object)]

#K18C129
volumetricData[,"K18C129"] <- LSdata$K18C129$Volume[match(rownames(volumetricData), LSdata$K18C129$Object)]

#K18C142
volumetricData[,"K18C142"] <- LSdata$K18C142$Volume[match(rownames(volumetricData), LSdata$K18C142$Object)]

#K18C143
volumetricData[,"K18C143"] <- LSdata$K18C143$Volume[match(rownames(volumetricData), LSdata$K18C143$Object)]

#K18C146
volumetricData[,"K18C146"] <- LSdata$K18C146$Volume[match(rownames(volumetricData), LSdata$K18C146$Object)]

#K18C149
volumetricData[,"K18C149"] <- LSdata$K18C149$Volume[match(rownames(volumetricData), LSdata$K18C149$Object)]

#K18C171
volumetricData[,"K18C171"] <- LSdata$K18C171$Volume[match(rownames(volumetricData), LSdata$K18C171$Object)]

#K18C172
volumetricData[,"K18C172"] <- LSdata$K18C172$Volume[match(rownames(volumetricData), LSdata$K18C172$Object)]

#K18C175
volumetricData[,"K18C175"] <- LSdata$K18C175$Volume[match(rownames(volumetricData), LSdata$K18C175$Object)]

#K18C222
volumetricData[,"K18C222"] <- LSdata$K18C222$Volume[match(rownames(volumetricData), LSdata$K18C222$Object)]

#K18C245
volumetricData[,"K18C245"] <- LSdata$K18C245$Volume[match(rownames(volumetricData), LSdata$K18C245$Object)]

#K18D014
volumetricData[,"K18D014"] <- LSdata$K18D014$Volume[match(rownames(volumetricData), LSdata$K18D014$Object)]

#K18D075
volumetricData[,"K18D075"] <- LSdata$K18D075$Volume[match(rownames(volumetricData), LSdata$K18D075$Object)]

#K18D084
volumetricData[,"K18D084"] <- LSdata$K18D084$Volume[match(rownames(volumetricData), LSdata$K18D084$Object)]

#K18D110
volumetricData[,"K18D110"] <- LSdata$K18D110$Volume[match(rownames(volumetricData), LSdata$K18D110$Object)]

#K18D127
volumetricData[,"K18D127"] <- LSdata$K18D127$Volume[match(rownames(volumetricData), LSdata$K18D127$Object)]

#K18D181
volumetricData[,"K18D181"] <- LSdata$K18D181$Volume[match(rownames(volumetricData), LSdata$K18D181$Object)]

#K18D213
volumetricData[,"K18D213"] <- LSdata$K18D213$Volume[match(rownames(volumetricData), LSdata$K18D213$Object)]

#K18D231
volumetricData[,"K18D231"] <- LSdata$K18D231$Volume[match(rownames(volumetricData), LSdata$K18D231$Object)]

#K18G409
volumetricData[,"K18G409"] <- LSdata$K18G409$Volume[match(rownames(volumetricData), LSdata$K18G409$Object)]

#K18G445
volumetricData[,"K18G445"] <- LSdata$K18G445$Volume[match(rownames(volumetricData), LSdata$K18G445$Object)]

#K18G447
volumetricData[,"K18G447"] <- LSdata$K18G447$Volume[match(rownames(volumetricData), LSdata$K18G447$Object)]

#K18G449
volumetricData[,"K18G449"] <- LSdata$K18G449$Volume[match(rownames(volumetricData), LSdata$K18G449$Object)]

#K18G462
volumetricData[,"K18G462"] <- LSdata$K18G462$Volume[match(rownames(volumetricData), LSdata$K18G462$Object)]

#K18G464
volumetricData[,"K18G464"] <- LSdata$K18G464$Volume[match(rownames(volumetricData), LSdata$K18G464$Object)]

#K18G474
volumetricData[,"K18G474"] <- LSdata$K18G474$Volume[match(rownames(volumetricData), LSdata$K18G474$Object)]

#K18G475
volumetricData[,"K18G475"] <- LSdata$K18G475$Volume[match(rownames(volumetricData), LSdata$K18G475$Object)]

#K18G530
volumetricData[,"K18G530"] <- LSdata$K18G530$Volume[match(rownames(volumetricData), LSdata$K18G530$Object)]

#K18G553
volumetricData[,"K18G553"] <- LSdata$K18G553$Volume[match(rownames(volumetricData), LSdata$K18G553$Object)]

#K18G556
volumetricData[,"K18G556"] <- LSdata$K18G556$Volume[match(rownames(volumetricData), LSdata$K18G556$Object)]

#K18G592
volumetricData[,"K18G592"] <- LSdata$K18G592$Volume[match(rownames(volumetricData), LSdata$K18G592$Object)]

#K18G620
volumetricData[,"K18G620"] <- LSdata$K18G620$Volume[match(rownames(volumetricData), LSdata$K18G620$Object)]

#K18H533
volumetricData[,"K18H533"] <- LSdata$K18H533$Volume[match(rownames(volumetricData), LSdata$K18H533$Object)]

#K18I411
volumetricData[,"K18I411"] <- LSdata$K18I411$Volume[match(rownames(volumetricData), LSdata$K18I411$Object)]

#K18I412
volumetricData[,"K18I412"] <- LSdata$K18I412$Volume[match(rownames(volumetricData), LSdata$K18I412$Object)]

#K18I431
volumetricData[,"K18I431"] <- LSdata$K18I431$Volume[match(rownames(volumetricData), LSdata$K18I431$Object)]

#K18I477
volumetricData[,"K18I477"] <- LSdata$K18I477$Volume[match(rownames(volumetricData), LSdata$K18I477$Object)]

#K18I481
volumetricData[,"K18I481"] <- LSdata$K18I481$Volume[match(rownames(volumetricData), LSdata$K18I481$Object)]

#K18I484
volumetricData[,"K18I484"] <- LSdata$K18I484$Volume[match(rownames(volumetricData), LSdata$K18I484$Object)]

#K18I513
volumetricData[,"K18I513"] <- LSdata$K18I513$Volume[match(rownames(volumetricData), LSdata$K18I513$Object)]

#K18I563
volumetricData[,"K18I563"] <- LSdata$K18I563$Volume[match(rownames(volumetricData), LSdata$K18I563$Object)]

#K18I565
volumetricData[,"K18I565"] <- LSdata$K18I565$Volume[match(rownames(volumetricData), LSdata$K18I565$Object)]

#K18I589
volumetricData[,"K18I589"] <- LSdata$K18I589$Volume[match(rownames(volumetricData), LSdata$K18I589$Object)]

#K18J419
volumetricData[,"K18J419"] <- LSdata$K18J419$Volume[match(rownames(volumetricData), LSdata$K18J419$Object)]

#K18J488
volumetricData[,"K18J488"] <- LSdata$K18J488$Volume[match(rownames(volumetricData), LSdata$K18J488$Object)]

#K18J490
volumetricData[,"K18J490"] <- LSdata$K18J490$Volume[match(rownames(volumetricData), LSdata$K18J490$Object)]

#K18J566
volumetricData[,"K18J566"] <- LSdata$K18J566$Volume[match(rownames(volumetricData), LSdata$K18J566$Object)]

#K18J569
volumetricData[,"K18J569"] <- LSdata$K18J569$Volume[match(rownames(volumetricData), LSdata$K18J569$Object)]

#K18J571
volumetricData[,"K18J571"] <- LSdata$K18J571$Volume[match(rownames(volumetricData), LSdata$K18J571$Object)]

#K18J572
volumetricData[,"K18J572"] <- LSdata$K18J572$Volume[match(rownames(volumetricData), LSdata$K18J572$Object)]

#K18I616--added late for similar sample sizes
volumetricData[,"K18I616"] <- LSdata$K18I616$Volume[match(rownames(volumetricData), LSdata$K18I616$Object)]

#K18G623--added late for similar sample sizes
volumetricData[,"K18G623"] <- LSdata$K18G623$Volume[match(rownames(volumetricData), LSdata$K18G623$Object)]

#K18J517--accidentally excluded the first time, traced with others though
volumetricData[,"K18J517"] <- LSdata$K18J517$Volume[match(rownames(volumetricData), LSdata$K18J517$Object)]


#Checking to see what 
  colnames(volumetricData)[colSums(is.na(volumetricData)) > 0]
  #character(0)
  
tvol <-t(volumetricData)
volumetricData <- as.data.frame(cbind(id=rownames(tvol), tvol))
  
write.table(volumetricData, "../../Bimp_MaternalBrainVolumetric/VolumetricData_scinotation.csv", sep = ",",row.names = FALSE, col.names = TRUE)








###Want to combine the volumetric data with TX data 

TXdata <- read.csv("G:/USUFiles_02072020/Documents/PHD/Studies/Bombus_impatiens/B.imp_MaternalBrain_CH2/Bimp_MomBrain_Treatment_Data_CSV.csv",
                   header = TRUE, sep=",")
#Remove extra columns and rows
TXdata <- TXdata[1:613,1:19]




#Subset Treatment data for only samples used
TXdata_samplesused <- subset(TXdata, id %in% volumetricData$id)
rownames(TXdata_samplesused) <- TXdata_samplesused$id

MB_volumetricdata <- merge(TXdata_samplesused, volumetricData)

write.table(MB_volumetricdata, "../../Bimp_MaternalBrainVolumetric/MomBrain_volumetricdata_Scinotation.csv", sep = ",",row.names = FALSE, col.names = TRUE)





