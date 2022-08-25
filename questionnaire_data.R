# R code for "Modelling dynamical changes in brain states with resting-state fMRI: A comparison of methods"
# Author: Liesa Ravijts
#
# Behavioural data acquired in the MPI-Leipzig Mind-Brain-Body dataset:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VMJ6NV
# ses.01 = LEMON protocol (n=228)
# ses.02 = N&C protocol (n=199)
# 109 participants took part in both protocols (total n=318)
#
# Three mind wandering questionnaires were explored to see which one would be
# the most fitting for the current study
# See also Table 2 in https://www.nature.com/articles/sdata2018307
#
# + Selection of the subsample (below)


if(!is.null(dev.list())) dev.off() # clear plots
rm(list=ls()); cat("\014") # clear workspace and console

# set working directory
setwd("D:/Users/Liesa/OneDrive/Universiteit Gent/Statistical Data Analysis/THESIS/Code")
write <- FALSE # set to TRUE to write files to working directory


### New York Cognition Questionnaire (NYC-Q)
NYCQ <- read.table("Data/NYCQ.tsv",header=TRUE,na.string="n/a")
variable.names(NYCQ)
# participant_id
# NYCQ_post.ses.01_1 -> 31 (23 items content of thoughts + 8 items form of thoughts)
# NYCQ_post.ses.02_1 -> 31 (23 items content of thoughts + 8 items form of thoughts)
# NYCQ_post.tasks_1 -> 23  (23 items content of thoughts; tasks: ETS and CCPT)

# sort participant ids
NYCQ <- NYCQ[order(NYCQ$participant_id),]
length(unique(NYCQ$participant_id))
# scored on a 9-point Likert scale
NYCQ[,-1][NYCQ[,-1] < 1] <- NA
NYCQ[,-1][NYCQ[,-1] > 9] <- NA
all(NYCQ[,-1]%%1==0,na.rm=TRUE)

# check available data
idx <- which(apply(NYCQ[,2:32],1,function(x) !all(is.na(x))))
NYCQ.ses01.post.scan <- NYCQ[idx,1] # 228 subjects
idx <- which(apply(NYCQ[,33:63],1,function(x) !all(is.na(x))))
NYCQ.ses02.post.scan <- NYCQ[idx,1] # 188 subjects
idx <- which(apply(NYCQ[,64:86],1,function(x) !all(is.na(x))))
NYCQ.ses02.post.tasks <- NYCQ[idx,1]

#write.csv(NYCQ,"Data/NYCQ_cleaned.csv",row.names=FALSE,quote=FALSE)


### Short Version of the New York Cognition Questionnaire (Short-NYC-Q)
SNYCQ <- read.table("Data/SNYCQ.tsv",header=TRUE,fill=TRUE,na.string="n/a")
variable.names(SNYCQ)
# participant_id
# SNYCQ_pre.ses.02_positive/.../intrusive (12 items)
# SNYCQ_post.ses.02.run.01.acq.AP_positive/.../intrusive (12 items)
# SNYCQ_post.ses.02.run.01.acq.PA_positive/.../intrusive (12 items)
# SNYCQ_post.ses.02.run.02.acq.AP_positive/.../intrusive (12 items)
# SNYCQ_post.ses.02.run.02.acq.PA_positive/.../intrusive (12 items)
# SNYCQ_post.ses.02.task.ETS_positive/.../intrusive (12 items)
# participant_id.1
# SNYCQ_post.ses.02.task.CPTS_positive/.../intrusive (12 items) (CPTS = CCPT)

# second column with participant ids
idcol2 <- which(colnames(SNYCQ)=="participant_id.1") 
# two dataframes pasted together, split them
SNYCQ.1 <- SNYCQ[,1:idcol2-1]
SNYCQ.2 <- SNYCQ[!SNYCQ$participant_id.1=="",idcol2:ncol(SNYCQ)]
# and merge them again
colnames(SNYCQ.2)[1] <- "participant_id"
SNYCQ.new <- merge(SNYCQ.1,SNYCQ.2,by="participant_id",all=TRUE) # full outer

# sort participant ids
SNYCQ.new <- SNYCQ.new[order(SNYCQ.new$participant_id),]
length(unique(SNYCQ.new$participant_id))
# scored on a scale ranging from 0% to 100% with increments of 5%
SNYCQ.new[,-1][SNYCQ.new[,-1] < 0] <- NA
SNYCQ.new[,-1][SNYCQ.new[,-1] > 100] <- NA
all(SNYCQ.new[,-1]%%5==0,na.rm=TRUE) # pre.scan & post.ETS are different versions
all(SNYCQ.new[,14:25]%%5==0,na.rm=TRUE) # check for in.scan1

# check available data
idx <- which(apply(SNYCQ.new[,2:13],1,function(x) !all(is.na(x))))
SNYCQ.ses02.pre.scan <- SNYCQ.new[idx,1]
idx <- which(apply(SNYCQ.new[,14:25],1,function(x) !all(is.na(x))))
SNYCQ.ses02.in.scan1.AP.run01 <- SNYCQ.new[idx,1] # 175 subjects
idx <- which(apply(SNYCQ.new[,26:37],1,function(x) !all(is.na(x))))
SNYCQ.ses02.in.scan2.PA.run01 <- SNYCQ.new[idx,1] # 174 subjects
idx <- which(apply(SNYCQ.new[,38:49],1,function(x) !all(is.na(x))))
SNYCQ.ses02.in.scan3.AP.run02 <- SNYCQ.new[idx,1] # 174 subjects
idx <- which(apply(SNYCQ.new[,50:61],1,function(x) !all(is.na(x))))
SNYCQ.ses02.in.scan4.PA.run02 <- SNYCQ.new[idx,1] # 170 subjects
idx <- which(apply(SNYCQ.new[,62:73],1,function(x) !all(is.na(x))))
SNYCQ.ses02.post.ETS <- SNYCQ.new[idx,1]
idx <- which(apply(SNYCQ.new[,74:85],1,function(x) !all(is.na(x))))
SNYCQ.ses02.post.CPTS <- SNYCQ.new[idx,1]

if(write) write.csv(SNYCQ.new,"Data/SNYCQ_cleaned.csv",row.names=FALSE,quote=FALSE)


### Spontaneous and Deliberate Mind Wandering (S-D-MW)
SDMW <- read.table("Data/SDMW.tsv",header=TRUE,na.string="n/a")
variable.names(SDMW)
# participant_id
# S.D.MW_delib_mean (mean score of 4 items on deliberate (intentional) mind wandering)
# S.D.MW_spont_mean (mean score of 4 items on spontaneous (unintentional) mind wandering)

# sort participant ids
SDMW <- SDMW[order(SDMW$participant_id),]
length(unique(SDMW$participant_id))
# scored on a 5-point Likert scale (mean score reported)
SDMW[,-1][SDMW[,-1] < 1] <- NA
SDMW[,-1][SDMW[,-1] > 5] <- NA
all(SDMW[,-1]%%0.25==0,na.rm=TRUE) # divided by 3?

# check available data
idx <- which(apply(SDMW[2:3],1,function(x) !all(is.na(x))))
SDMW.ses02.post.scan <- SDMW[idx,1]

#write.csv(SDMW,"Data/SDMW_cleaned.csv",row.names=FALSE,quote=FALSE)


rm(idx,idcol2)

##########################################################################
# Selection of participants from N&C with age <= 35 years and
# for which all SNYCQ scores after first rsfMRI run (in.scan1) are available

# read participant info (from https://openneuro.org/datasets/ds000221/versions/1.0.0)
participant_info <- read.table("Data/participants.tsv",header=TRUE,na.string="n/a")
participant_info <- participant_info[order(participant_info$participant_id),]
length(unique(participant_info$participant_id))
table(participant_info$gender,useNA="ifany"); table(participant_info$age,useNA="ifany")

# adjust "missing" data (info from https://www.nature.com/articles/sdata2018308#Sec72)
participant_info[participant_info$participant_id=="sub-010270",]$age <- "20-25"
participant_info[participant_info$participant_id=="sub-010293",]$age <- "20-25"

# select participants
participant_info$age <- factor(participant_info$age,ordered=TRUE)
NC <- scan("sub_ids_N&C.txt",what="",sep="\n",quiet=TRUE) # file from sort_participants.m
selected <- participant_info$participant_id[
    participant_info$age <= "30-35" & # age cut-off at 35 years
    participant_info$participant_id %in% SNYCQ.ses02.in.scan1.AP.run01 &
    apply(SNYCQ.new[,14:25],1,function(x) all(!is.na(x))) & # only complete cases
    participant_info$participant_id %in% NC # participants should be in N&C
    ] # note: all available cases are completes

# gender and age of selected participants
selected_info <- participant_info[participant_info$participant_id %in% selected,]
print(table(selected_info$gender)); print(table(droplevels(selected_info$age)))

# participants removed after MRI preprocessing (see main_script.m)
removed <- c("sub-010138","sub-010021","sub-010038","sub-010064",
             "sub-010065","sub-010067","sub-010071","sub-010075",
             "sub-010078","sub-010080","sub-010094","sub-010144",
             "sub-010151","sub-010171","sub-010217","sub-010221")
# sub-010138: file error fmap phasediff
# sub-010021 > sub-010217: high motion subjects
# sub-010221: time series extraction
selected <- selected[-which(selected %in% removed)]
removed_info <- selected_info[which(selected_info$participant_id %in% removed),]

if(write) write.table(selected,"sub_ids_selected.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

##########################################################################
# SNYCQ descriptives (N=119)
print(sprintf("N = %s", length(selected)))

# get scores of selected participants (after subject removal)
selected_SNYCQ <- SNYCQ.new[SNYCQ.new$participant_id %in% selected,c(1,14:25)]
selected_SNYCQ[,-1] <- lapply(selected_SNYCQ[,-1],function(x) factor(x,levels=seq(0,100,5),ordered=TRUE))
colnames(selected_SNYCQ)[-1] <- c("Content: Positive","Content: Negative","Content: Future",
                                  "Content: Past","Content: Myself","Content: Others",
                                  "Content: Surrounding","Vigilance","Form: Images",
                                  "Form: Words","Form: Specific","Form: Intrusive")

# five-number summaries + bar charts (loop over SNYCQ items)
if(write) png("Figures/A.2_SNYCQ_barplots.png",width=1800,height=2375,res=300)
par(mfrow=c(4,3),cex.main=1.1,cex.lab=1.1,mgp=c(3,0.8,0)) # some graphical parameters
x_ticks <- c("0",rep("",3),"20",rep("",3),"40",rep("",3),"60",rep("",3),"80",rep("",3),"100") # tick labels
for(i in 2:13){
    print(colnames(selected_SNYCQ)[i]) # item name
    print(quantile(selected_SNYCQ[,i],type=1),max.levels=0) # five-number summary
    # bar chart
    barplot(table(selected_SNYCQ[,i]),main=colnames(selected_SNYCQ)[i],ylim=c(0,16),
            names.arg=x_ticks,las=2) # for grid
    grid(nx=NA,ny=NULL) # display grid behind the bars, only horizontal lines
    barplot(table(selected_SNYCQ[,i]),main=colnames(selected_SNYCQ)[i],ylim=c(0,16),
            names.arg=x_ticks,las=2,add=TRUE) # add bars on top
    title(xlab="Score (%)",line=2.4); title(ylab="Number of subjects",line=2.4) # axis titles
}
if(write) dev.off()

rm(i,x_ticks)

