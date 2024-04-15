##### Lizard Virome Analysis ABBREVIATED #####
##### version 20240412 ######

#------------------------------------------------------------------------------------------------------#
#set up workspace

setwd("/Users/mwille/Desktop/UniMelb/2024_Jackie_Lizards/20240412")

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(multcomp)

#if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
#BiocManager::install("phyloseq")

library(vegan)
library(phyloseq)

#------------------------------------------------------------------------------------------------------#
#Notes: Some of the library names and taxonomy names are slightly different in the R script/dataset and the final manuscript:
#Articulavirales >> Amnoonviridae
#Flavivirus >> Orthoflavivirus
#Ccyg_M >> Cmet_M
#Gaus_M >> Garn_M

#------------------------------------------------------------------------------------------------------#
#read in data, basic data cleaning

dat <- read.csv("Lizard_viruses_for_alphaDiv_final.csv",na.strings=c("", "NA"), header=TRUE,sep=",")

#filter out unwanted coloumns i.e. suspected index hopping
unique(dat$Remove)
dat[is.na(dat)] <- "keep"
dat2 <-filter(dat, Remove== "keep")

#------------------------------------------------------------------------------------------------------#
# Overview of dataset

aggregate(expected_count ~ Library + Host_infraorder  + Host_family+Host_genus+Habitat, data=dat, FUN=sum)

# Overview Biologically Relevant vs Other Viruses

#Aggregate data
lizard.host<-aggregate(expected_count ~ Library + broad_host + Non.host.rRNA_reads, data=dat2, FUN=sum)
#Calculate "abundance" which is the proportion of total reads
lizard.host$Abundance<-(lizard.host$expected_count/lizard.host$Non.host.rRNA_reads)*100

#reorder for plotting
#https://r-graphics.org/recipe-dataprep-factor-reorder#RECIPE-DATAPREP-FACTOR-REORDER
lizard.host$broad_host2<-factor(lizard.host$broad_host, levels=c("Cam_M", "Cmun_M ", "Csex_M", "Ccyg_M",
                                                                 "Gaus_M", "Gnan_A", "Gnan_M",
                                                                 "Hbin_A", "Hbin_M", "Hplan_A",
                                                                 "Omar_M"))

p.overview<-ggplot(lizard.host, aes(x=Library, y=Abundance, fill=broad_host)) + 
  geom_bar(position="stack", stat="identity")  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = c("goldenrod", "#39993A"))+theme_classic()
ggsave("1_Vert_Non_Vert_20230725.pdf", height = 5, width = 7, units = c("in"), dpi = 600) 

#------------------------------------------------------------------------------------------------------#
# further cleaning

#Remove non-biologically relevant viruses from the dataset
lizards<-filter(dat2, broad_host=="Vertebrates")

# Remove Csex_M it’s the only library found on riparian habitat
# Remove Diplodactylidae as its the sole member of its family
lizards<-filter(lizards, Library != "Csex_M" & Host_family != "Diplodactylidae")

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#Abundance

#aggregate to get sum of expected_count for each library with all relevant metadata
LizardAbund<-aggregate(expected_count ~Library+Host_infraorder +Host_family + Host_genus+ Biome+ Habitat+ Origins+No._individuals+ Non.host.rRNA_reads , data=lizards, FUN=sum)
#calculate abundance
LizardAbund$Abundance<-(LizardAbund$expected_count/LizardAbund$Non.host.rRNA_reads)*100

## Check effect of differing number of individuals in library
m.No._individuals<-lm(Abundance ~ No._individuals, data=LizardAbund)
print(summary(m.No._individuals))
print(anova(m.No._individuals,test="Chisq"))

#plot for later
p6<-ggplot(LizardAbund, aes(No._individuals, Abundance, fill=Host_family)) + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_colour_manual(values=c( "#5D90FC","#FF9899")) +theme(legend.position = "none")

#plot abundance by family
p2.1<-ggplot(LizardAbund, aes(Host_family, Abundance, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

p4H<-ggplot(LizardAbund, aes(Habitat, Abundance, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Abundance") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

# Let's do some GLMS
# Generalized linear model (GLM) is a generalization of ordinary linear regression that allows for response variables that have error distribution models other than a normal distribution like Gaussian distribution.
# The Gaussian family is is the default for a glm()

## Host family ####
m.1<-glm(Abundance ~ Host_family+Habitat+Biome, data=LizardAbund)
m.1.1<-glm(Abundance ~ Host_family+Habitat+Biome+No._individuals, data=LizardAbund)
m.2<-glm(Abundance ~ Host_family+Habitat, data=LizardAbund)
m.3<-glm(Abundance ~ Host_family+Biome, data=LizardAbund)
m.4<-glm(Abundance ~ Host_family, data=LizardAbund)
m.5<-glm(Abundance ~ Habitat, data=LizardAbund)
m.6<-glm(Abundance ~ Biome, data=LizardAbund)

AIC(m.1, m.1.1, m.2, m.3, m.4, m.5, m.6)
#  If a model is more than 2 AIC units lower than another, then it is considered significantly better than that model.
# Here the AIC’s between our models are very similar. Here model 1, 2, 5 are within 2 AIC. 

print(anova(m.1,test="Chisq"))
print(anova(m.2,test="Chisq"))
print(anova(m.5,test="Chisq"))

# As requested by the reviewer:
# The Poisson family is better for small sample sizes and can be used for count data but herein we don't have counts
# Use quasi-Poisson regression to use Poisson regression on non whole numbers
# Cannot use AIC's to compare models with quasipoisson distributions

m.1<-glm(Abundance ~ Host_family+Habitat+Biome, family = quasipoisson, data=LizardAbund)
m.1.1<-glm(Abundance ~ Host_family+Habitat+Biome+No._individuals, family = quasipoisson, data=LizardAbund)
m.2<-glm(Abundance ~ Host_family+Habitat, family = quasipoisson, data=LizardAbund)
m.3<-glm(Abundance ~ Host_family+Biome, family = quasipoisson, data=LizardAbund)
m.4<-glm(Abundance ~ Host_family, family = quasipoisson, data=LizardAbund)
m.5<-glm(Abundance ~ Habitat, family = quasipoisson, data=LizardAbund)
m.6<-glm(Abundance ~ Biome, family = quasipoisson, data=LizardAbund)

print(anova(m.1,test="Chisq"))
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
## Alpha Diversity

#------------------------------------------------------------------------------------------------------#
# Plots: virus families per sample

### aggregate by viral family ###
Lizard.Family<-aggregate(expected_count ~Library+Host_family + Biome+ Habitat+ Family+ Genus+ No._individuals+ Non.host.rRNA_reads , data=lizards, FUN=sum)
#Drop everything for which abundance is zero?
Lizard.Family[Lizard.Family==0] <- NA
Lizard.Family<-Lizard.Family[complete.cases(Lizard.Family),]

#Calculate "abundance" which is the proportion of total reads
Lizard.Family$Abundance<-(Lizard.Family$expected_count/Lizard.Family$Non.host.rRNA_reads)*100
Lizard.Family$Abundance2<-Lizard.Family$Abundance*10000
#because some metrics require that values be >1, we need to calculate this thing called Abundance2

Lizard.Family2$Abundance<-as.numeric(Lizard.Family2$Abundance)

unique(Lizard.Family$Family)
colors<-cbind("royalblue1",
              "skyblue1",
              "gray85",
              "salmon",
              "red",
              "orange",
              "plum1",
              "#6A3D9A",
              "#33A02C")

names(colors) <- levels(Lizard.Family$Family)
colScale <- scale_fill_manual(name = "Family",values = colors)

#New coloumn which is the concatenations of family and genus
Lizard.Family$Genus2<-paste(Lizard.Family$Family, Lizard.Family$Genus, sep= "-")

#colous
unique(Lizard.Family$Genus)
colors2<-cbind("royalblue1",
               "skyblue1",
               "gray85",
               "salmon",
               "red",
               "#FDBF6F",
               "#FF7F00",
               "plum1",
               "#CAB2D6",
               "#6A3D9A",
               "#B2DF8A", 
               "#33A02C")

names(colors2) <- levels(Lizard.Family$Genus2)
colScale2 <- scale_fill_manual(name = "Genus",values = colors2)


Lizard.Family$Host_family <- factor(Lizard.Family$Host_family, 
                                     levels = c("Scincidae", "Gekkonidae"))

Lizard.Family$Library <- factor(Lizard.Family$Library, levels = c("Cam_M", "Cmun_M", "Ccyg_M", "Gnan_A", "Gnan_M", "Gaus_M", "Hbin_A", "Hbin_M", "Hplan_A"))

p8<-ggplot(Lizard.Family, aes(x=Library, y=Abundance, fill=Family)) + 
  geom_bar(stat="identity")  +
  colScale+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_classic()

p9<-ggplot(Lizard.Family, aes(x=Library, y=Abundance, fill=Genus2)) + 
  geom_bar(stat="identity")  +
  colScale2+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_classic()

p10<-ggplot(Lizard.Family, aes(x=Library, y=Abundance, fill=Family)) + 
  geom_bar(stat="identity", position="fill")  +
  colScale+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_classic()

p11<-ggplot(Lizard.Family, aes(x=Library, y=Abundance, fill=Genus2)) + 
  geom_bar(stat="identity", position="fill")  +
  colScale2+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_classic()

ggarrange(p8, p10, p9, p11, nrow=2, ncol=2, align="hv")
ggsave("AlphaDiv_RelativeAbundance_20221110.pdf", height = 7, width = 13, units = c("in"), dpi = 600)  

#------------------------------------------------------------------------------------------------------#
## Read in and run formulas to calculate a variety of alpha diversity measures
#This script comes from the updated Rhea script set with some minor modification 
# https://github.com/Lagkouvardos/Rhea/blob/master/2.Alpha-Diversity/Alpha-Diversity.R 

# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >1 normalized counts (normalization produces real values for counts)
  count=sum(x[x>1]^0)
  return(count)
}

#' The abundance filtering cutoff 
eff.cutoff <- 0.0025 # this is the default value for Effective Richness (0.25%)

# Calculate the Effective species richness in each individual sample
Eff.Species.richness <- function(x)
{
  # Count only the OTUs that are present more than the set proportion
  total=sum(x)
  count=sum(x[x/total>eff.cutoff]^0)
  return(count)
}

#' The normalized depth cutoff
norm.cutoff <- 1000 # this is the default value for Standard Richness (1000)

# Calculate the Normalized species richness in each individual sample
Norm.Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  # Given a fixed Normalization reads depth
  total=sum(x)
  count=sum(x[norm.cutoff*x/total>0.5]^0)
  return(count)
}

# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
{
  total=sum(x)
  se=-sum(x[x>0]/total*log(x[x>0]/total))
  return(se)
}

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
{
  total=sum(x)
  se=round(exp(-sum(x[x>0]/total*log(x[x>0]/total))),digits =2)
  return(se)
}

# Calculate the Simpson diversity index
Simpson.concentration <- function(x)
{
  total=sum(x)
  si=sum((x[x>0]/total)^2)
  return(si)
}

# Calculate the effective number of species for Simpson
Simpson.effective <- function(x)
{
  total=sum(x)
  si=round(1/sum((x[x>0]/total)^2),digits =2)
  return(si)
}

#------------------------------------------------------------------------------------------------------#
# Run calculations at Virus genus level
# https://github.com/Lagkouvardos/Rhea/blob/master/2.Alpha-Diversity/Alpha-Diversity.R 

#remake otu table with Abundance2
otu_table<- dcast(Lizard.Family, Genus2~Library, value.var = "Abundance2", sum)
otu_table$Family<-NULL #remove all coloumns wiht non-numeric information.. only numbers.
otu_table[is.na(otu_table)] <- 0 ## replace NA with 0

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

# Order and transpose OTU-table
my_otu_table <- otu_table[,order(names(otu_table))] 
my_otu_table <-data.frame(t(my_otu_table))

# Apply diversity functions to table
otus_div_stats<-data.frame(my_otu_table[,0])
otus_div_stats$Richness<-apply(my_otu_table,1,Species.richness)
otus_div_stats$Normalized.Richness<-apply(my_otu_table,1,Norm.Species.richness)
otus_div_stats$Effective.Richness<-apply(my_otu_table,1,Eff.Species.richness)
otus_div_stats$Shannon.Index<-apply(my_otu_table,1,Shannon.entropy)
otus_div_stats$Shannon.Effective<-apply(my_otu_table,1,Shannon.effective)
otus_div_stats$Simpson.Index<-apply(my_otu_table,1,Simpson.concentration)
otus_div_stats$Simpson.Effective<-apply(my_otu_table,1,Simpson.effective)
otus_div_stats$Evenness <- otus_div_stats$Shannon.Index/log(otus_div_stats$Richness,2)

otus_div_stats

#Add in MetaData!
#remake the Lizard.Family metadat table
Lizard.Family<-aggregate(expected_count ~Library+Host_infraorder +Host_family + Host_genus+ Biome+ Habitat+ No._individuals , data=lizards, FUN=sum)
#remove the expectd count colomn
Lizard.Family$expected_count<-NULL
#Lizard.Family<-Lizard.Family[order(Lizard.Family$Library),] #Sort by library name

#the row names here are the library names
df<-tibble::rownames_to_column(otus_div_stats, "Library")
df2<-merge(x = Lizard.Family, y = df, by = "Library", all = TRUE) #merge to get the metadata we want
#convert NA to 0
df2[is.na(df2)] <- 0
write.table(df2, file="AlphaDiversity_RawTable_20221110.txt")

dflong<-melt(df2, id.vars=c("Library","Host_infraorder","Host_family","Host_genus","Biome","Habitat","No._individuals"))

#------------------------------------------------------------------------------------------------------#

## Check No._individuals to appreciate whether this needs to be included in final models. 

m.Pool.Richness<-lm(Richness ~ No._individuals, data=df2)
print(summary(m.Pool.Richness))

m.Pool.Shannon.Index<-lm(Shannon.Index ~ No._individuals, data=df2)
print(summary(m.Pool.Shannon.Index))

m.Pool.Simpson.Index<-lm(Simpson.Index ~ No._individuals, data=df2)
print(summary(m.Pool.Simpson.Index))

m.Pool.Shannon.Effective<-lm(Shannon.Effective ~ No._individuals, data=df2)
print(summary(m.Pool.Shannon.Effective))

m.Pool.Simpson.Effective<-lm(Simpson.Effective ~ No._individuals, data=df2)
print(summary(m.Pool.Simpson.Effective))

# Yes, number of individuals per library is important and should be included in all models.

#RICHNESS
p25<-ggplot(dflong[which(dflong$variable=="Richness"),], aes(No._individuals, value, fill=Host_family)) +  
  geom_point(colour = "black", size = 3.5)+ 
  geom_point(aes(colour = Host_family), size=3)  +
  labs(y = "Richness") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +
  theme(legend.position = "none")+ 
  annotate("text",x = 10, y=4.5,label=c("R^2= 0.4143 , p=0.03644"))

#SHANNON
p26<-ggplot(dflong[which(dflong$variable=="Shannon.Index"),], aes(No._individuals, value, fill=Host_family)) +  
  geom_point(colour = "black", size = 3.5)+ 
  geom_point(aes(colour = Host_family), size=3)  +
  labs(y = "Shannon.Index") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +
  theme(legend.position = "none")+
  annotate("text",x = 10, y=1,label=c("R^2= 0.4596 , p=0.0.02678"))

#SIMPOSON
p27<-ggplot(dflong[which(dflong$variable=="Simpson.Index"),], aes(No._individuals, value, fill=Host_family)) + 
  geom_point(colour = "black", size = 3.5)+ 
  geom_point(aes(colour = Host_family), size=3)  +
  labs(y = "Simpson.Index") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +
  theme(legend.position = "none")+
  annotate("text",x = 10, y=1.1,label=c("R^2= 0.441 , p=0.03047"))

#SHANNON EFFECTIVE
p28<-ggplot(dflong[which(dflong$variable=="Shannon.Effective"),], aes(No._individuals, value, fill=Host_family)) + 
  geom_point(colour = "black", size = 3.5)+ 
  geom_point(aes(colour = Host_family), size=3)  +
  labs(y = "Shannon.Effective ") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +
  theme(legend.position = "none")+
  annotate("text",x = 10, y=3,label=c("R^2= 0.4443 , p=0.02979"))
#Multiple R-squared:  0.3365,	Adjusted R-squared:  0.2628 
#F-statistic: 4.565 on 1 and 9 DF,  p-value: 0.06137

#SIMPSON EFFECTIVE
p29<-ggplot(dflong[which(dflong$variable=="Simpson.Effective"),], aes(No._individuals, value, fill=Host_family)) +  
  geom_point(colour = "black", size = 3.5)+ 
  geom_point(aes(colour = Host_family), size=3)  +
  labs(y = "Simpson.Effective") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +
  theme(legend.position = "none")+
  annotate("text",x = 10, y=2.5,label=c("R^2= 0.4229 , p=0.03443"))

#ABUNDANCE FROM ABOVE
p6.1<-p6+
  annotate("text",x = 10, y=0.25,label=c("R^2=0.1215, p=0.19"))

ggarrange(p6.1, p25, p26, p27, p28, p29, nrow=3, ncol=2, align="hv")
ggsave("Correlation_PoolSize_Abund_Diversity_Effective_202400415.pdf", height = 8, width = 7, units = c("in"), dpi = 600)

#------------------------------------------------------------------------------------------------------#

##PLotting alpha diversity

#Taxonomy
#RICHNESS
p11<-ggplot(dflong[which(dflong$variable=="Richness"),], aes(Host_family, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Richness") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

#SHANNON
p14<-ggplot(dflong[which(dflong$variable=="Shannon.Index"),], aes(Host_family, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Shannon Index") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

#SIMPSON
p17<-ggplot(dflong[which(dflong$variable=="Simpson.Index"),], aes(Host_family, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Simpson Index") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

#SHANNON EFFECTIVE
p20<-ggplot(dflong[which(dflong$variable=="Shannon.Effective"),], aes(Host_family, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Shannon Effective") +
  scale_colour_manual(values=c( "#5D90FC","#FF9899")) +theme(legend.position = "none")

#SIMPSON
p23<-ggplot(dflong[which(dflong$variable=="Simpson.Effective"),], aes(Host_family, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Simpson Effective") +
  scale_colour_manual(values=c( "#5D90FC","#FF9899")) +theme(legend.position = "none")

#export 14x10
ggarrange(p2.1, p11, p14, p17, p20, p23,  
          nrow=3, ncol=2, align="hv")
ggsave("AlphaDiversity_HostTaxonomy_Effectives_20240415.pdf", height = 12, width = 9, units = c("in"), dpi = 600)  


#RICHNESS
#By Habitat
p30<-ggplot(dflong[which(dflong$variable=="Richness"),], aes(Habitat, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Richness") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

#SHANNON
#By Habitat
p32<-ggplot(dflong[which(dflong$variable=="Shannon.Index"),], aes(Habitat, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Shannon.Index") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

#SIMPSON
#By Habitat
p34<-ggplot(dflong[which(dflong$variable=="Simpson.Index"),], aes(Habitat, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Simpson.Index") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

#SHANNONeffective
#By Habitat
p32.1<-ggplot(dflong[which(dflong$variable=="Shannon.Effective"),], aes(Habitat, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Shannon.Effective") +
  scale_colour_manual(values=c( "#5D90FC","#FF9899")) +theme(legend.position = "none")

#SIMPSONeffective
#By Habitat
p34.1<-ggplot(dflong[which(dflong$variable=="Simpson.Effective"),], aes(Habitat, value, fill=Host_family)) + geom_boxplot(alpha=0.7, fill="white") + geom_point(colour = "black", size = 3.5)+ geom_point(aes(colour = Host_family), size=3)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(y = "Simpson.Effective") +
  scale_colour_manual(values=c("#5D90FC","#FF9899")) +theme(legend.position = "none")

ggarrange(p4H, p30, p32,p34,p32.1,p34.1,
          nrow=3, ncol=2, align="hv")
ggsave("AlphaDiversity_Habitat_20240415.pdf", height = 9, width = 12.17, units = c("in"), dpi = 600) 

ggarrange(p12,p30, nrow=2, ncol=1, align="hv")
ggsave("tmp.pdf", height = 9, width =4.5 , units = c("in"), dpi = 600)

#------------------------------------------------------------------------------------------------------#
## Stats model selection

df2$Host_family<-as.factor(df2$Host_family)
df2$Habitat<-as.factor(df2$Habitat)
df2$Biome<-as.factor(df2$Biome)

#RICHNESS
m.1.Richness<-glm(Richness ~ Host_family+Biome+Habitat+No._individuals,data=df2)
m.2.Richness<-glm(Richness ~ Host_family+Biome+No._individuals, data=df2)
m.3.Richness<-glm(Richness ~ Host_family+Habitat+No._individuals, data=df2)
m.4.Richness<-glm(Richness ~ Host_family+No._individuals, data=df2)
m.5.Richness<-glm(Richness ~ Habitat+No._individuals, data=df2)
m.6.Richness<-glm(Richness ~ Biome+No._individuals, data=df2)
m.7.Richness<-glm(Richness ~ Host_family, data=df2)
m.8.Richness<-glm(Richness ~ Habitat, data=df2)
m.9.Richness<-glm(Richness ~ Biome, data=df2)
m.10.Richness<-glm(Richness ~ No._individuals, data=df2)

AIC(m.1.Richness, m.2.Richness, m.3.Richness, m.4.Richness, m.5.Richness, m.6.Richness, m.7.Richness, m.8.Richness, m.9.Richness, m.10.Richness) 

print(anova(m.3.Richness,test="Chisq"))
summary(glht(m.3.Richness,linfct = mcp(Host_family="Tukey")))
summary(glht(m.3.Richness,linfct = mcp(Habitat="Tukey")))

# Rerun with poisson distribution to account for small sample size at request of reviewer
m.1.Richness<-glm(Richness ~ Host_family+Biome+Habitat+No._individuals,family=quasipoisson,data=df2)
m.2.Richness<-glm(Richness ~ Host_family+Biome+No._individuals, family=quasipoisson,data=df2)
m.3.Richness<-glm(Richness ~ Host_family+Habitat+No._individuals,family=quasipoisson, data=df2)
m.4.Richness<-glm(Richness ~ Host_family+No._individuals, family=quasipoisson,data=df2)
m.5.Richness<-glm(Richness ~ Habitat+No._individuals, family=quasipoisson,data=df2)
m.6.Richness<-glm(Richness ~ Biome+No._individuals, family=quasipoisson,data=df2)
m.7.Richness<-glm(Richness ~ Host_family, family=quasipoisson,data=df2)
m.8.Richness<-glm(Richness ~ Habitat, family=quasipoisson,data=df2)
m.9.Richness<-glm(Richness ~ Biome, family=quasipoisson,data=df2)
m.10.Richness<-glm(Richness ~ No._individuals, family=quasipoisson,data=df2)

print(anova(m.3.Richness,test="Chisq"))

m.3.Richness<-glm(Richness ~ Host_family+Habitat+No._individuals,family=quasipoisson, data=df2)
print(anova(m.3.Richness,test="Chisq"))
summary(glht(m.3.Richness,linfct = mcp(Host_family="Tukey")))
summary(glht(m.3.Richness,linfct = mcp(Habitat="Tukey")))

#SHANNON
m.1.Shannon.Index<-glm(Shannon.Index ~ Host_family+Biome+Habitat+No._individuals,data=df2)
m.2.Shannon.Index<-glm(Shannon.Index ~ Host_family+Biome+No._individuals, data=df2)
m.3.Shannon.Index<-glm(Shannon.Index ~ Host_family+Habitat+No._individuals, data=df2)
m.4.Shannon.Index<-glm(Shannon.Index ~ Host_family+No._individuals, data=df2)
m.5.Shannon.Index<-glm(Shannon.Index ~ Habitat+No._individuals, data=df2)
m.6.Shannon.Index<-glm(Shannon.Index ~ Biome+No._individuals, data=df2)
m.7.Shannon.Index<-glm(Shannon.Index ~ Host_family, data=df2)
m.8.Shannon.Index<-glm(Shannon.Index ~ Habitat, data=df2)
m.9.Shannon.Index<-glm(Shannon.Index ~ Biome, data=df2)
m.10.Shannon.Index<-glm(Shannon.Index ~ No._individuals, data=df2)

AIC(
  m.1.Shannon.Index, 
  m.2.Shannon.Index, 
  m.3.Shannon.Index, 
  m.4.Shannon.Index, 
  m.5.Shannon.Index, 
  m.6.Shannon.Index, 
  m.7.Shannon.Index, 
  m.8.Shannon.Index, 
  m.9.Shannon.Index, 
  m.10.Shannon.Index)

print(anova(m.3.Shannon.Index,test="Chisq"))
summary(glht(m.3.Shannon.Index,linfct = mcp(Host_family="Tukey")))
summary(glht(m.3.Shannon.Index,linfct = mcp(Habitat="Tukey")))

# Rerun with poisson distribution to account for small sample size at request of reviewer
m.1.Shannon.Index<-glm(Shannon.Index ~ Host_family+Biome+Habitat+No._individuals,family=quasipoisson, data=df2)
m.2.Shannon.Index<-glm(Shannon.Index ~ Host_family+Biome+No._individuals, family=quasipoisson,data=df2)
m.3.Shannon.Index<-glm(Shannon.Index ~ Host_family+Habitat+No._individuals, family=quasipoisson,data=df2)
m.4.Shannon.Index<-glm(Shannon.Index ~ Host_family+No._individuals, family=quasipoisson,data=df2)
m.5.Shannon.Index<-glm(Shannon.Index ~ Habitat+No._individuals, family=quasipoisson,data=df2)
m.6.Shannon.Index<-glm(Shannon.Index ~ Biome+No._individuals, family=quasipoisson,data=df2)
m.7.Shannon.Index<-glm(Shannon.Index ~ Host_family, family=quasipoisson,data=df2)
m.8.Shannon.Index<-glm(Shannon.Index ~ Habitat, family=quasipoisson,data=df2)
m.9.Shannon.Index<-glm(Shannon.Index ~ Biome, family=quasipoisson,data=df2)
m.10.Shannon.Index<-glm(Shannon.Index ~ No._individuals, family=quasipoisson,data=df2)

print(anova(m.3.Shannon.Index,test="Chisq"))
summary(glht(m.3.Shannon.Index,linfct = mcp(Host_family="Tukey")))
summary(glht(m.3.Shannon.Index,linfct = mcp(Habitat="Tukey")))

#SIMPSON
m.1.Simpson.Index<-glm(Simpson.Index ~ Host_family+Biome+Habitat+No._individuals,data=df2)
m.2.Simpson.Index<-glm(Simpson.Index ~ Host_family+Biome+No._individuals, data=df2)
m.3.Simpson.Index<-glm(Simpson.Index ~ Host_family+Habitat+No._individuals, data=df2)
m.4.Simpson.Index<-glm(Simpson.Index ~ Host_family+No._individuals, data=df2)
m.5.Simpson.Index<-glm(Simpson.Index ~ Habitat+No._individuals, data=df2)
m.6.Simpson.Index<-glm(Simpson.Index ~ Biome+No._individuals, data=df2)
m.7.Simpson.Index<-glm(Simpson.Index ~ Host_family, data=df2)
m.8.Simpson.Index<-glm(Simpson.Index ~ Habitat, data=df2)
m.9.Simpson.Index<-glm(Simpson.Index ~ Biome, data=df2)
m.10.Simpson.Index<-glm(Simpson.Index ~ No._individuals, data=df2)

AIC(
  m.1.Simpson.Index, 
  m.2.Simpson.Index, 
  m.3.Simpson.Index, 
  m.4.Simpson.Index, 
  m.5.Simpson.Index, 
  m.6.Simpson.Index, 
  m.7.Simpson.Index, 
  m.8.Simpson.Index, 
  m.9.Simpson.Index, 
  m.10.Simpson.Index)

print(anova(m.3.Simpson.Index,test="Chisq"))
summary(glht(m.3.Simpson.Index,linfct = mcp(Habitat="Tukey")))

#rerun with quasipoisson

m.1.Simpson.Index<-glm(Simpson.Index ~ Host_family+Biome+Habitat+No._individuals,family=quasipoisson,data=df2)
m.2.Simpson.Index<-glm(Simpson.Index ~ Host_family+Biome+No._individuals, family=quasipoisson,data=df2)
m.3.Simpson.Index<-glm(Simpson.Index ~ Host_family+Habitat+No._individuals, family=quasipoisson,data=df2)
m.4.Simpson.Index<-glm(Simpson.Index ~ Host_family+No._individuals, family=quasipoisson,data=df2)
m.5.Simpson.Index<-glm(Simpson.Index ~ Habitat+No._individuals, family=quasipoisson,data=df2)
m.6.Simpson.Index<-glm(Simpson.Index ~ Biome+No._individuals,family=quasipoisson, data=df2)
m.7.Simpson.Index<-glm(Simpson.Index ~ Host_family, family=quasipoisson,data=df2)
m.8.Simpson.Index<-glm(Simpson.Index ~ Habitat, family=quasipoisson,data=df2)
m.9.Simpson.Index<-glm(Simpson.Index ~ Biome, family=quasipoisson,data=df2)
m.10.Simpson.Index<-glm(Simpson.Index ~ No._individuals, family=quasipoisson,data=df2)

m.3.Simpson.Index<-glm(Simpson.Index ~ Host_family+Habitat+No._individuals, family=quasipoisson, data=df2)
print(anova(m.3.Simpson.Index,test="Chisq"))

#SHANNON EFFECTIVE
m.1.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+Biome+Habitat+No._individuals,data=df2)
m.2.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+Biome+No._individuals, data=df2)
m.3.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+Habitat+No._individuals, data=df2)
m.4.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+No._individuals, data=df2)
m.5.Shannon.Effective<-glm(Shannon.Effective ~ Habitat+No._individuals, data=df2)
m.6.Shannon.Effective<-glm(Shannon.Effective ~ Biome+No._individuals, data=df2)
m.7.Shannon.Effective<-glm(Shannon.Effective ~ Host_family, data=df2)
m.8.Shannon.Effective<-glm(Shannon.Effective ~ Habitat, data=df2)
m.9.Shannon.Effective<-glm(Shannon.Effective ~ Biome, data=df2)
m.10.Shannon.Effective<-glm(Shannon.Effective ~ No._individuals, data=df2)

AIC(
  m.1.Shannon.Effective, 
  m.2.Shannon.Effective, 
  m.3.Shannon.Effective, 
  m.4.Shannon.Effective, 
  m.5.Shannon.Effective, 
  m.6.Shannon.Effective, 
  m.7.Shannon.Effective, 
  m.8.Shannon.Effective, 
  m.9.Shannon.Effective, 
  m.10.Shannon.Effective)

print(anova(m.5.Shannon.Effective,test="Chisq"))
summary(glht(m.5.Shannon.Effective,linfct = mcp(Habitat="Tukey")))

# Rerun with poisson distribution to account for small sample size at request of reviewer
m.1.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+Biome+Habitat+No._individuals,family=quasipoisson, data=df2)
m.2.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+Biome+No._individuals, family=quasipoisson,data=df2)
m.3.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+Habitat+No._individuals, family=quasipoisson,data=df2)
m.4.Shannon.Effective<-glm(Shannon.Effective ~ Host_family+No._individuals, family=quasipoisson,data=df2)
m.5.Shannon.Effective<-glm(Shannon.Effective ~ Habitat+No._individuals, family=quasipoisson,data=df2)
m.6.Shannon.Effective<-glm(Shannon.Effective ~ Biome+No._individuals, family=quasipoisson,data=df2)
m.7.Shannon.Effective<-glm(Shannon.Effective ~ Host_family, family=quasipoisson,data=df2)
m.8.Shannon.Effective<-glm(Shannon.Effective ~ Habitat, family=quasipoisson,data=df2)
m.9.Shannon.Effective<-glm(Shannon.Effective ~ Biome, family=quasipoisson,data=df2)
m.10.Shannon.Effective<-glm(Shannon.Effective ~ No._individuals, family=quasipoisson,data=df2)

print(anova(m.5.Shannon.Effective,test="Chisq"))
summary(glht(m.5.Shannon.Effective,linfct = mcp(Habitat="Tukey")))

#SIMPSON
m.1.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+Biome+Habitat+No._individuals,data=df2)
m.2.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+Biome+No._individuals, data=df2)
m.3.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+Habitat+No._individuals, data=df2)
m.4.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+No._individuals, data=df2)
m.5.Simpson.Effective<-glm(Simpson.Effective ~ Habitat+No._individuals, data=df2)
m.6.Simpson.Effective<-glm(Simpson.Effective ~ Biome+No._individuals, data=df2)
m.7.Simpson.Effective<-glm(Simpson.Effective ~ Host_family, data=df2)
m.8.Simpson.Effective<-glm(Simpson.Effective ~ Habitat, data=df2)
m.9.Simpson.Effective<-glm(Simpson.Effective ~ Biome, data=df2)
m.10.Simpson.Effective<-glm(Simpson.Effective ~ No._individuals, data=df2)

AIC(
  m.1.Simpson.Effective, 
  m.2.Simpson.Effective, 
  m.3.Simpson.Effective, 
  m.4.Simpson.Effective, 
  m.5.Simpson.Effective, 
  m.6.Simpson.Effective, 
  m.7.Simpson.Effective, 
  m.8.Simpson.Effective, 
  m.9.Simpson.Effective, 
  m.10.Simpson.Effective)

print(anova(m.5.Simpson.Effective,test="Chisq"))
summary(glht(m.5.Simpson.Effective,linfct = mcp(Habitat="Tukey")))

#rerun with quasipoisson

m.1.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+Biome+Habitat+No._individuals,family=quasipoisson,data=df2)
m.2.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+Biome+No._individuals, family=quasipoisson,data=df2)
m.3.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+Habitat+No._individuals, family=quasipoisson,data=df2)
m.4.Simpson.Effective<-glm(Simpson.Effective ~ Host_family+No._individuals, family=quasipoisson,data=df2)
m.5.Simpson.Effective<-glm(Simpson.Effective ~ Habitat+No._individuals, family=quasipoisson,data=df2)
m.6.Simpson.Effective<-glm(Simpson.Effective ~ Biome+No._individuals,family=quasipoisson, data=df2)
m.7.Simpson.Effective<-glm(Simpson.Effective ~ Host_family, family=quasipoisson,data=df2)
m.8.Simpson.Effective<-glm(Simpson.Effective ~ Habitat, family=quasipoisson,data=df2)
m.9.Simpson.Effective<-glm(Simpson.Effective ~ Biome, family=quasipoisson,data=df2)
m.10.Simpson.Effective<-glm(Simpson.Effective ~ No._individuals, family=quasipoisson,data=df2)

print(anova(m.5.Simpson.Effective,test="Chisq"))
summary(glht(m.5.Simpson.Effective,linfct = mcp(Habitat="Tukey")))


#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
# Beta Diversity

#we are working at virus family level here

#Prepare data and make phyloseq object

#Make basic map and lizard files
Lizard.tmp<-aggregate(expected_count ~Library+Host_infraorder +Host_family + Host_genus+ Biome+ Habitat+Family+No._individuals+ Non.host.rRNA_reads , data=lizards, FUN=sum)

#remove for all where expected count is zero
#Drop everything for which abundance is zero?
Lizard.tmp[Lizard.tmp==0] <- NA
Lizard.tmp<-Lizard.tmp[complete.cases(Lizard.tmp),]

#Calculate "abundance" which is the proportion of total reads
Lizard.tmp$Abundance<-(Lizard.tmp$expected_count/Lizard.tmp$Non.host.rRNA_reads)*100
Lizard.tmp$Abundance2<-Lizard.tmp$Abundance*100000

map2<-aggregate(expected_count ~Host_infraorder +Host_family + Host_genus+ Biome+ Habitat+Range_size+Origins+ No._individuals+Library, data=lizards, FUN=length)
map2$expected_count<-NULL

### Ordination with Phyloseq, which is EXTREMELY ANNOYING! #####
#Need to create a phyloseq object, which is aweful and is a long and annnoying process
#justphyloseqthings
#https://github.com/joey711/phyloseq?tab=readme-ov-file

#need to create a metadata table
#pay careful attention to the order in which teh libraries appear so that it all correponsonds with "sa1-sa11"
map3 <- sample_data(map3)##make meta data into phyloseq format

#create an OTU table & Tax table
tmp2<- dcast(Lizard.tmp, Family~Library, value.var = "Abundance2") #tip over the table
tmp2[is.na(tmp2)] <- 0
tax2 <- tmp2$Family # pull out virus family names for the tax table below

#OTU table
tmp2$Family<-NULL #remove the virus family names
as.matrix(tmp2) #convert to a matrix
#rename col names as defined in tax2
colnames(tmp2)<-c("sa1", "sa2", "sa3", "sa4", "sa5", "sa6", "sa7", "sa8", "sa9")

OTU = otu_table(tmp2, taxa_are_rows = TRUE)

#tax table
TAX<-as.matrix(tax2)
TAX2 = tax_table(TAX)

#make phyloseq object
physeq = phyloseq(OTU, TAX2, map3)

## ordinate
lizard_bray <- ordinate(
  physeq = physeq,
  method = "NMDS",
  distance = "bray"
)

plot_ordination(physeq, lizard_bray, color="Library") + geom_point(alpha = 1, size = 5, stroke = 2) 
#Cmun_M is clearly an outlier here throwing off the entire thing. The library contains only bornavirus, which is absent in all other libraries

#------------------------------------------------------------------------------------------------------#
## nMDS for dataset with Cmun_M removed as its clearly an outlier here throwing off the entire thing

map3<-filter(map2, Library!="Cmun_M")
map3 <- sample_data(map3)##make meta data into phyloseq format

#create an OTU table & Tax table
Lizard.tmp2<-filter(Lizard.tmp, Library!="Cmun_M")
tmp2<- dcast(Lizard.tmp2, Family~Library, value.var = "Abundance2") #tip over the table
tmp2[is.na(tmp2)] <- 0
tax2 <- tmp2$Family # pull out virus family names for the tax table below

#OTU table
tmp2$Family<-NULL #remove the virus family names
as.matrix(tmp2) #convert to a matrix
#rename col names as defined in tax2
colnames(tmp2)<-c("sa1", "sa2", "sa3", "sa4", "sa5", "sa6", "sa7", "sa8")

OTU = otu_table(tmp2, taxa_are_rows = TRUE)

#tax table
TAX<-as.matrix(tax2)
TAX2 = tax_table(TAX)

#make phyloseq object
physeq = phyloseq(OTU, TAX2, map3)

## ordinate
lizard_bray <- ordinate(
  physeq = physeq,
  method = "NMDS",
  distance = "bray"
)

plot_ordination(physeq, lizard_bray, color="Library") + geom_point(alpha = 1, size = 5, stroke = 2) 
#CCyg which is not visitble is under Gnan_A
#------------------------------------------------------------------------------------------------------#

#Host_family
p.nMDS.Family<-plot_ordination(physeq, lizard_bray, color="Host_family") + geom_point(alpha = 1, size = 5, stroke = 2)+ scale_color_manual(values=c("#5D90FC","#FF9899")) + theme_bw() 
#CCyg which is not visible is under Gnan_A. Manually fix

#Habitat
p.nMDS.Habitat<-plot_ordination(physeq, lizard_bray, color="Habitat") + geom_point(alpha = 1, size = 5, stroke = 2)+ scale_color_manual(values=c("#AFE1AF", "#50C878","#097969")) + theme_bw() 
#CCyg which is not visible is under Gnan_A. Manualy fix.

ggarrange(p.nMDS.Family, p.nMDS.Habitat, nrow=1, ncol=2, align="hv")

##stats
lizard_bray <- phyloseq::distance(physeq, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq))

#Of note is that you cannot compare the outputs of PERMANOVAs using AICs
#So, let's build a full model and drop terms that are not signficaint. 

vegan::adonis2(lizard_bray ~ Host_family+Habitat+Biome, data = sampledf)
vegan::adonis2(lizard_bray ~ Host_family+Habitat, data = sampledf)

#------------------------------------------------------------------------------------------------------#


