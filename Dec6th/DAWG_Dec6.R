#### Setup ####

### Clear workspace ###
rm(list=ls())

### Install Packages ###
#Corncob Tutorial and more info: https://rdrr.io/github/bryandmartin/corncob/f/vignettes/corncob-intro.Rmd
#DESeq2 Tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html

# if you don't have devtools installed, install it

# devtools for mac/linux
devtools::install_github("hadley/devtools")
# devtools for windows

library(devtools)
build_github_devtools()

# restart R #

devtools::install_github("bryandmartin/corncob",force=TRUE)

#Selbal paper: https://doi.org/10.1101/219386
#Selbal vignettes: https://htmlpreview.github.io/?https://github.com/UVic-omics/selbal/blob/master/vignettes/vignette.html
devtools::install_url(url="https://github.com/UVic-omics/selbal/archive/master.zip",
                      INSTALL_opt= "--no-multiarch")
library("selbal")
library("grid")
library("phyloseq")

BiocManager::install("DESeq2")
#OR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("corncob","phyloseq", "ggplot2","DESeq2", "selbal", "grid")
ipak(packages)

### Prep Phyloseq Objects ###
#genus level
#ps.gen<-taxa_level(ps.adj,"Genus") #microbiomeSeq package

#Separate by sample type
ps.gen.stool<-subset_samples(ps.gen, Sample =="stool")
ps.gen.trach<-subset_samples(ps.gen, Sample =="tracheal")

#subset for appropriate processing
ps.gen.pair<-subset_samples(ps.gen, Paired == 1)
ps.gen.stool.death<-subset_samples(ps.gen.stool, Death !="NA")


#### Comparing Corncob to DESeq2####
### Corncob ###
set.seed(1)
da_Stool <- differentialTest(formula = ~ Sepsis, #the DA formula - if you have a controlling factor it must be included. 
                                phi.formula = ~ 1, #the DV formula - if you have a controlling factor it must be included
                                formula_null = ~ 1, # DA controlling factors
                                phi.formula_null = ~ 1, #DV controlling factors
                                test = "Wald", boot = FALSE, #wald test is "standard"
                                data = ps.gen.stool, #PS object
                                fdr_cutoff = 0.05) #pval false discovery value
#explore significant taxa
da_Stool$significant_taxa
da_Stool$p_fdr
plot(da_Stool)

### DEseq2 ###
diagdds = phyloseq_to_deseq2(ps.gen.stool, ~ Sepsis) #model statement - here is where you can add additional controlling factors

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

#Filter results for p-values<0.05
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), "matrix")
sigtab$genus<-rownames(sigtab)
head(sigtab)

#order by genus
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
#plot
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
ggplot(sigtab, aes(x=genus, y=log2FoldChange,color=baseMean)) + 
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#compare results
da_Stool$significant_taxa
sigtab$genus

#### Corncob: test Death, controlling for Sepsis ####
#by controlling sepsis, we my be able to tease out genera that lead to poor health outside of sepsis
set.seed(1)
da_Stool_death <- differentialTest(formula = ~ Death+Sepsis, #the DA formula - if you have a controlling factor it must be included. 
                             phi.formula = ~ Death+Sepsis, #the DV formula - if you have a controlling factor it must be included
                             formula_null = ~ Sepsis, # DA controlling factors
                             phi.formula_null = ~ Sepsis, #DV controlling factors
                             test = "Wald", boot = FALSE, #wald test is "standard"
                             data = ps.gen.stool.death, #PS object
                             fdr_cutoff = 0.05) #pval false discovery value
#explore significant taxa
da_Stool_death$significant_taxa
da_Stool_death$significant_models
da_Stool_death$p_fdr
plot(da_Stool_death)

####Relative abunance of Sepsis Taxa####
#Descriptive plots of taxa that are more or less abundant based on Sepsis
#Taxa chosen based on DA analysis and/or ones known to be associated with sepsis
ps.perc<-transform_sample_counts(ps.adj, function(x) x / sum(x)) #This calculates the percentage of each read in each sample
meta<-sample_data(ps.adj)

Pseudo <- subset_taxa(ps.perc, Genus == "Pseudomonas") #many species, P. aeruginosa is second most common infection in hospitals
Pseudo.sum<-rowSums(otu_table(Pseudo))
Sub <- subset_taxa(ps.perc, Genus == "Subdoligranulum") #found in stool/common in gut
Sub.sum<-rowSums(otu_table(Sub))
Kleb <- subset_taxa(ps.perc, Genus == "Klebsiella") # opportunistic pathogens: can lead to a wide range of disease states, notably pneumonia, urinary tract infections, sepsis, meningitis, diarrhea, and soft tissue infections.
Kleb.sum<-rowSums(otu_table(Kleb))
Esch <- subset_taxa(ps.perc, Genus == "Escherichia/Shigella") #E. coli, potentially pathogenic (not all Escherichia are)
Esch.sum<-rowSums(otu_table(Esch))
Bac<- subset_taxa(ps.perc, Genus == "Bacteroides") #important gut bacteria (obligate?) for digestion
Bac.sum<-rowSums(otu_table(Bac))
Aci <- subset_taxa(ps.perc, Genus == "Acinetobacter") #key source of infection in hospitalized patients - pneumonia, UTI, bacteremia, secondary menegitis
Aci.sum<-rowSums(otu_table(Aci))
Ent <- subset_taxa(ps.perc, Genus == "Enterococcus") #frequently antibiotic resistant, associated with several diseseases (UTI,Meningitis)
Ent.sum<-rowSums(otu_table(Ent))
Blau <- subset_taxa(ps.perc, Genus == "Blautia") #common and potentially beneficial gut bacteria
Blau.sum<-rowSums(otu_table(Blau))

gen.mat<-cbind(Pseudo.sum,Sub.sum,Kleb.sum,Esch.sum,Bac.sum,Aci.sum,Ent.sum,Blau.sum)*100

###Sepsis Relative Abundance
#Aggregate by sepsis
gen.ag<-aggregate(gen.mat~meta$Sepsis, FUN=mean) #change parameter in this line for other plots

#Add rownames
rownames(gen.ag)<-gen.ag[,1]
gen.ag.2<-as.matrix(gen.ag[,-1])
gen.ag.names<-rownames(gen.ag)


#Plot creation
quartz() #or windows()
par(mar=par()$mar+c(0,0,0,10))
barplot(apply(t(gen.ag.2),2,rev),names=gen.ag.names,
        col=c("gainsboro","midnightblue","mediumturquoise","mediumslateblue","darkmagenta",
              "firebrick2","gold","forestgreen"),
        horiz=FALSE, ylab="Relative abundance (%)",xlab="Sepsis",
        main="Differential Bacteria Genera")
legend(3,20,xpd=TRUE,
       legend = c("Pseudomonas","Subdoligranulum","Klebsiella","Escherichia/Shigella",
                  "Bacteroides","Acinetobacter","Enterococcus","Blautia"), 
       fill = c("forestgreen","gold","firebrick2","darkmagenta",
                "mediumslateblue","mediumturquoise","midnightblue","gainsboro"))

####Relative abunance of Death Taxa####
#Descriptive plots of taxa that are more or less abundant for Death when controlled by "sepsis"
Pedo <- subset_taxa(ps.perc, Genus == "Pedobacter") #common DNA extraction kit contaminant
Pedo.sum<-rowSums(otu_table(Pedo))
Sub <- subset_taxa(ps.perc, Genus == "Subdoligranulum") 
Sub.sum<-rowSums(otu_table(Sub))
Trep <- subset_taxa(ps.perc, Genus == "Treponema_2") #Associated with syphilis, bejel, and yaws
Trep.sum<-rowSums(otu_table(Trep))
Esch <- subset_taxa(ps.perc, Genus == "Escherichia/Shigella") 
Esch.sum<-rowSums(otu_table(Esch))
Sel<- subset_taxa(ps.perc, Genus == "Selenomonas") #common in gastrointestinal tracts of animals
Sel.sum<-rowSums(otu_table(Sel))
Fre <- subset_taxa(ps.perc, Genus == "Fretibacterium") #one known species, associated with periodontitis
Fre.sum<-rowSums(otu_table(Fre))
Allo <- subset_taxa(ps.perc, Genus == "Alloprevotella") #common in gut, some in it's shared family (Prevotella) may be opportunistic pathogen
Allo.sum<-rowSums(otu_table(Allo))
Scar <- subset_taxa(ps.perc, Genus == "Scardovia") #potential dental pathogen
Scar.sum<-rowSums(otu_table(Scar))

gen.mat.de<-cbind(Pedo.sum,Sub.sum,Trep.sum,Esch.sum,Sel.sum,Fre.sum,Allo.sum,Scar.sum)*100

###Sepsis Relative Abundance
#Aggregate by sepsis
gen.ag.de<-aggregate(gen.mat~meta$Death_cat, FUN=mean) #change parameter in this line for other plots

#Add rownames
rownames(gen.ag.de)<-gen.ag.de[,1]
gen.ag.2.de<-as.matrix(gen.ag.de[,-1])
gen.ag.names.de<-rownames(gen.ag.de)


#Plot creation
quartz() #or windows()
par(mar=par()$mar+c(0,0,0,10))
barplot(apply(t(gen.ag.2.de),2,rev),names=gen.ag.names.de,
        col=c("gainsboro","midnightblue","mediumturquoise","mediumslateblue","darkmagenta",
              "firebrick2","gold","forestgreen"),
        horiz=FALSE, ylab="Relative abundance (%)",xlab="Alive >50 Days",
        main="Differential Bacteria Genera")
legend(2.5,15,xpd=TRUE,
       legend = c("Pedobacter","Subdoligranulum","Treponema_2","Escherichia/Shigella",
                  "Selenomonas","Fretibacterium","Alloprevotella","Scardovia"), 
       fill = c("forestgreen","gold","firebrick2","darkmagenta",
                "mediumslateblue","mediumturquoise","midnightblue","gainsboro"))

####IF TIME - Explore "paired" samples####
#what if you don't control for variability?
set.seed(1)
da_analysis <- differentialTest(formula = ~ Date+Sample, #the DA formula - if you have a controlling factor it must be included. 
                                phi.formula = ~ Date+Sample, #the DV formula - if you have a controlling factor it must be included
                                formula_null = ~Sample, # DA controlling factors
                                phi.formula_null = ~ Sample, #DV controlling factors
                                test = "Wald", boot = FALSE, #wald test is "standard"; it is unlikely you will bootstrap
                                data = ps.gen.pair, #PS object
                                fdr_cutoff = 0.10) #pval false discovery value
#explore significant taxa
da_analysis$significant_taxa
#actual p_fdr (adjusted p values) for all tests
da_analysis$p_fdr
plot(da_analysis)

#### Selbal analysis ####

# Rarefy
set.seed(1234)
ps.rarefied <- rrarefy(otu_table(ps.pruned), min(rowSums(otu_table(ps.pruned))))

map <- sample_data(ps.pruned)
map$month <- as.factor(map$month) #binary response variable must be a factor, not numeric
map$monthChar <- map$month
levels(map$monthChar) <- c("Survival","Mortality")
tax <- tax_table(ps.pruned)
ps <- merge_phyloseq(sample_data(map), tax_table(tax), otu_table(ps.rarefied, taxa_are_rows=FALSE))

genus_level <- taxa_level(ps, "Genus") #(554 genera)

# Keep genera present in at least 10% of the samples
# (This is necesary because Selbal uses geometric Bayesian multiplicative replacemnt 
# to eliminate zeros in the OTU table, which doesn't work if there are too many zeros to start with)
few.gen <- otu_table(genus_level)[,which(colSums(otu_table(genus_level)!=0)>(205*.1))]
dim(few.gen) # 93 genera

plot(jitter(map$APACHE, factor=0.5),jitter(as.numeric(map$month)-1,factor=0.5),
     xlab="APACHE", ylab="0=survival, 1=mortality in 1st month") # some correlation
plot(jitter(map$ARDS, factor=0.5),jitter(as.numeric(map$month)-1,factor=0.5),
     xlab="ARDS", ylab="0=survival, 1=mortality in 1st month") 
plot(jitter(map$Shock, factor=0.2),jitter(as.numeric(map$month)-1,factor=0.5),
     xlab="Shock", ylab="0=survival, 1=mortality in 1st month") 
plot(jitter(map$Age, factor=0.2),jitter(as.numeric(map$month)-1,factor=0.5),
     xlab="Age", ylab="0=survival, 1=mortality in 1st month") # some correlation
#boxplots
plot(map$monthChar, map$APACHE)
plot(map$monthChar, map$Age)

### stool only ###
few.stool <- few.gen[which(map$Sample=="stool"),] # 102 samples
## month (categorical) ##
month.stool <- map$month[which(map$Sample=="stool")]
# controling for confounding factors (i.e. APACHE score and Age)
z <- cbind(map$Age, map$APACHE)
colnames(z)<-c("Age","APACHE") # this is important (no column names gives an error in the selbal.cv function)
# what is the optimal number of taxa to be included in the balance? cross-validation
start_time <- Sys.time()
CVbalStool <- selbal.cv(few.stool, month.stool, 
                        covar = z[which(map$Sample=="stool"),], seed = 1234) # you can run the function without covariates as well
end_time <- Sys.time()
end_time - start_time # Time difference 17.9716 mins

CVbalStool$accuracy.nvar # check the smallest # of variables above a minimum Accuracy: 2
CVbalStool$var.barplot # check for which taxa are included in numerator/denominator
# Denoninator: Parasutterella
# Numerator: Haemophilus
plot(CVbalStool$global.plot)
plot(CVbalStool$ROC.plot)
# check how robust the balance variable selection is
plot.tab(CVbalStool$cv.tab)
# the global balance (Haemophilus/Parasutterella) is most frequently slected (30%)
summary(CVbalStool$cv.accuracy) #mean=cross validation accuracy
CVbalStool$glm

## obtain "global balance" ##
start_time <- Sys.time()
gbStool <- selbal(x=few.stool, y=month.stool, covar = z[which(map$Sample=="stool"),], maxV = 2)
end_time <- Sys.time()
end_time - start_time # Time difference 1.891175 mins

plot(gbStool$global.plot) # (Haemophilus/Parasutterella)
plot(gbStool$ROC.plot) # 0.758

## death (continuous) ##
death.stool <- map$Death[which(map$Sample=="stool")]
start_time <- Sys.time()
CVbalStool2 <- selbal.cv(few.stool[which(!is.na(death.stool)),], 
                         death.stool[which(!is.na(death.stool))], 
                         covar = z[which(map$Sample=="stool"),][which(!is.na(death.stool)),], 
                         seed = 1234) 
end_time <- Sys.time()
end_time - start_time # Time difference 39.2334 mins 

CVbalStool2$accuracy.nvar # check for how many variables gives the highest AUC: 2
CVbalStool2$var.barplot # check for which taxa are included in numerator/denominator
# Denoninator: Desulfovibrio
# Numerator: Subdoligranulum
plot(CVbalStool2$global.plot)
# check how robust the balance variable selection is
plot.tab(CVbalStool2$cv.tab) 
# the global balance (Subdoligranulum/Corynebacterium) is not most frequently slected
# Subdoligranulum/Desulfovibrio is most frequently selected (12%)
summary(CVbalStool2$cv.accuracy) #mean=cross validation accuracy
CVbalStool2$glm






