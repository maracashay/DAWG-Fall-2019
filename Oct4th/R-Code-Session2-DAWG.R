# code developed by DAWG team leaders for Fall 2019 workshops

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("vegan", "RVAideMemoire", "Hmisc", "ecodist", "ggplot2", "FSA", "corrplot", "nlme", "mgcv", "ggpubr", "digest", "backports")
ipak(packages)

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

BiocAnyway <- function(name_of_package){
  if (substr(version$version.strin, 11, 15) >= 3.5){
    install.packages("BiocManager")
    BiocManager::install(name_of_package)}
  else {source('http://bioconductor.org/biocLite.R')
    biocLite(name_of_package)}
}

BiocAnyway("phyloseq")

l.packages <- c("phyloseq", "vegan", "RVAideMemoire", "Hmisc", "ecodist", "ggplot2", "FSA", "corrplot", "nlme", "mgcv", "ggpubr", "digest")

lapply(l.packages, library, character.only = TRUE)

dir <- setwd("E:/")
biom.file <- paste(dir, "otu_table.biom",  sep="")
biom.table <- import_biom(biom.file)
metadata <- read.table("Judie-50-subset-meta.txt", row.names=1, header=TRUE)
map <- sample_data(metadata)


otu.table <- data.frame(otu_table(biom.table))
otu.t <- t(otu.table)
otu <- otu_table(otu.t, taxa_are_rows = FALSE)

ps <- merge_phyloseq(otu, map, tax_table(biom.table))
head(ps@tax_table)
colnames(tax_table(ps)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(ps@tax_table)

ps = subset_taxa(ps, Kingdom %in% c("k__Archaea", "k__Bacteria"))


ntaxa(ps)
# there were 2,367 OTUs in my dataset

ps.1 <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
# this removes all singleton OTUs
ntaxa(ps.1)
# now there are 1,303 OTUs, a decrease of ~ 1,000 OTUs


#### Test for correlations between metadata variables ####
# requires numerical values
metadata
metadata.1 = metadata
metadata.1 <- data.frame(lapply(metadata.1, function(x) as.numeric(as.character(x))))

map <- metadata[,c(2,8,12, 13, 14)]

cor.r <- cor(map, method = "pearson")

corrplot(cor.r, type = "upper", 
         tl.col = "black", tl.srt = 45, na.label="-")


alpha.div<-estimate_richness(ps, measures=c("Shannon", "Simpson"))
alpha.div

map$Shannon <- alpha.div$Shannon
map$Simpson <- alpha.div$Simpson

cor.r <- cor(map, method="pearson")

corrplot(cor.r, type = "upper", tl.col = "black", tl.srt = 45, na.label="-")



ggscatter(map, x = "Age", y = "Death", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coeff.args = list(method="pearson", label.x=55))


ggscatter(map, x = "Shannon", y = "Vent", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coeff.args = list(method="pearson", label.x=1))


#### Heatmaps ####
# the function below allows a user to compress OTUs to a taxonomic level and was developed as part of the microbiome package
taxa_level <- function(physeq,which_level){
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  OTU <- otu_table(physeq)
  SAM <- sample_data(physeq)
  OTU_taxonomy <- tax_table(physeq)
  new_abund_table<-NULL
  if(which_level=="Otus"){
    OTU_tree <- phy_tree(physeq)
    new_abund_table<-OTU
  } else {
    list<-na.omit(unique(OTU_taxonomy[,which_level]))
    new_abund_table<-NULL
    for(i in list){
      rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
      tmp<-data.frame(rowSums(OTU[,rt]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  OTU<-as.data.frame(as(new_abund_table,"matrix"))
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(SAM)
  #reconstruct the phyloseq object
  physeq<-NULL
  if(which_level=="Otus"){
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }
  return(physeq)
}

ps.phylum<- taxa_level(ps.1, "Phylum") #this calculates the abunddances of each phylum across every sample into a new phyloseq object
plot_heatmap(ps.phylum, "PCoA", "bray", "Gender")
plot_heatmap(ps.phylum, "PCoA", "bray", "ARDS")
plot_heatmap(ps.phylum, "PCoA", "bray", "Location")


# Use this to look at genera # 
ps.genera<- taxa_level(ps.1, "Genus")
plot_heatmap(ps.genera, "PCoA", "bray", "Gender")
# this doesn't show us very much because it is too much data 
ps.gen.filtered <- filter_taxa(ps.genera, function(x) sum(x > 1) > (0.50*length(x)), TRUE)
# the above command will keep genera that have more than 1 sequence count in at least 50% of the samples
plot_heatmap(ps.gen.filtered, "PCoA", "bray", "Gender")



plot_heatmap(ps.gen.filtered, "PCoA", "bray", "month")

#### Beta diversity ####
# Bray curtis # 
bc.dist <- phyloseq::distance(otu_table(ps.1), "bray")
adonis(bc.dist ~ ARDS + Gender + inhosp + Location, data = metadata)
#Time is important but there is not an interaction between line and time
# since there is only two time points, there is no need to perform a pairwise comparison

rt.pcoa = ordinate(ps.1, method="PCoA", distance=bc.dist)
plot_scree(rt.pcoa, "Screen plot")
# PCoA ordination
plot_ordination(ps.1, rt.pcoa, "samples", color = "inhosp", shape = "ARDS") + geom_point(size = 7)

ps.1@sam_data$ARDS <- as.factor(ps.1@sam_data$ARDS)
ps.1@sam_data$inhosp <- as.factor(ps.1@sam_data$inhosp)

plot_ordination(ps.1, rt.pcoa, "samples", color = "inhosp", shape = "ARDS") + geom_point(size = 7)


# Identify important variables for explaining Biological variation ####
meta.df <- data.frame(sample_data(ps.1))
otu.table <- data.frame(otu_table(ps.1))

ccamodel <- cca(otu.table~., meta.df, na.action = na.omit)

finalmodel <- ordistep(ccamodel, scope=formula(ccamodel))


finalmodel
# Constrained = 100% 
# Unconstrained = 0%

# some of these variables have VIF > 10 and need to be removed from the dataset for this particular analysis
meta.df$Severe <- NULL
meta.df$Shock <-NULL
meta.df$HAP <- NULL
meta.df$HCAP <- NULL
meta.df$year <- NULL

# make new model
ccamodel <- cca(otu.table~., meta.df, na.action = na.omit)

finalmodel <- ordistep(ccamodel, scope=formula(ccamodel))
finalmodel
# constrained now explains  76.6 % of intertia
# unconstrained explaines 23.4% of inertia

vif.cca(finalmodel)
# some variables are > 10, so need to remove a few more
meta.df$Vent <- NULL
meta.df$inhosp <- NULL
meta.df$ICU <- NULL
ccamodel <- cca(otu.table~., meta.df, na.action = na.omit)

finalmodel <- ordistep(ccamodel, scope=formula(ccamodel))
finalmodel

vif.cca(finalmodel)


mod0<- cca(otu.table ~ 1, meta.df)
mod <- step(mod0, scope=formula(finalmodel), test="perm", perm.max=99)

anova(mod)

#### db - rda ####

ord <- capscale(vegdist(otu.table) ~., meta.df)
ord

anova(ord, by ="margin", perm.max=100)


#### Mantel Tests (EcoDist) ####
xdist.b <- vegdist(otu.table, method = "bray")
ydist.b <- vegdist(meta.df, method = "euclid")
meta.df
map.mantel <- meta.df[,c(2,6,8)]
ydist.b <- vegdist(map.mantel, method = "euclid")

vegan::mantel(xdist.b, ydist.b, method = "spear")


ydist.b <- vegdist(map.mantel$Death, method = "euclid")
vegan::mantel(xdist.b, ydist.b, method = "spear")

ydist.b <- vegdist(map.mantel$APACHE, method = "euclid")
vegan::mantel(xdist.b, ydist.b, method = "spear")


#### Assess Correlations between Continuous variables and Taxa ####
ps.filtered <- filter_taxa(ps.1, function(x) sum(x > 1) > (0.75*length(x)), TRUE)
otu.filtered <- data.frame(otu_table(ps.filtered))

cor.kendall <- cor(otu.filtered, map.mantel, method = "kendall")
cor.kendall
# 3 OTUs have correlation > 0.50 with APACHE


ps.gen.filtered <- filter_taxa(ps.genera, function(x) sum(x > 1) > (0.5*length(x)), TRUE)
gen.filtered <- data.frame(otu_table(ps.gen.filtered))

cor.kendall.gen <- cor(gen.filtered, map.mantel, method = "kendall")
cor.kendall.gen


#### Assessing homogeneity of group dispersions - PermDisp ####
dis <- vegdist(otu.table, "bray")
disp <- betadisper(dis, meta.df$ARDS)
disp

anova(disp)
# since p > 0.05, fail to reject hypothesis that the variances within groups are different between groups

