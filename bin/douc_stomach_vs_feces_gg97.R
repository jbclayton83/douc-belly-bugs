#### Pre-Processing ####
library(vegan)
library(ggplot2)
library(ape)
source('lib/pcoa_helper.R') # for pcoa plots, will have to respond n
library(phyloseq)
library(ggsignif)
library(gplots)
library(RColorBrewer)
library(robCompositions) # Composition magic
library(polycor)
library(beeswarm)
library(reshape2)
library(ltm)

#### DATA WRANGLING ####
map = read.delim('data/gg97/douc_stomach_vs_feces_matching_mapfile_082417.txt', row.names = 1) # Grab the map
### Laod and Filter OTU and Taxa tables with less than 1000 sequences ###
taxa = read.delim('data/gg97/douc_stomach_vs_feces_taxatable_gg97.txt',row=1,as.is=T)
taxa = taxa[,rownames(map)]  # sync the sample names for the Taxa table
otu = read.delim('data/gg97/douc_stomach_vs_feces_otutable_gg97.txt',row=1,as.is=T)
otu = otu[,rownames(map)]  # sync the sample names for the OTU table
dim(otu) # 4502 otus identified across 12 samples
sum(otu) # total reads is under 1 million so singleton threshold is 1
otu <- otu[(rowSums(otu) > 1), ]  # Remove singletons

dim(otu) # 3379 otus identified across 12 samples
depths = colSums(otu)
sort(depths)  # All are good, lowest is 58484.

dim(taxa) # 371 taxa identified across 12 samples
taxa <- taxa[(rowSums(taxa) > 1), ]  # Remove singletons
dim(taxa) # 318 taxa identified across 12 samples
depths.t = colSums(taxa)
sort(depths.t)  # All are good, lowest is 58552.

# Remove rare otus/taxa
#otu = otu[rowMeans(otu > 0) >= 0.1, ]  # Remove rare otu in only one of samples in this case
#dim(otu)
#taxa = taxa[rowMeans(taxa > 0) >= 0.1, ]  # Remove rare taxa, <10% prevalence = if only in one sample in this case
#dim(taxa)
otu.rare <- data.frame(t(rrarefy(t(otu), 58480)))
taxa.rare <- data.frame(t(rrarefy(t(taxa), 58550)))

colSums(taxa.rare)
colSums(otu.rare)

#### F/B Ratio ####
pdf("results/gg97_stomach_feces/ReadDepth_closed_bodysite_stomachvsfeces_gg97.pdf",width=6,height=5.5)
plot(colSums(taxa) ~ map$Bodysite,xlab="Bodysite",ylab="Read depth (closed)")
dev.off()
isFirmicutes = grepl('p__Firmicutes',rownames(taxa))     # Save "trues" for Firmicutes, false otherwise
isBacteroides = grepl('p__Bacteroidetes',rownames(taxa)) # Like above for Bacteroidetes
FBratio = log(colSums(taxa[isFirmicutes,])/colSums(taxa[isBacteroides,])) # Firmicutes/Bacteroidetes log-ratio, vector of lrs named by sampleID

#### Statistical Tests - Bodysite ####
# Independent 2-group Mann Whitney U Test (assumes independent samples)
wilcox.test(FBratio ~ map$Bodysite)
# Paired Wilcoxon (assumes related samples)
wilcox.test(FBratio[map$Bodysite == 'Foregut'], FBratio[map$Bodysite == 'Hindgut'], paired = TRUE, alternative = "two.sided")
pdf("results/gg97_stomach_feces/FBratio_bodysite_stomachvsfeces_gg97.pdf",width=6,height=5.5)
plot(FBratio ~ map$Bodysite, xlab="Body-site", ylab="Log F:B ratio")
dev.off()
df = data.frame(FBratio = colSums(taxa[isFirmicutes,])/colSums(taxa[isBacteroides,]), Bodysite = map$Bodysite)     # Split into groups by bodysite
tapply(FBratio, map$Bodysite, mean)  # Get the mean FB ratio per group
tapply(FBratio, map$Bodysite, sd)    # Gets the standard devs per group
# Plotting Number of Taxa: Bacteroides vs. Firmicutes
# plot(colSums(taxa[isBacteroides,]) ~ map$Bodysite)
# plot(colSums(taxa[isFirmicutes,]) ~ map$Bodysite)
# New Plot of F:B ratio
ggplot(df, aes(x=Bodysite, y=FBratio, fill=Bodysite)) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(pch = 21, alpha = 0.8, stroke = 0.9, width = 0.2) +
  theme_classic() +
  labs(x = "Body-site", y = "Log F:B ratio") +
  guides(fill = guide_legend(title = "Body-site"))
# note: subjects 4 and 5 are driving top quartile of foregut F:B ratio
subj45f <- which(map$Subject %in% c("4a","5a"))
FBratio[subj45f]
colSums(taxa[isFirmicutes,])[subj45f]  # Firmicutes total - subject4 is high here compared to other foreguts
colSums(taxa[isBacteroides,])[subj45f] # Bacteroides total - both lower than expected
mean(colSums(taxa[isFirmicutes,])[-subj45f]) #group mean excluding these two foregut samples
mean(colSums(taxa[isBacteroides,])[-subj45f]) #group mean excluding these two foregut samples


#### PCoA plots - Bodysite #### (requires map and otu table loaded)
tree = read_tree_greengenes('data/gg97/gg97.tre')
otu.s = as.matrix(otu)        # Matrix form of otu table 
bray = vegdist(t(otu.s))      # Get some bray curtis distances
pdf("results/gg97_stomach_feces/bray_bodysite_stomachvsfeces_gg97.pdf",width=6,height=4.75); plot_pcoa(bray,map,category='Bodysite');
pdf("results/gg97_stomach_feces/uuf_bodysite_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.u = plot_unifrac(otu.s,map,tree,category='Bodysite',weight=F); 
pdf("results/gg97_stomach_feces/wuf_bodysite_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.w = plot_unifrac(otu.s,map,tree,category='Bodysite',weight=T); 
graphics.off()
adonis(pcoa.u ~ map$Bodysite)       # Do stats for clustering (unweighted)
adonis(pcoa.w ~ map$Bodysite)       # Do stats for clustering (weighted)
adonis(bray ~ map$Bodysite)         # Do stats for bray-curtis too, why not


#### PCoA plots - Dead or Alive #### (requires map and otu table loaded) 
# tree = read_tree_greengenes('data/gg97/gg97.tre')
# otu.s = as.matrix(otu)      # Matrix form of otu table
# bray = vegdist(t(otu.s))    # Get some bray curtis distances
# pdf("results/gg97_stomach_feces/bray_alive_or_deceased_stomachvsfeces_gg97.pdf",width=6,height=4.75); plot_pcoa(bray,map,category='Alive_or_Deceased_as_of_091814');
# pdf("results/gg97_stomach_feces/uuf_alive_or_deceased_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.u = plot_unifrac(otu.s,map,tree,category='Alive_or_Deceased_as_of_091814',weight=F); 
# pdf("results/gg97_stomach_feces/wuf_alive_or_deceased_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.w = plot_unifrac(otu.s,map,tree,category='Alive_or_Deceased_as_of_091814',weight=T); 
# graphics.off()
# adonis(pcoa.u ~ map$Alive_or_Deceased_as_of_091814)       # Do stats for clustering (unweighted)
# adonis(pcoa.w ~ map$Alive_or_Deceased_as_of_091814)       # Do stats for clustering (weighted)
# adonis(bray ~ map$Alive_or_Deceased_as_of_091814)         # Do stats for bray-curtis too, why not


#### Alpha diversity violin plots - Bodysite ####
(mindepth = min(colSums(otu)))    # Minimum depth is 58484
otu.r = rrarefy(t(otu),mindepth)  # Rarefy to the minimum sample depth - otu.r is transposed
div.shannon = diversity(otu.r,"shannon")   # Shannon index
div.isimp = diversity(otu.r,"invsimpson")  # Simpson (inverse) index
plot(div.isimp~map$Bodysite, xlab="Bodysite",ylab="Shannon Diversity")
plot(div.shannon~map$Bodysite, xlab="Bodysite",ylab="Inverse Simpson Diversity")
plot(rowSums(otu.r > 0) ~ map$Bodysite, xlab="Bodysite",ylab="Number of OTUs")

otu.ad = data.frame(Div=div.shannon, Bodysite=map$Bodysite)
grps = levels(map$Bodysite)
lab = "Alpha Diversity (Shannon)" #paste0("Alpha Diversity (",dix,")")
pdf(paste0("results/gg97_stomach_feces/AlphaDiv_bodysite_doucvsfeces_gg97.pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Bodysite,y=Div,fill=Bodysite)) + ylab(lab) +xlab("Bodysite") + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
       geom_jitter(aes(color=Bodysite),position=position_jitter(0.2),size=2) +
       theme(panel.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size = 14))  )
dev.off()
tapply(otu.ad$Div, otu.ad$Bodysite, mean)    # Gets the mean Shan index per group
tapply(otu.ad$Div, otu.ad$Bodysite, sd)      # Gets the standard devs per group

otu.ad = data.frame(Div=rowSums(otu.r > 0), Bodysite=map$Bodysite)
grps = levels(map$Bodysite)
lab = "Alpha Diversity (Number of OTUs)" 
pdf(paste0("results/gg97_stomach_feces/numOTU_bodysite_stomachvsfeces_gg97.pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Bodysite,y=Div,fill=Bodysite)) + ylab(lab) + xlab("Bodysite") + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
       geom_jitter(aes(color=Bodysite),position=position_jitter(0.2),size=2) +
       theme(panel.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size = 14))  )
dev.off()
tapply(otu.ad$Div, otu.ad$Bodysite, mean)    # Gets the mean num OTUs per group
tapply(otu.ad$Div, otu.ad$Bodysite, sd)      # Gets the standard devs per group


#### Alpha diversity violin plots - Alive or Deceased####
# (mindepth = min(colSums(otu)))
# otu.r = rrarefy(t(otu),mindepth) # Rarefy to the minimum sample depth, otu.r is transposed
# div.shannon = diversity(otu.r,"shannon")
# div.isimp = diversity(otu.r,"invsimpson")
# plot(div.isimp~map$Alive_or_Deceased_as_of_091814, xlab="Alive or Dead",ylab="Shannon Diversity")
# plot(div.shannon~map$Alive_or_Deceased_as_of_091814, xlab="Alive or Dead",ylab="Inverse Simpson Diversity")
# plot(rowSums(otu.r > 0) ~ map$Alive_or_Deceased_as_of_091814, xlab="Alive or Dead",ylab="Number of OTUs")
# 
# otu.ad = data.frame(Div=div.shannon, Status=map$Alive_or_Deceased_as_of_091814)
# grps = levels(map$Alive_or_Deceased_as_of_091814)
# lab = "Alpha Diversity (Shannon)" #paste0("Alpha Diversity (",dix,")")
# pdf(paste0("results/gg97_stomach_feces/AlphaDiv_alive_or_deceased_doucvsfeces_gg97.pdf"),width=6,height=5.5)
# plot(ggplot(otu.ad,aes(x=Status,y=Div,fill=Status)) + ylab(lab) +xlab("Alive or Dead") + geom_violin(alpha=0.3) + 
#        geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
#        geom_jitter(aes(color=Status),position=position_jitter(0.2),size=2) +
#        theme(panel.background = element_blank())  )
# dev.off()
# otu.ad = data.frame(Div=rowSums(otu.r > 0), Status=map$Alive_or_Deceased_as_of_091814)
# grps = levels(map$Alive_or_Deceased_as_of_091814)
# lab = "Alpha Diversity (Number of OTUs)" 
# pdf(paste0("results/gg97_stomach_feces/numOTU_alive_or_deceased_stomachvsfeces_gg97.pdf"),width=6,height=5.5)
# plot(ggplot(otu.ad,aes(x=Status,y=Div,fill=Status)) + ylab(lab) + xlab("Alive or Dead") + geom_violin(alpha=0.3) + 
#        geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
#        geom_jitter(aes(color=Status),position=position_jitter(0.2),size=2) +
#        theme(panel.background = element_blank())  )
# dev.off()


#### Differential taxa testing - by Bodysite ####
bT = c(2,4,6,7)  # Levels to consider for the bugs
lscolors = c("brown1","deepskyblue")
for (L in 1:length(bT)) {
  # Massage the taxa names
  split = strsplit(rownames(taxa),";")  # Split by semicolon into levels
  taxaStrings = sapply(split,function(x) paste(x[1:bT[L]],collapse=";"))          # Create collapsed names
  for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)  # Clean tips
  otu.t = rowsum(taxa,taxaStrings)  # Collapse table by name
  
  # Filter out bugs that don't appear in enough samples
  select = rowSums(otu.t > 0) > min(table(map$Bodysite))/2  # reasonable general min
  otu.t = otu.t[select,]                                    # Apply drop mask
  otu.t = otu.t[,rownames(map)]                             # Sync with map samp order
  otu.n = sweep(otu.t,2,colSums(otu.t),'/');                # Normalize to relative abundance
  
  ### CLR transformation ###
  all.otus.c <- t(otu.t); eps <- 0.5
  all.otus.c <- all.otus.c * (1 - rowSums(all.otus.c==0) * eps / rowSums(all.otus.c))
  all.otus.c[all.otus.c == 0] <- eps
  all.otus.c <- sweep(all.otus.c, 1, rowSums(all.otus.c), '/')
  ls <- log(all.otus.c)
  all.otus.c <- t(ls - rowMeans(ls))
  all.otus.c <- all.otus.c[, !is.nan(colSums(all.otus.c))]
  otu.t <- as.data.frame(all.otus.c)
  
  # Go through each taxon and test for significance w/group
  otu.t = as.matrix(otu.t)
  ntax = nrow(otu.t)
  BS = map$Bodysite
  Grp.Pvals=rep(1,ntax)
  Grp.Corrs=rep(0,ntax)
  TT.Pvals=rep(1,ntax)
  for (m.ix in 1:ntax) {  # Loop through all the rows (taxa)
    try({ # Because some correlations may be inadmissable
      # For >2 levels
      #ps = polyserial(otu.t[m.ix,],map$Bodysite,ML=T,std.err = T)
      #if (is.na(pchisq(ps$chisq, ps$df))) next    # Invalid correlation
      #Grp.Corrs[m.ix] = ps$rho                    # Find intensity of correlation
      #Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
      # For 2 levels (Paired Wilcoxon)
      fg <- otu.t[m.ix,map$Bodysite == "Foregut"]  # Only foregut vals- should be ordered by patient
      hg <- otu.t[m.ix,map$Bodysite == "Hindgut"]  # Only hindgut vals- should be ordered by patient
      pw <- wilcox.test(fg, hg, paired = TRUE)     # Paired difference
      Grp.Pvals[m.ix] <- pw$p.value
      Grp.Corrs[m.ix] <- biserial.cor(otu.t[m.ix,], map$Bodysite, use = "complete.obs")
    },silent=T)
    TT.Pvals[m.ix] <- t.test(otu.t[m.ix,map$Bodysite == "Foregut"], otu.t[m.ix,map$Bodysite == "Hindgut"], paired = TRUE, conf.level = 0.95)$p.value  # Paired difference - ttest
  }
  
  ## Taxa barplots -- Top 15 most abundant (t-test sig. + other?)
  otu.m = otu.n
  otu.m = sweep(sqrt(otu.m),2,colSums(sqrt(otu.m)),'/')
  meanAb = apply(otu.m,1,FUN=function(x) tapply(x, map$Bodysite, mean)) # Group mean
  ranked = order(apply(meanAb,2,max),decreasing=T)
  otu.m = otu.m[ranked,]
  otu.m = otu.m[!grepl('c__Chloroplast', rownames(otu.m), fixed = T), ]  # Remove chloroplast if present
  
  # Truncate names to last 2 informative levels
  split = strsplit(rownames(otu.m),";")        # Split by semicolon into levels
  Taxa = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  lim = 25
  if (nrow(otu.m) > lim) Taxa[lim:nrow(otu.m)] = "Other"
  otu.m = rowsum(otu.m,Taxa)
  
  # Sort by average abundance for display
  byAbundance = rownames(otu.m)[order(rowMeans(otu.m),decreasing=T)]
  
  # Acrobatics for ggplot2
  otu.m = data.frame(t(otu.m),check.names=F)      # Flip table
  otu.m$SampleID = rownames(otu.m)  # Add a column for the sample IDs
  map$SampleID = rownames(map)      # Add a column for the sample IDs
  
  # The following separates taxa abundances per sample, then splices in the column of interest
  otu.m = melt(otu.m, id.vars = "SampleID", variable.name = "Taxa", value.name = "RelativeAbundance")
  otu.m = merge(otu.m, map[,c("SampleID","Bodysite")], by="SampleID")
  otu.m$Taxa = factor(otu.m$Taxa,levels=byAbundance,ordered=T) # add Taxa column
  
  ## Plot according to Bodysite, sorted by abundance
  pdf(paste0("results/gg97_stomach_feces/TaxaSummary_stomachvsfeces_gg97_L",bT[L],".pdf"), width=8, height=7) # Make room for legend
  plot(ggplot(otu.m, aes(x = Bodysite, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Bodysite",y="Root Relative Abundance") +
         guides(fill=guide_legend(ncol=1)) +
         theme(panel.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size = 14)) +
         scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
                                    "green4",
                                    "#6A3D9A", # purple
                                    "#FF7F00", # orange
                                    "black","gold1",
                                    "skyblue2","#FB9A99", # lt pink
                                    "palegreen2",
                                    "#CAB2D6", # lt purple
                                    "#FDBF6F", # lt orange
                                    "gray70", "khaki2",
                                    "maroon","orchid1","deeppink1","blue1","steelblue4",
                                    "darkturquoise","green1","yellow4","yellow3",
                                    "darkorange4","brown"))) 
  dev.off()
  
  ## Display differential taxa/stats
  # Adjust for multiple tests, sort by significance
  gpb = Grp.Pvals; wpb = TT.Pvals;
  Grp.Pvals = p.adjust(gpb,method = "BH")
  TT.Pvals = p.adjust(wpb,method = "BH")
  res = data.frame(TT.Pvals, Grp.Pvals, Grp.Corrs,row.names=rownames(otu.t))
  res = res[order(res$TT.Pvals),]
  
  # Add bivariate filter
  sig = 0.05
  selection = res$TT.Pvals < sig
  
  # Display all significant with p < 0.05
  num_sig = sum(selection, na.rm = T) # Count how many are significant
  res = res[selection,]
  
  # Truncate names to last 2 informative levels
  split = strsplit(rownames(res),";")        # Split by semicolon into levels
  res$short = sapply(split,function(x) paste(tail(x,2),collapse=";"))
  
  pdf(paste0("results/gg97_stomach_feces/TaxaSwarms_stomachvsfeces_gg97_L",bT[L],".pdf"),width = 6.5,height=6.5)
  sink(paste0("results/gg97_stomach_feces/Taxa_Significance_stomachvsfeces_gg97_L",bT[L],".txt"))   # Get ready to write the significant ones
  cat("Taxon\tPWilcoxon_P\tBiserial_Cor\tBodySite_P\n")  # Print header
  if (num_sig) for (i in 1:num_sig) {
    taxon = rownames(res)[i]
    cat(res[taxon,]$short,'\t',res$Grp.Pvals[i],'\t',-res$Grp.Corrs[i],'\t',res$TT.Pvals[i],'\t','\n',sep='')
    beeswarm(otu.t[taxon,] ~ map$Bodysite, xlab="Bodysite",ylab="CLR Relative Abundance",main=res[taxon,]$short,
             col=alpha(lscolors,0.7),
             cex.axis=1.3,cex.main=1.4,cex.lab=1.3,cex=1.1,corral="random",pch=16)
    bxplot(otu.t[taxon,] ~ map$Bodysite, add = TRUE)
  }
  sink(NULL)
  dev.off()
  
  ## Heatmap ##
  if (num_sig < 2) next
  # Need to have created the clr taxa table as a matrix
  my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299) # Sebastian Raschka
  gl = map$Bodysite
  glpos = c(grep("Foregut",gl),grep("Hindgut",gl))
  gl = gl[glpos]
  mat = otu.t[rownames(res),glpos]
  
  # Truncate names to last 2 informative levels
  split = strsplit(as.character(rownames(mat)),";")        # Split by semicolon into levels
  rownames(mat) = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  
  levels(gl)= c("brown1","deepskyblue") # lscolors
  png(paste0("results/gg97_stomach_feces/Taxa_heatmap_stomachvsfeces_gg97_L",bT[L],".png"),  # Create PNG for the heat map        
      width = 8*300,                          # 8 x 300 pixels
      height = 6*300,
      res = 300,                              # 300 pixels per inch
      pointsize = 10)                         # Smaller font size
  heatmap.2(mat,
            #cellnote = mat,       # Same dataset for cell labels
            #notecol="black",      # Change font color of cell labels to black
            main = "", # heat map title
            density.info="none",   # Turns off density plot inside color legend
            trace="none",          # Turns off trace lines inside the heat map
            margins = c(2,22),     # Widens margins around plot
            col=my_palette,        # Use on color palette defined earlier
            #breaks=col_breaks,    # Enable color transition at specified limits
            ColSideColors = as.character(gl),
            dendrogram="row",      # Only draw a row dendrogram
            lhei=c(1,4.7), lwid=c(1,5),
            labCol = "",
            cexRow = 0.80,         # Change font size of row labels
            hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete"),
            Colv="NA"              # Turn off column clustering
  )
  par(lend = 1)           # Square line ends for the color legend
  legend("topright",      # Location of the legend on the heatmap plot
         inset=c(.1,-0),  # Adjust placement upward
         legend = levels(map$Bodysite), # Category labels
         col = levels(gl),  # Color key
         lty= 1,            # Line style
         lwd = 10,          # Line width
         cex = 0.65,
         xpd=TRUE  # allow drawing outside
  )
  dev.off()
}


#### Differential taxa testing - by Subject ####
bT = c(2,4,6,7)  # Levels to consider for the bugs
lscolors = c("brown1","deepskyblue")
for (L in 1:length(bT)) {
  # Massage the taxa names
  split = strsplit(rownames(taxa),";")  # Split by semicolon into levels
  taxaStrings = sapply(split,function(x) paste(x[1:bT[L]],collapse=";"))           # Create collapsed names
  for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)   # Clean tips
  otu.t = rowsum(taxa,taxaStrings) # Collapse table by name
  
  # Filter out bugs that don't appear in enough samples
  select = rowSums(otu.t > 0) > min(table(map$Subject))/2 # reasonable general min
  otu.t = otu.t[select,]                                    # Apply drop mask
  otu.t = otu.t[,rownames(map)]                             # Sync with map samp order
  otu.n = sweep(otu.t,2,colSums(otu.t),'/');                # Normalize to relative abundance
  
  ### CLR transformation ###
  all.otus.c <- t(otu.t); eps <- 0.5
  all.otus.c <- all.otus.c * (1 - rowSums(all.otus.c==0) * eps / rowSums(all.otus.c))
  all.otus.c[all.otus.c == 0] <- eps
  all.otus.c <- sweep(all.otus.c, 1, rowSums(all.otus.c), '/')
  ls <- log(all.otus.c)
  all.otus.c <- t(ls - rowMeans(ls))
  all.otus.c <- all.otus.c[, !is.nan(colSums(all.otus.c))]
  otu.t <- as.data.frame(all.otus.c)
  
  # Calculate Foregut/Hindgut Relationships (FF, HH, FH) - Distances/Permutations
  # We want to measure the observed distances among all foreguts, among all hindguts,
  # and between foregut/hindgut for each individual. To validate these values, we can
  # conduct a permutation of these distances (changing the labels in the otu tree).
  otu.d <- as.matrix(otu)
  phyobj <- phyloseq(otu_table(otu.d, taxa_are_rows=T), tree)
  obs <- UniFrac(phyobj, weighted=F) # observed unweighted UniFrac distances
  obs.m <- as.matrix(obs)  # view as matrix for calculating sums
  obs.ff <- sum(obs.m[seq(from=1, to=12, by=2), seq(from=2, to=12, by=2)])/30 # foreguts across all subjects (mean)
  obs.hh <- sum(obs.m[seq(from=2, to=12, by=2), seq(from=2, to=12, by=2)])/30 # hindguts across all subjects (mean)
  obs.fh <- sum(c(obs.m[2,1], obs.m[4,3], obs.m[6,5], obs.m[8,7], obs.m[10,9], obs.m[12,11]))/6 # foregut-hindgut for each subject (mean)
  ### Foregut/Hindgut
  adonis(obs ~ map$Bodysite, permutations=999)   # permutation to determine validity of distances
  
  ## Taxa barplots -- Top 15 most abundant
  otu.m = otu.n # Normalized
  otu.m = sweep(sqrt(otu.m),2,colSums(sqrt(otu.m)),'/')
  meanAb = apply(otu.m,1,FUN=function(x) tapply(x, map$Subject, mean)) # Group mean
  ranked = order(apply(meanAb,2,max),decreasing=T)
  otu.m = otu.m[ranked,]
  
  # Truncate names to last 2 informative levels
  split = strsplit(rownames(otu.m),";")        # Split by semicolon into levels
  Taxa = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  lim = 25
  if (nrow(otu.m) > lim) Taxa[lim:nrow(otu.m)] = "Other"
  otu.m = rowsum(otu.m,Taxa)
  
  # Sort by average abundance for display
  byAbundance = rownames(otu.m)[order(rowMeans(otu.m),decreasing=T)]
  
  # Acrobatics for ggplot2 (yes, this is inane)
  otu.m = data.frame(t(otu.m),check.names=F)      # flip table
  otu.m$SampleID = rownames(otu.m)  # add a column for the sample IDs
  map$SampleID = rownames(map)      # add a column for the sample IDs
  
  # The following separates taxa abundances per sample, then splices in the column of interest
  otu.m = melt(otu.m, id.vars = "SampleID", variable.name = "Taxa", value.name = "RelativeAbundance")
  otu.m = merge(otu.m, map[,c("SampleID","Subject")], by="SampleID")
  otu.m$Taxa = factor(otu.m$Taxa,levels=byAbundance,ordered=T) # add Taxa column
  otu.m$Subject = gsub("a", "f", otu.m$Subject) # change name to #f for foregut
  otu.m$Subject = gsub("b", "h", otu.m$Subject) # change name to #h for hindgut
  
  ## Plot according to Subject, sorted by abundance
  pdf(paste0("results/gg97_stomach_feces/TaxaSummary_stomachvsfeces_subject_gg97_L",bT[L],".pdf"),width = 8,height=7) # Make room for legend
  plot(ggplot(otu.m, aes(x = Subject, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Subject",y="Root Relative Abundance", title="Relative Abundance by Subject") +
         guides(fill=guide_legend(ncol=1)) +
         theme(panel.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size = 14)) +
         scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
                                    "green4",
                                    "#6A3D9A", # purple
                                    "#FF7F00", # orange
                                    "black","gold1",
                                    "skyblue2","#FB9A99", # lt pink
                                    "palegreen2",
                                    "#CAB2D6", # lt purple
                                    "#FDBF6F", # lt orange
                                    "gray70", "khaki2",
                                    "maroon","orchid1","deeppink1","blue1","steelblue4",
                                    "darkturquoise","green1","yellow4","yellow3",
                                    "darkorange4","brown"))) 
  dev.off()
  
  ## Plot according to Subject, sorted by abundance (Foregut only)
  otu.mf <- otu.m[grep("f", otu.m$Subject),]
  otu.mf$Subject <- gsub("f", "", otu.mf$Subject)
  pdf(paste0("results/gg97_stomach_feces/TaxaSummary_stomachvsfeces_subject_gg97_L",bT[L],"_F.pdf"),width = 8,height=7) # Make room for legend
  plot(ggplot(otu.mf, aes(x = Subject, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Subject",y="Root Relative Abundance", title="Foregut Relative Abundance") +
         guides(fill=guide_legend(ncol=1)) +
         theme(panel.background = element_blank(), axis.text = element_text(size=14), axis.title = element_text(size = 16), plot.title = element_text(size = 18)) +
         scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
                                    "green4",
                                    "#6A3D9A", # purple
                                    "#FF7F00", # orange
                                    "black","gold1",
                                    "skyblue2","#FB9A99", # lt pink
                                    "palegreen2",
                                    "#CAB2D6", # lt purple
                                    "#FDBF6F", # lt orange
                                    "gray70", "khaki2",
                                    "maroon","orchid1","deeppink1","blue1","steelblue4",
                                    "darkturquoise","green1","yellow4","yellow3",
                                    "darkorange4","brown"))) 
  dev.off()
  
  ## Plot according to Subject, sorted by abundance (Hindgut only)
  otu.mh <- otu.m[grep("h", otu.m$Subject),]
  otu.mh$Subject <- gsub("h", "", otu.mh$Subject)
  pdf(paste0("results/gg97_stomach_feces/TaxaSummary_stomachvsfeces_subject_gg97_L",bT[L],"_H.pdf"),width = 8,height=7) # Make room for legend
  plot(ggplot(otu.mh, aes(x = Subject, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Subject",y="Root Relative Abundance", title="Hindgut Relative Abundance") +
         guides(fill=guide_legend(ncol=1)) +
         theme(panel.background = element_blank(), axis.text = element_text(size=14), axis.title = element_text(size = 16), plot.title = element_text(size = 18)) +
         scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
                                    "green4",
                                    "#6A3D9A", # purple
                                    "#FF7F00", # orange
                                    "black","gold1",
                                    "skyblue2","#FB9A99", # lt pink
                                    "palegreen2",
                                    "#CAB2D6", # lt purple
                                    "#FDBF6F", # lt orange
                                    "gray70", "khaki2",
                                    "maroon","orchid1","deeppink1","blue1","steelblue4",
                                    "darkturquoise","green1","yellow4","yellow3",
                                    "darkorange4","brown"))) 
  dev.off()
  
  ## Heatmap ##
  # if (num_sig < 2) next
  # Need to have created the clr taxa table as a matrix
  my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299) # Sebastian Raschka
  gl = map$Subject
  glpos = c(grep("1a",gl),grep("1b",gl),grep("2a",gl),grep("2b",gl),grep("3a",gl),grep("3b",gl),grep("4a",gl),grep("4b",gl),grep("5a",gl),grep("5b",gl),grep("6a",gl),grep("6b",gl))
  gl = gl[glpos]
  mat = otu.t[,glpos] #could further restrict here to significant taxa only
  
  # Truncate names to last 2 informative levels
  split = strsplit(as.character(rownames(mat)),";")        # Split by semicolon into levels
  rownames(mat) = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  
  levels(gl)= rep(c("brown1","deepskyblue"), 6) # lscolors
  png(paste0("results/gg97_stomach_feces/Taxa_heatmap_stomachvsfeces_subject_gg97_L",bT[L],".png"),  # create PNG for the heat map        
      width = 8*300,                        # 8 x 300 pixels
      height = 6*300,
      res = 300,                              # 300 pixels per inch
      pointsize = 10)                          # smaller font size
  heatmap.2(as.matrix(mat),
            #cellnote = mat,  # same data set for cell labels
            main = "", # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins = c(2,22),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            ColSideColors = as.character(gl),
            dendrogram="row",     # only draw a row dendrogram
            lhei=c(1,4.7), lwid=c(1,5),
            labCol = "",
            hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete"),
            Colv="NA"            # turn off column clustering
  )
  par(lend = 1)           # square line ends for the color legend
  legend("topright",      # location of the legend on the heatmap plot
         inset=c(.1,-0), # adjust placement upward
         legend = levels(map$Subject), # category labels
         col = levels(gl),  # color key
         lty= 1,            # line style
         lwd = 10,          # line width
         cex = 0.65,
         xpd=TRUE  # allow drawing outside
  )
  dev.off()
}


#### PICRUST ####
# Read in the PICRUSt L3 summarized pathways (stage 3 output)
picrust = read.delim('data/gg97/douc_stomach_vs_feces_otutable_gg97_predictions_categorized_L3.txt',
                     skip=1, row.names = 1) #Grab picrust table, skipping bad first row
picrust = as.matrix(picrust[,rownames(map)]) # sync and drop extra

# Convert to relative abundance - CLR
picrust = t(picrust); eps = 0.2
picrust = picrust*(1 - rowSums(picrust==0)*eps/rowSums(picrust))
picrust[picrust==0]=eps
picrust = sweep(picrust,1,rowSums(picrust),'/');
ls = log(picrust)
picrust = t(ls - rowMeans(ls))

# CLR with simple substitutions of zeros
#picrust[picrust==0]=0.5
#picrust = sweep(picrust,2,colSums(picrust),'/')
#picrust.clr = cenLR(t(picrust))$x.clr
#picrust = t(picrust.clr)

# Just relative abundance (no CLR)
#picrust = sweep(picrust,2,colSums(picrust),'/')
#picrust = sweep(sqrt(picrust),2,colSums(sqrt(picrust)),'/')

# Go through each picrust pathway and test for significance w/group
npaths = nrow(picrust)
BS = map$Bodysite
Grp.Pvals=rep(1,npaths)
Grp.Corrs=rep(0,npaths)
TT.Pvals=rep(1,npaths)
for (m.ix in 1:npaths) {  # Loop through all the rows (taxa)
  try({ # Because some correlations may be inadmissable
    #ps = polyserial(picrust[m.ix,],map$Bodysite,ML=T,std.err = T)
    #if (is.na(pchisq(ps$chisq, ps$df))) next # Invalid correlation
    #Grp.Corrs[m.ix] = ps$rho             # Find intensity of correlation
    #Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
    # Paired Wilcoxon
    fg <- picrust[m.ix,map$Bodysite == "Foregut"]
    hg <- picrust[m.ix,map$Bodysite == "Hindgut"]
    pw <- wilcox.test(fg, hg, paired = TRUE)
    Grp.Pvals[m.ix] <- pw$p.value
    Grp.Corrs[m.ix] <- biserial.cor(picrust[m.ix,], map$Bodysite, use = "complete.obs")
  },silent=T)
  TT.Pvals[m.ix] = t.test(picrust[m.ix,map$Bodysite == "Foregut"], picrust[m.ix,map$Bodysite == "Hindgut"], paired = TRUE, conf.level = 0.95)$p.value
}

# Adjust for multiple tests
Grp.Pvals = p.adjust(Grp.Pvals, method = "fdr")

TT.Pvals = p.adjust(TT.Pvals, method = "fdr")
res = data.frame(TT.Pvals, Grp.Pvals, Grp.Corrs,row.names=rownames(picrust))
res = res[order(res$TT.Pvals),]

# Add bivariate filter
sig = 0.05
selection = res$TT.Pvals < sig

# Display all significant with p < 0.05
num_sig = sum(selection, na.rm = T) # Count how many are significant
res = res[selection,]

pdf("results/gg97_stomach_feces/PicrustSwarms_gg97.pdf",width = 6.5,height=6.5)
sink("results/gg97_stomach_feces/Picrust_Significance_gg97.txt")   # Get ready to write the significant ones
cat("Pathway\tPWilcoxon_Q\tBiserial_Cor\tBodySite_Q\n")  # Print header

if (num_sig) for (i in 1:num_sig) {
  pathway = rownames(res)[i]
  cat(pathway,'\t',res$Grp.Pvals[i],'\t',-res$Grp.Corrs[i],'\t',res$TT.Pvals[i],'\n',sep='')
  beeswarm(picrust[pathway,] ~ map$Bodysite, xlab="Bodysite",ylab="Pathway Abundance",main=pathway,
           col=alpha(lscolors,0.7),
           cex.axis=1.3,cex.main=1.4,cex.lab=1.3,cex=1.1,corral="random",pch=16)
  bxplot(picrust[pathway,] ~ map$Bodysite, add = TRUE)
}
sink(NULL)
dev.off()

# PICRUSt heatmap too, why not
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299) # Sebastian Raschka
gl = map$Bodysite
glpos = c(grep("Foregut",gl),grep("Hindgut",gl))
gl = gl[glpos]
mat = picrust[rownames(res[abs(res$Grp.Corrs) > 0.75,]),glpos]
mat = sweep(mat,1,rowSums(abs(mat)),'/')                      # Normalize to relative abundance
mat = sweep(mat,1,max(apply(mat,1,max),apply(mat,1,min)),'/') # Constrain extrema to [-1, 1]

levels(gl)= lscolors #c("brown1", "deepskyblue")

png("results/gg97_stomach_feces/PiMap_gg97.png",  # create PNG for the heat map        
    width = 8*300,                        # 5 x 300 pixels
    height = 6*300,
    res = 300,                              # 300 pixels per inch
    pointsize = 11)                          # smaller font size
heatmap.2(mat,
          #cellnote = mat,  # same data set for cell labels
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(2,22),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          ColSideColors = as.character(gl),
          dendrogram="row",     # only draw a row dendrogram
          lhei=c(1,4), lwid=c(1,4),
          labCol = "",
          hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete"),
          Colv="NA"            # turn off column clustering
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       inset=c(.1,-0), # adjust placement upward
       legend = levels(map$Bodysite), # category labels
       col = levels(gl),  # color key
       lty= 1,            # line style
       lwd = 10,          # line width
       cex = 0.75,
       xpd=TRUE  # allow drawing outside
)
dev.off()

