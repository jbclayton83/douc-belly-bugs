#### DATA WRANGLING ####
#setwd('/data/analyses/monkey/') # Sets the working directory
map = read.delim('data/gg97/douc_stomach_vs_feces_matching_mapfile_082417.txt', row.names = 1) # Grab the map
# map$Bodysite = factor(levels = c("Foregut","Hindgut"),ordered = T)

### Laod and Filter OTU and Taxa tables with less than 1000 sequences ###
taxa = read.delim('data/gg97/douc_stomach_vs_feces_taxatable_gg97.txt',row=1,as.is=T)
taxa = taxa[,rownames(map)]  # sync the sample names for the Taxa table
otu = read.delim('data/gg97/douc_stomach_vs_feces_otutable_gg97.txt',row=1,as.is=T)
otu = otu[,rownames(map)]  # sync the sample names for the OTU table


#### F/B Ratio ####
pdf("results/gg97_stomach_feces/ReadDepth_closed_bodysite_stomachvsfeces_gg97.pdf",width=6,height=5.5)
plot(colSums(taxa) ~ map$Bodysite,xlab="Bodysite",ylab="Read depth (closed)")
dev.off()
isFirmicutes = grepl('p__Firmicutes',rownames(taxa)) # Save "trues" for Firmicutes, false otherwise
isBacteroides = grepl('p__Bacteroidetes',rownames(taxa)) # Like above for Bacteroidetes
FBratio = log(colSums(taxa[isFirmicutes,])/colSums(taxa[isBacteroides,])) # Firmicutes/Bacteroidetes log-ratio, vector of lrs named by sampleID

# Kruskal Test - nonparametric, one-way ANOVA with rank (unnecessary for only two levels)
#kruskal.test(FBratio ~ map$Bodysite)
# 1 Sample Wilcoxon - nonparametric, similar to 1 sample t-test (loop is unnecessary for two levels)
#for (i in 1:length(levels(map$Bodysite))) for (j in i:length(levels(map$Bodysite))) {
  #if (i == j) next
  #p = wilcox.test(FBratio[map$Bodysite==levels(map$Bodysite)[i]],FBratio[map$Bodysite==levels(map$Bodysite)[j]])
  #cat(levels(map$Bodysite)[i]," vs ",levels(map$Bodysite)[j]," p = ",p$p.value,"\n",sep='')
#}
# Independent 2-group Mann Whitney U Test (assumes independent samples)
wilcox.test(FBratio ~ map$Bodysite)
# Paired Wilcoxon (this test assumes related samples, as in multiple measurements from same patients)
wilcox.test(FBratio[map$Bodysite == 'Foregut'], FBratio[map$Bodysite == 'Hindgut'], paired = TRUE, alternative = "two.sided")
pdf("results/gg97_stomach_feces/FBratio_bodysite_stomachvsfeces_gg97.pdf",width=6,height=5.5)
plot(FBratio ~ map$Bodysite, xlab="Bodysite", ylab="Log F:B ratio")
dev.off()
df = data.frame(FBratio, map$Bodysite)     # Split manually into groups
tapply(FBratio, map$Bodysite, mean)  # Get the means per group
tapply(FBratio, map$Bodysite, sd)    # Gets the standard devs per group

#plot(colSums(taxa[isBacteroides,]) ~ map$Bodysite)
#plot(colSums(taxa[isFirmicutes,]) ~ map$Bodysite)


#### PCoA plots - Bodysite #### (requires map and otu table loaded) 
source('lib/pcoa_helper.R') # This gives us our nice pcoa functions
library(phyloseq)
tree = read_tree_greengenes('data/gg97/gg97.tre')
otu.s = as.matrix(otu)  # Rip of the taxonomy column, as it is not needed 
bray = vegdist(t(otu.s))             # Get some bray curtis distances
pdf("results/gg97_stomach_feces/bray_bodysite_stomachvsfeces_gg97.pdf",width=6,height=4.75); plot_pcoa(bray,map,category='Bodysite');
pdf("results/gg97_stomach_feces/uuf_bodysite_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.u = plot_unifrac(otu.s,map,tree,category='Bodysite',weight=F); 
pdf("results/gg97_stomach_feces/wuf_bodysite_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.w = plot_unifrac(otu.s,map,tree,category='Bodysite',weight=T); 
graphics.off()
adonis(pcoa.u ~ map$Bodysite)       # Do stats for clustering (unweighted)
adonis(pcoa.w ~ map$Bodysite)       # Do stats for clustering (weighted)
adonis(bray ~ map$Bodysite)         # Do stats for bray-curtis too, why not


#### PCoA plots - Dead or Alive #### (requires map and otu table loaded) 
source('lib/pcoa_helper.R') # This gives us our nice pcoa functions
library(phyloseq)
tree = read_tree_greengenes('data/gg97/gg97.tre')
otu.s = as.matrix(otu)  # Rip of the taxonomy column, as it is not needed 
bray = vegdist(t(otu.s))             # Get some bray curtis distances
pdf("results/gg97_stomach_feces/bray_alive_or_deceased_stomachvsfeces_gg97.pdf",width=6,height=4.75); plot_pcoa(bray,map,category='Alive_or_Deceased_as_of_091814');
pdf("results/gg97_stomach_feces/uuf_alive_or_deceased_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.u = plot_unifrac(otu.s,map,tree,category='Alive_or_Deceased_as_of_091814',weight=F); 
pdf("results/gg97_stomach_feces/wuf_alive_or_deceased_stomachvsfeces_gg97.pdf",width=6,height=4.75); pcoa.w = plot_unifrac(otu.s,map,tree,category='Alive_or_Deceased_as_of_091814',weight=T); 
graphics.off()
adonis(pcoa.u ~ map$Alive_or_Deceased_as_of_091814)       # Do stats for clustering (unweighted)
adonis(pcoa.w ~ map$Alive_or_Deceased_as_of_091814)       # Do stats for clustering (weighted)
adonis(bray ~ map$Alive_or_Deceased_as_of_091814)         # Do stats for bray-curtis too, why not
# ~Dropbox/douc-belly-bugs/src/

#### Alpha diversity violin plots - Bodysite ####
library(vegan)
(mindepth = min(colSums(otu))) # for row, -grep("Chloroplast",otu$taxonomy)?
otu.r = rrarefy(t(otu),mindepth) # Rarefy to the minimum sample depth
div.shannon = diversity(otu.r,"shannon")
div.isimp = diversity(otu.r,"invsimpson")
plot(div.isimp~map$Bodysite, xlab="Body Site",ylab="Shannon Diversity")
plot(div.shannon~map$Bodysite, xlab="Body Site",ylab="Inverse Simpson Diversity")
plot(rowSums(otu.r > 0) ~ map$Bodysite, xlab="Body Site",ylab="Number of OTUs")

library(ggsignif)
otu.ad = data.frame(Div=div.shannon, Body_site=map$Bodysite)
grps = levels(map$Bodysite)
lab = "Alpha Diversity (Shannon)" #paste0("Alpha Diversity (",dix,")")
pdf(paste0("results/gg97_stomach_feces/AlphaDiv_bodysite_doucvsfeces_gg97.pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Body_site,y=Div,fill=Body_site)) + ylab(lab) +xlab("Body Site") + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
       geom_jitter(aes(color=Body_site),position=position_jitter(0.2),size=2) +
       theme(panel.background = element_blank())  )
dev.off()
otu.ad = data.frame(Div=rowSums(otu.r > 0), Body_site=map$Bodysite)
grps = levels(map$Bodysite)
lab = "Alpha Diversity (Number of OTUs)" 
pdf(paste0("results/gg97_stomach_feces/numOTU_bodysite_stomachvsfeces_gg97.pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Body_site,y=Div,fill=Body_site)) + ylab(lab) + xlab("Body Site") + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
       geom_jitter(aes(color=Body_site),position=position_jitter(0.2),size=2) +
       theme(panel.background = element_blank())  )
dev.off()


#### Alpha diversity violin plots - Alive or Deceased####
library(vegan)
(mindepth = min(colSums(otu))) # for row, -grep("Chloroplast",otu$taxonomy)?
otu.r = rrarefy(t(otu),mindepth) # Rarefy to the minimum sample depth
div.shannon = diversity(otu.r,"shannon")
div.isimp = diversity(otu.r,"invsimpson")
plot(div.isimp~map$Alive_or_Deceased_as_of_091814, xlab="Alive or Dead",ylab="Shannon Diversity")
plot(div.shannon~map$Alive_or_Deceased_as_of_091814, xlab="Alive or Dead",ylab="Inverse Simpson Diversity")
plot(rowSums(otu.r > 0) ~ map$Alive_or_Deceased_as_of_091814, xlab="Alive or Dead",ylab="Number of OTUs")

library(ggsignif)
otu.ad = data.frame(Div=div.shannon, Status=map$Alive_or_Deceased_as_of_091814)
grps = levels(map$Alive_or_Deceased_as_of_091814)
lab = "Alpha Diversity (Shannon)" #paste0("Alpha Diversity (",dix,")")
pdf(paste0("results/gg97_stomach_feces/AlphaDiv_alive_or_deceased_doucvsfeces_gg97.pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Status,y=Div,fill=Status)) + ylab(lab) +xlab("Alive or Dead") + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
       geom_jitter(aes(color=Status),position=position_jitter(0.2),size=2) +
       theme(panel.background = element_blank())  )
dev.off()
otu.ad = data.frame(Div=rowSums(otu.r > 0), Status=map$Alive_or_Deceased_as_of_091814)
grps = levels(map$Alive_or_Deceased_as_of_091814)
lab = "Alpha Diversity (Number of OTUs)" 
pdf(paste0("results/gg97_stomach_feces/numOTU_alive_or_deceased_stomachvsfeces_gg97.pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Status,y=Div,fill=Status)) + ylab(lab) + xlab("Alive or Dead") + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)]), test='t.test', map_signif_level = T) + 
       geom_jitter(aes(color=Status),position=position_jitter(0.2),size=2) +
       theme(panel.background = element_blank())  )
dev.off()


#### Differential taxa testing - by Bodysite ####
library(gplots)
library(RColorBrewer)
library(robCompositions) # Composition magic
library(polycor)
library(beeswarm)
library(reshape2)
library(ltm)

bT = c(2,4,6,7)  # Levels to consider for the bugs
lscolors = c("red","green")
for (L in 1:length(bT)) {
  # Massage the taxa names
  split = strsplit(rownames(taxa),";")  # Split by semicolon into levels
  taxaStrings = sapply(split,function(x) paste(x[1:bT[L]],collapse=";"))           # Create collapsed names
  for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)   # Clean tips
  otu.t = rowsum(taxa,taxaStrings) # Collapse table by name
  
  # Filter out bugs that don't appear in enough samples
  select = rowSums(otu.t > 0) > min(table(map$Bodysite))/2 # reasonable general min
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
  
  # otu.c = impRZilr(t(otu.t)+0.0,dl=rep(1,nrow(otu.t)),maxit = 3,verbose = T,method = "lm") # zeros
  # otu.c = t(cenLR(otu.c$x)$x.clr)    # Centered log-ratio transform for compositions
  # colnames(otu.c) = colnames(otu.t)  # Because this gets rid of the names...
  # otu.t = otu.c                      # otu.t is our active table; give it the CLR
  
  # Go through each taxon and test for significance w/group
  otu.t = as.matrix(otu.t)
  ntax = nrow(otu.t)
  BS = map$Bodysite
  Grp.Pvals=rep(1,ntax)
  Grp.Corrs=rep(0,ntax)
  #KW.Pvals=rep(1,ntax)
  TT.Pvals=rep(1,ntax)
  for (m.ix in 1:ntax) {  # Loop through all the rows (taxa)
    try({ # Because some correlations may be inadmissable
      # For >2 levels
      #ps = polyserial(otu.t[m.ix,],map$Bodysite,ML=T,std.err = T)
      #if (is.na(pchisq(ps$chisq, ps$df))) next # Invalid correlation
      #Grp.Corrs[m.ix] = ps$rho             # Find intensity of correlation
      #Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
      # For 2 levels (Paired Wilcoxon)
      fg <- otu.t[m.ix,map$Bodysite == "Foregut"] #only foregut vals- should be ordered by patient
      hg <- otu.t[m.ix,map$Bodysite == "Hindgut"] #only hindgut vals- should be ordered by patient
      pw <- wilcox.test(fg, hg, paired = TRUE) #just looks for difference, can also test "greater" or "less"
      Grp.Pvals[m.ix] <- pw$p.value
      Grp.Corrs[m.ix] <- biserial.cor(otu.t[m.ix,], map$Bodysite, use = "complete.obs")
    },silent=T)
    #KW.Pvals[m.ix] = kruskal.test(otu.t[m.ix,] ~ map$Bodysite)$p.val
    TT.Pvals[m.ix] <- t.test(otu.t[m.ix,map$Bodysite == "Foregut"], otu.t[m.ix,map$Bodysite == "Hindgut"], paired = TRUE, conf.level = 0.95)$p.value
    # add Chi-Squared Test as well (? - assumes independence of variables)
  }
  
  ## Taxa barplots -- Top 15 most abundant (t-test sig. + other?)
  otu.m = otu.n
  otu.m = sweep(sqrt(otu.m),2,colSums(sqrt(otu.m)),'/')
  meanAb = apply(otu.m,1,FUN=function(x) tapply(x, map$Bodysite, mean)) # group mean
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
  otu.m = merge(otu.m, map[,c("SampleID","Bodysite")], by="SampleID")
  otu.m$Taxa = factor(otu.m$Taxa,levels=byAbundance,ordered=T) # add Taxa column
  
  ## Plot according to Bodysite, sorted by abundance
  pdf(paste0("results/gg97_stomach_feces/TaxaSummary_stomachvsfeces_gg97_L",bT[L],".pdf"),width = 8,height=7) # Make room for legend
  plot(ggplot(otu.m, aes(x = Bodysite, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Bodysite",y="Root Relative Abundance") +
         guides(fill=guide_legend(ncol=1)) +
         theme(panel.background = element_blank()) +
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
    beeswarm(otu.t[taxon,] ~ map$Bodysite, xlab="Body Site",ylab="CLR Relative Abundance",main=res[taxon,]$short,
             col=alpha(lscolors,0.7),
             cex.axis=1.1,cex.main=1,cex=1.1,corral="random",pch=19)
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
  
  levels(gl)= c("red","green") # lscolors
  png(paste0("results/gg97_stomach_feces/Taxa_heatmap_stomachvsfeces_gg97_L",bT[L],".png"),  # create PNG for the heat map        
      width = 8*300,                        # 8 x 300 pixels
      height = 6*300,
      res = 300,                              # 300 pixels per inch
      pointsize = 10)                          # smaller font size
  heatmap.2(mat,
            #cellnote = mat,  # same data set for cell labels
            #notecol="black",      # change font color of cell labels to black
            main = "", # heat map title
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins = c(2,22),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            ColSideColors = as.character(gl),
            dendrogram="row",     # only draw a row dendrogram
            lhei=c(1,4.7), lwid=c(1,5),
            labCol = "",
            cexRow = 0.80,         # change font size of row labels
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
         cex = 0.65,
         xpd=TRUE  # allow drawing outside
  )
  dev.off()
}


#### Differential taxa testing - by Subject ####
library(gplots)
library(RColorBrewer)
library(robCompositions) # Composition magic
library(polycor)
library(beeswarm)
library(reshape2)

bT = c(2,4,6,7)  # Levels to consider for the bugs
lscolors = c("red","green")
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
  
  # otu.c = impRZilr(t(otu.t)+0.0,dl=rep(1,nrow(otu.t)),maxit = 3,verbose = T,method = "lm") # zeros
  # otu.c = t(cenLR(otu.c$x)$x.clr)    # Centered log-ratio transform for compositions
  # colnames(otu.c) = colnames(otu.t)  # Because this gets rid of the names...
  # otu.t = otu.c                      # otu.t is our active table; give it the CLR
  
  # Go through each taxon and test for significance w/group
  otu.t = as.matrix(otu.t)
  ntax = nrow(otu.t)
  BS = map$Subject
  Grp.Pvals=rep(1,ntax)
  Grp.Corrs=rep(0,ntax)
  #KW.Pvals=rep(1,ntax)
  TT.Pvals=rep(1,ntax)
  for (m.ix in 1:ntax) {  # Loop through all the rows (taxa)
    try({ # Because some correlations may be inadmissable
      #ps = polyserial(otu.t[m.ix,],map$Subject,ML=T,std.err = T)
      #if (is.na(pchisq(ps$chisq, ps$df))) next # Invalid correlation
      #Grp.Corrs[m.ix] = ps$rho             # Find intensity of correlation
      #Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
      # Wilcoxon Test
      one <- otu.t[m.ix,grep("a",map$Subject,value=F)] # Only subject id's #a
      two <- otu.t[m.ix,grep("b",map$Subject,value=F)] # Only subject id's #b
      Grp.Pvals[m.ix] <- wilcox.test(one, two)$p.value
      Grp.Corrs[m.ix] <- biserial.cor(otu.t[m.ix,], gsub('[ab]','',map$Subject), use = "complete.obs")
    },silent=T)
    #KW.Pvals[m.ix] = kruskal.test(otu.t[m.ix,] ~ map$Subject)$p.val
    TT.Pvals[m.ix] = t.test(otu.t[m.ix,], gsub('[ab]','',map$Subject), conf.level = 0.95)$p.value
  }
  
  ## Taxa barplots -- Top 15 most abundant (kruskal sig. + other?)
  otu.m = otu.n
  otu.m = sweep(sqrt(otu.m),2,colSums(sqrt(otu.m)),'/')
  meanAb = apply(otu.m,1,FUN=function(x) tapply(x, map$Subject, mean)) # group mean
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
  
  ## Plot according to Subject, sorted by abundance
  pdf(paste0("results/gg97_stomach_feces/TaxaSummary_stomachvsfeces_subject_gg97_L",bT[L],".pdf"),width = 8,height=7) # Make room for legend
  plot(ggplot(otu.m, aes(x = Subject, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Subject",y="Root Relative Abundance") +
         guides(fill=guide_legend(ncol=1)) +
         theme(panel.background = element_blank()) +
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
  res = res[order(res$KW.Pvals),]
  
  # Add bivariate filter
  sig = 0.05
  selection = res$TT.Pvals < sig
  
  # Display all significant with p < 0.05
  num_sig = sum(selection, na.rm = T) # Count how many are significant
  res = res[selection,]
  
  # Truncate names to last 2 informative levels
  split = strsplit(rownames(res),";")        # Split by semicolon into levels
  res$short = sapply(split,function(x) paste(tail(x,2),collapse=";"))
  
  pdf(paste0("results/gg97_stomach_feces/TaxaSwarms_stomachvsfeces_subject_gg97_L",bT[L],".pdf"),width = 6.5,height=6.5)
  sink(paste0("results/gg97_stomach_feces/Taxa_Significance_stomachvsfeces_subject_gg97_L",bT[L],".txt"))   # Get ready to write the significant ones
  cat("Taxon\tWilcoxon_P\tBiserial_Cor\tSubject_P\n")  # Print header
  if (num_sig) for (i in 1:num_sig) {
    taxon = rownames(res)[i]
    cat(res[taxon,]$short,'\t',res$Grp.Pvals[i],'\t',-res$Grp.Corrs[i],'\t',res$KW.Pvals[i],'\t','\n',sep='')
    beeswarm(otu.t[taxon,] ~ map$Subject, xlab="Subject",ylab="CLR Relative Abundance",main=res[taxon,]$short,
             col=alpha(lscolors,0.7),
             cex.axis=1.1,cex.main=1,cex=1.1,corral="random",pch=19)
    bxplot(otu.t[taxon,] ~ map$Subject, add = TRUE)
  }
  sink(NULL)
  dev.off()
  
  ## Heatmap ##
  if (num_sig < 2) next
  # Need to have created the clr taxa table as a matrix
  my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299) # Sebastian Raschka
  gl = map$Subject
  glpos = c(grep("1a",gl),grep("1b",gl),grep("2a",gl),grep("2b",gl),grep("3a",gl),grep("3b",gl),grep("4a",gl),grep("4b",gl),grep("5a",gl),grep("5b",gl),grep("6a",gl),grep("6b",gl))
  gl = gl[glpos]
  mat = otu.t[rownames(res),glpos]
  
  # Truncate names to last 2 informative levels
  split = strsplit(as.character(rownames(mat)),";")        # Split by semicolon into levels
  rownames(mat) = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  
  levels(gl)= c("red","green") # lscolors
  png(paste0("results/gg97_stomach_feces/Taxa_heatmap_stomachvsfeces_subject_gg97_L",bT[L],".png"),  # create PNG for the heat map        
      width = 8*300,                        # 8 x 300 pixels
      height = 6*300,
      res = 300,                              # 300 pixels per inch
      pointsize = 10)                          # smaller font size
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
library(polycor)
library(robCompositions)
library(beeswarm)

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
#KW.Pvals=rep(1,npaths)
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
  #KW.Pvals[m.ix] = kruskal.test(picrust[m.ix,] ~ map$Bodysite)$p.val
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
  beeswarm(picrust[pathway,] ~ map$Bodysite, xlab="Body Site",ylab="Pathway Abundance",main=pathway,
           col=alpha(lscolors,0.7),
           cex.axis=1.1,cex.main=1,cex=1.1,corral="random",pch=19)
  bxplot(picrust[pathway,] ~ map$Bodysite, add = TRUE)
}
sink(NULL)
dev.off()

# PICRUSt heatmap too, why not
library(gplots)
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299) # Sebastian Raschka
gl = map$Bodysite
glpos = c(grep("Foregut",gl),grep("Hindgut",gl))
gl = gl[glpos]
mat = picrust[rownames(res[abs(res$Grp.Corrs) > 0.75,]),glpos]
mat = sweep(mat,1,rowSums(abs(mat)),'/')                      # Normalize to relative abundance
mat = sweep(mat,1,max(apply(mat,1,max),apply(mat,1,min)),'/') # Constrain extrema to [-1, 1]

levels(gl)= lscolors #c("red","orange","yellow")

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

