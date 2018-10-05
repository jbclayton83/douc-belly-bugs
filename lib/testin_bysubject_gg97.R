#### Differential taxa testing - by Subject ####
library(gplots)
library(RColorBrewer)
library(robCompositions) # Composition magic
library(polycor)
library(beeswarm)
library(reshape2)
library(phyloseq)
library(gtools)

bT = c(2,4,6,7)  # Levels to consider for the bugs
lscolors = c("red","green")
L=3 # assume 6th level for now

# For Loop contents
# Load OTU data for calculating distances
otu.d <- as.matrix(otu) # get matrix
# Filter out bugs that don't appear in enough samples
select = rowSums(otu.d > 0) > min(table(map$Subject))   # Reasonable general min
otu.d = otu.d[select,]                                  # Apply drop mask
otu.d = otu.d[,rownames(map)]                           # Sync with map order
otu.dn = sweep(otu.d, 2, colSums(otu), '/');             # Normalize to relative abundance

# Convert to relative abundance - CLR
eps = 0.5
otu.d = otu.d * (1 - rowSums(otu.d==0) * eps / rowSums(otu.d))
otu.d[otu.d==0] <- eps
otu.d = sweep(otu.d,1,rowSums(otu.d),'/');
ls = log(otu.d) # 
otu.d = t(ls - rowMeans(ls))
otu.d <- otu.d[, !is.nan(colSums(otu.d))]

# Calculate Foregut/Hindgut Relationship per Subject
phyobj <- phyloseq(otu_table(otu.d, taxa_are_rows=F), tree)
obs <- UniFrac(phyobj, weighted=F) #observed unweighted UniFrac distances --note:Is computing in parallel possible on this platform?
# permuation to determine validity of distance
nperm = 1000 # Number of permuations desired - technically should be done num.samples*num.otus times, but that is too large
samples = rownames(otu.d)
matnames = rep("permuation",1000)
permvals <- array(, dim=c(12,12,nperm),dimnames=list(samples,samples,matnames))
otu.temp <- otu.d
for (i in 1:nperm) {
  #randomly assign community labels to OTUs
  otu.temp <- sample(0:1, size=12*ncol(otu.d) , replace=T)
  #create phyloseq object
  phyobj.temp <- phyloseq(otu_table(otu.temp, taxa_are_rows=F), tree)
  #measure Unifrac distances & update permvals
  permvals[,,i] <- as.matrix(UniFrac(phyobj, weighted=F))
}
# get p-values for distances
obs.ff <- obs[]# foreguts across all subjects
obs.hh <- obs[]# hindguts across all subjects
obs.fh <- obs[]# foregut-hindgut for each subject






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
otu.t = as.matrix(otu.t)

# Go through each taxon and test for significance w/group
ntax = nrow(otu.t)
BS = map$Subject
Grp.Pvals=rep(1,ntax)
Grp.Corrs=rep(0,ntax)
KW.Pvals=rep(1,ntax)
#TT.Pvals=rep(1,ntax)
for (m.ix in 1:ntax) {  # Loop through all the rows (taxa)
  try({ # Because some correlations may be inadmissable
    ps = polyserial(otu.t[m.ix,],map$Subject,ML=T,std.err = T)
    if (is.na(pchisq(ps$chisq, ps$df))) next # Invalid correlation
    Grp.Corrs[m.ix] = ps$rho             # Find intensity of correlation
    Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
    # Wilcoxon Test
    #one <- otu.t[m.ix,grep("a",map$Subject,value=F)] # Only subject id's #a
    #two <- otu.t[m.ix,grep("b",map$Subject,value=F)] # Only subject id's #b
    #Grp.Pvals[m.ix] <- wilcox.test(one, two)$p.value
    #Grp.Corrs[m.ix] <- biserial.cor(otu.t[m.ix,], gsub('[ab]','',map$Subject), use = "complete.obs")
  },silent=T)
  KW.Pvals[m.ix] = kruskal.test(otu.t[m.ix,] ~ map$Subject)$p.val
  #TT.Pvals[m.ix] = t.test(otu.t[m.ix,], gsub('[ab]','',map$Subject), conf.level = 0.95)$p.value
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
otu.m$Subject = gsub("a", "f", otu.m$Subject) # change name to #f for foregut
otu.m$Subject = gsub("b", "h", otu.m$Subject) # change name to #h for hindgut

# Final Taxa Summary Plots would go here, they are already working. so I did not include them.