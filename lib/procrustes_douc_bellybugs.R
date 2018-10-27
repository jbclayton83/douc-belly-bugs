# Generate a more logical plot for procrustes comparisons
### By Robin Shields-Cutler

library(ape)
library(vegan)
library(ggplot2)

# Metadata needs three columns - 
# 1. sample IDs (e.g. "P21_pre_stool")
# 2. the unifying participant/unit ID (e.g."P21") 
# 3. the binary metadata group, corresponding to the two distance input matrices (e.g. "pre" or "post")
meta <- read.delim('../data/gg97/procrustes_map.txt',
                   header=1, row.names=1, check.names = F)

thing <- 'AnimalID_name'  # The meta header for the participant ID column
group <- 'Bodysite'  # The meta header for the metadata group category (which of the two distance matrices)
groups <- c('Foregut','Hindgut')  # Names of the two categories in the group
metaA <- meta[meta[,group] == groups[1],]  # Split the metadata into the two groups
metaB <- meta[meta[,group] == groups[2],]

metaA <- metaA[order(metaA[,thing]),]  # Sort the metadata by unifying ID
metaB <- metaB[order(metaB[,thing]),]

# CRITICAL:
# Ensure that the original distance matrices are in the same order by participant
# Uses the order generated from the metadata dataframe
run_protest <- function(intable_A, intable_B, metaA=metaA, metaB=metaB) {
  betatable_A <- intable_A[rownames(metaA),rownames(metaA)]
  betatable_B <- intable_B[rownames(metaB),rownames(metaB)]
  
  # Get the principal coordinates
  pcoa_A <- pcoa(betatable_A)$vectors
  for(c in 1:ncol(pcoa_A)){
    colnames(pcoa_A)[c] <- paste0("PC",c)
  }
  pcoa_B <- pcoa(betatable_B)$vectors
  for(c in 1:ncol(pcoa_B)){
    colnames(pcoa_B)[c] <- paste0("PC",c)
  }
  
  # crusty <- procrustes(pcoa_A, pcoa_B, symmetric = T)  # Run Procrustes
  crust_test_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 999))$signif
  cat(crust_test_p)
}

# Read in beta diversity tables
bray_A <- read.delim('../data/gg97/beta_div_stomach/unifrac_douc_stomach_otutable_gg97.txt',
                     header=1, row.names = 1, check.names = F)
bray_B <- read.delim('../data/gg97/beta_div_feces/bray_curtis_douc_feces_otutable_gg97.txt',
                     header=1, row.names = 1, check.names = F)


unw_A <- read.delim('../data/gg97/beta_div_stomach/unifrac_douc_stomach_otutable_gg97.txt',
                    header=1, row.names = 1, check.names = F)
unw_B <- read.delim('../data/gg97/beta_div_feces/unifrac_douc_feces_otutable_gg97.txt',
                    header=1, row.names = 1, check.names = F)


wuni_A <- read.delim('../data/gg97/beta_div_stomach/weighted_unifrac_douc_stomach_otutable_gg97.txt',
                     header=1, row.names = 1, check.names = F)
wuni_B <- read.delim('../data/gg97/beta_div_feces/weighted_unifrac_douc_feces_otutable_gg97.txt',
                     header=1, row.names = 1, check.names = F)


run_protest(bray_A, bray_B, metaA, metaB)
run_protest(unw_A, unw_B, metaA, metaB)
run_protest(wuni_A, wuni_B, metaA, metaB)

pvals <- c(0.132,.076,0.861)
p.adjust(pvals, method = 'fdr')






betatable_A <- read.delim('../data/gg97/beta_div_stomach/unifrac_douc_stomach_otutable_gg97.txt',
                    header=1, row.names = 1, check.names = F)
betatable_B <- read.delim('../data/gg97/beta_div_feces/unifrac_douc_feces_otutable_gg97.txt',
                    header=1, row.names = 1, check.names = F)


# Metadata needs three columns - 
# 1. sample IDs (e.g. "P21_pre_stool")
# 2. the unifying participant/unit ID (e.g."P21") 
# 3. the binary metadata group, corresponding to the two distance input matrices (e.g. "pre" or "post")
meta <- read.delim('../data/gg97/procrustes_map.txt',
                   header=1, row.names=1, check.names = F)

thing <- 'AnimalID_name'  # The meta header for the participant ID column
group <- 'Bodysite'  # The meta header for the metadata group category (which of the two distance matrices)
groups <- c('Foregut','Hindgut')  # Names of the two categories in the group
metaA <- meta[meta[,group] == groups[1],]  # Split the metadata into the two groups
metaB <- meta[meta[,group] == groups[2],]

metaA <- metaA[order(metaA[,thing]),]  # Sort the metadata by unifying ID
metaB <- metaB[order(metaB[,thing]),]

# CRITICAL:
# Ensure that the original distance matrices are in the same order by participant
# Uses the order generated from the metadata dataframe
betatable_A <- betatable_A[rownames(metaA),rownames(metaA)]
betatable_B <- betatable_B[rownames(metaB),rownames(metaB)]

# Get the principal coordinates
pcoa_A <- pcoa(betatable_A)$vectors
for(c in 1:ncol(pcoa_A)){
  colnames(pcoa_A)[c] <- paste0("PC",c)
}
pcoa_B <- pcoa(betatable_B)$vectors
for(c in 1:ncol(pcoa_B)){
  colnames(pcoa_B)[c] <- paste0("PC",c)
}

crusty <- procrustes(pcoa_A, pcoa_B, symmetric = T)  # Run Procrustes
crust_test_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 999))$signif