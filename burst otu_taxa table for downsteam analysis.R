#### Loading ####
setwd("/data/analyses/");                         # Set file location
taxa = read.delim("taxa.txt",row=1);              # Read taxa table

#### Summarizing ####
L = 5;                                            # What level to summarize at (editable)
split = strsplit(rownames(taxa),";");             # Split and rejoin at desired level
taxaStrings = sapply(split,function(x) paste(x[1:L],collapse=";"));
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
taxa = rowsum(taxa,taxaStrings);                  # Collapse by taxonomy name
# Uncomment the next line to add NA's up to 'L' levels; replace L with desired level
#rownames(taxa) = sapply(strsplit(rownames(taxa),";"),function(x) paste(x[1:L],collapse=";"));

#### Massaging ####
taxa = sweep(taxa,2,colSums(taxa),'/');           # Normalize
taxa = taxa[order(rowMeans(taxa),decreasing=T),]; # Sort by avg. abundance (for looks)
taxa = taxa[rowMeans(taxa) >= 0.001,];            # Drop rare taxa (abundance)
taxa = taxa[rowSums(taxa > 0) > 3,];              # Drop rare taxa (prevalence)

#### Saving ####
sink("taxa_processed.txt"); cat("#Taxonomy\t");
write.table(taxa,file="taxa_processed.txt",quote=F,sep="\t",append = T);
sink(NULL)