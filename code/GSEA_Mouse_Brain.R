#####Gene set enrichment
cytoscapePing()
BiocManager::install("RCy3")
library(RCy3)
installApp('WikiPathways') 
installApp('CyTargetLinker') 
installApp('stringApp') 
library(clusterProfiler)

### GO
topgenes_gini = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']

sc_norm_exp <- 2^(sc_matrix2norm@norm_expr)-1

ExprSubset<-sc_norm_exp[as.character(topgenes_gini$genes),]

length(unique(rownames(ExprSubset)))

Sig2 <- Sig
colnames(Sig2) <- c("Oligos_5", "Neurons_48", "Vascular_67", "Astrocytes_41", "Oligos_0", "Neurons_27", 
                    "Oligos_1", "Oligos_14", "Ependymal_47", "Neurons_60", "Oligos_53", "Neurons_26", "Vascular_68", "Neurons_59",
                    "Astrocytes_14", "Neurons_49", "Neurons_12", "Neurons_58", "Vascular_70", "Astrocytes_40", "Immune_35", 
                    "Astrocytes_42", "Neurons_28", "Neurons_25", "Neurons_21", "Neurons_52", "Neurons_62", "Neurons_51", "Vascular_14", 
                    "Neurons_61", "Excluded_38", "Astrocytes_44", "Neurons_30", "Excluded_6", "Immune_32", "Immune_34", "Neurons_14", 
                    "Neurons_18", "Immune_14", 
                    "Neurons_15", "Neurons_11", "Neurons_63", "Blood_73", "Neurons_24", "Excluded_44", "Vascular_69", "Excluded_30", 
                    "Neurons_20", "Neurons_54", 
                    "Neurons_23", "Vascular_46", "Neurons_17", "Neurons_10", "Neurons_16", "Neurons_22", "Neurons_19")

Astrocytes_14 <- gini_markers_subclusters[(gini_markers_subclusters$cluster=="Astrocytes_14"),]
Blood_73 <- gini_markers_subclusters[(gini_markers_subclusters$cluster=="Blood_73"),]
Ependymal_47<- gini_markers_subclusters[(gini_markers_subclusters$cluster=="Ependymal_47"),]
Neurons_25 <- gini_markers_subclusters[(gini_markers_subclusters$cluster=="Neurons_25"),]
Oligos_5 <- gini_markers_subclusters[(gini_markers_subclusters$cluster=="Oligos_5"),]
Vascular_70 <- gini_markers_subclusters[(gini_markers_subclusters$cluster=="Vascular_70"),]

bkgd.genes_mouse <- gini_markers_subclusters$genes
Astrocytes_14_GO_up <- Astrocytes_14[Astrocytes_14$expression > 1.5 & Astrocytes_14$comb_score >0.05, 1] %>% drop_na()
Blood_73_GO_up <- Blood_73[Blood_73$expression > 1.5 & Blood_73$comb_score>0.05, 1] %>% drop_na()
Ependymal_47_GO_up <- Ependymal_47[Ependymal_47$expression > 1.5 & Ependymal_47$comb_score > 0.05, 1] %>% drop_na()
Neurons_25_GO_up <- Neurons_25[Neurons_25$expression > 1.5 & Neurons_25$comb_score > 0.05, 1] %>% drop_na()
Oligos_5_GO_up <- Oligos_5[Oligos_5$expression > 1.5 & Oligos_5$comb_score > 0.05, 1] %>% drop_na()
Vascular_70_GO_up <- Vascular_70[Vascular_70$expression > 1.5 & Vascular_70$comb_score > 0.05, 1] %>% drop_na()
view(Neurons_25)

up.genes.entrez_astro14_mouse <- clusterProfiler::bitr(Astrocytes_14_GO_up$genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

up.genes.entrez_Blood_73_mouse <- clusterProfiler::bitr(Blood_73_GO_up$genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

up.genes.entrez_Ependymal_47_mouse <- clusterProfiler::bitr(Ependymal_47_GO_up$genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

up.genes.entrez_Neurons_25_mouse <- clusterProfiler::bitr(Neurons_25_GO_up$genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

up.genes.entrez_Oligos_5_mouse <- clusterProfiler::bitr(Oligos_5_GO_up$genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

up.genes.entrez_Vascular_70_mouse <- clusterProfiler::bitr(Vascular_70_GO_up$genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

bkgd.genes_entrez_mouse <- clusterProfiler::bitr(bkgd.genes_mouse, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)


h_gene_sets = msigdbr::msigdbr(species = "mouse", category = "C8")


hs_gene_sets = msigdbr::msigdbr(species = "human", category = "C8")

h_gene_sets$gs_name == hs_gene_sets$gs_name


m_t2g <- msigdbr(species = "Mus musculus", category = "C8") %>% 
  dplyr::select(gs_name, entrez_gene) #this is not mouse, this calls human gene sets

install.packages("rjson")
install.packages("tidyjson")

library(rjson)
library(tidyjson)
myData <- fromJSON(file="../Downloads/download_file.jsp")

#https://bioconductor.org/packages/release/data/experiment/vignettes/msigdb/inst/doc/msigdb.html
BiocManager::install("msigdb")
library(msigdb)
library(ExperimentHub)
library(GSEABase)

package.version("msigdb")

eh = ExperimentHub()
query(eh , 'msigdb')

mm <- eh[['EH6779']]

msigdb.mm = getMsigdb(org = 'mm', id = 'EZID', version = '7.4')

msigdb.mm = appendKEGG(msigdb.mm)

listCollections(msigdb.mm)

table(sapply(lapply(msigdb.mm, collectionType), bcCategory))

listSubCollections(msigdb.mm)

mm_t2g<- subsetCollection(msigdb.mm, 'c8') #call mouse msigdb

m_msigdb_ids = geneIds(mm_t2g)


dfs <- lapply(m_msigdb_ids, data.frame, stringsAsFactors = TRUE)

dfs_rbind <- as.data.frame(do.call(rbind,            # Convert nested list to data frame by column
                                   dfs))

aa <- dfs_rbind #rename it just in case

aa$gs_name <- rownames(aa) #add a new column
head(aa)

library(dplyr)
aa2<- aa %>%
  mutate(gs_name = gsub("(.*?)\\.\\d+", "\\1", gs_name)) #eliminate the extra numbers at the end of the column

aa2 <- aa2[,c(2,1)] #change the column ordering

colnames(aa2) <- c("gs_name", "entrez_gene") #relabel the columns

aa2 <- read.csv("../Desktop/Nature_Paper/mouse_oct2023_v7.5.1_c8_geneLists.csv")

aa3 <- aa2[,-1]

require(RColorBrewer)
brewer.pal(5, "Dark2")

library(clusterProfiler)

astro <- enricher(up.genes.entrez_astro14_mouse$ENTREZID, TERM2GENE = aa3)
blood <- enricher(up.genes.entrez_Blood_73_mouse$ENTREZID, TERM2GENE = aa3)
ependymal <- enricher(up.genes.entrez_Ependymal_47_mouse$ENTREZID, TERM2GENE = aa3)
neuron <- enricher(up.genes.entrez_Neurons_25_mouse$ENTREZID, TERM2GENE = aa3)
oligo <- enricher(up.genes.entrez_Oligos_5_mouse$ENTREZID, TERM2GENE = aa3)
vascular <- enricher(up.genes.entrez_Vascular_70_mouse$ENTREZID, TERM2GENE = aa3)


write_csv(as.data.frame(up.genes.entrez_astro14_mouse$ENTREZID), "../Desktop/Nature_Paper/astro_upgenes.csv")
write_csv(as.data.frame(up.genes.entrez_Blood_73_mouse$ENTREZID), "../Desktop/Nature_Paper/blood_upgenes.csv")
write_csv(as.data.frame(up.genes.entrez_Ependymal_47_mouse$ENTREZID), "../Desktop/Nature_Paper/ependymal_upgenes.csv")
write_csv(as.data.frame(up.genes.entrez_Neurons_25_mouse$ENTREZID), "../Desktop/Nature_Paper/neuron_upgenes.csv")
write_csv(as.data.frame(up.genes.entrez_Oligos_5_mouse$ENTREZID), "../Desktop/Nature_Paper/oligo_upgenes.csv")
write_csv(as.data.frame(up.genes.entrez_Vascular_70_mouse$ENTREZID), "../Desktop/Nature_Paper/vascular_upgenes.csv")

colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854") #colorblind friendly colors

par(mfrow = c(1, 1))

dotplot(astro, showCategory = 15)

dotplot(astro, showCategory = 15) + scale_colour_gradient2(low = "#E41A1C", mid = "#D95F02", high = "#1B9E77")
dotplot(blood, showCategory = 15) + scale_colour_gradient2(low = "#E41A1C", mid = "#D95F02", high = "#1B9E77")
dotplot(ependymal, showCategory = 15) + scale_colour_gradient2(low = "#E41A1C", mid = "#D95F02", high = "#1B9E77")
dotplot(neuron, showCategory = 15) + scale_colour_gradient2(low = "#E41A1C",  mid = "#D95F02", high = "#1B9E77")
dotplot(oligo, showCategory = 15) + scale_colour_gradient2(low = "#E41A1C",  mid = "#D95F02", high = "#1B9E77")
dotplot(vascular, showCategory = 15) + scale_colour_gradient2(low = "#E41A1C",mid = "#D95F02", high = "#1B9E77")

###end GO

library(org.Mm.eg.db)
##up and downregulated genes
up.genes_ST_mouse <- counts.metaST_mouse_allclust[counts.metaST_mouse_allclust$avg_log2FC > 0 & counts.metaST_mouse_allclust$p_val_adj < 0.01, 7] 
dn.genes_ST_mouse <- counts.metaST_mouse_allclust[counts.metaST_mouse_allclust$avg_log2FC < 0 & counts.metaST_mouse_allclust$p_val_adj < 0.01, 7]
bkgd.genes_ST_mouse <- counts.metaST_mouse_allclust[,7]
#Enrichment

BiocManager::install("org.Mm.eg.db")
library("org.Hs.eg.db")
BiocManager::install("rWikiPathways")
library(magrittr)
library("rWikiPathways")
load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

write.csv(up.genes_ST_mouse, file = "up.genes_ST_mouse.csv") 
write.csv(dn.genes_ST_mouse, file = "dn.genes_ST_mouse.csv") 
write.csv(bkgd.genes_ST_mouse, file = "bkgd.genes_ST_mouse.csv") 

up.genes_ST_mouse_2 <- read.csv("up.genes_ST_mouse_3.csv",header=FALSE)
up.genes_ST_mouse_3 <- unlist(up.genes_ST_mouse_2, recursive = TRUE, use.names = TRUE)

dn.genes_ST_mouse_2 <- read.csv("dn.genes_ST_mouse_3.csv",header=FALSE)
dn.genes_ST_mouse_3 <- unlist(dn.genes_ST_mouse_2, recursive = TRUE, use.names = TRUE)

bkgd.genes_ST_mouse_2 <- read.csv("bkgd.genes_ST_mouse_3.csv",header=FALSE)
bkgd.genes_ST_mouse_3 <- unlist(bkgd.genes_ST_mouse_2, recursive = TRUE, use.names = TRUE)

up.genes.entrez_ST_mouse <- clusterProfiler::bitr(up.genes_ST_mouse_3, fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_ST_mouse <- clusterProfiler::bitr(dn.genes_ST_mouse_3, fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_ST_mouse <- clusterProfiler::bitr(bkgd.genes_ST_mouse_3, fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


##Gene ontology
egobp_ST_mouse <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_ST_mouse[[2]],
  universe = bkgd.genes.entrez_ST_mouse[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

head(egobp_ST_mouse,10)

barplot(egobp_ST_mouse, showCategory = 20)
dotplot(egobp_ST_mouse, showCategory = 20)
emapplot(egobp_ST_mouse, showCategory = 20)
goplot(egobp_ST_mouse)
ggplot(egobp_ST_mouse[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

##Downgenes
##Gene ontology
egobp_ST_mouse_dn <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_ST_mouse[[2]],
  universe = bkgd.genes.entrez_ST_mouse[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

head(egobp_ST_mouse_dn,10)

barplot(egobp_ST_mouse_dn, showCategory = 20)
dotplot(egobp_ST_mouse_dn, showCategory = 20)
emapplot(egobp_ST_mouse_dn, showCategory = 20)

goplot(egobp_ST_mouse_dn)

ggplot(egobp_ST_mouse_dn[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))


##KEGG
kegg_ST_mouse <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_ST_mouse[[2]],
  universe = bkgd.genes.entrez_ST_mouse[[2]],
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  use_internal_data = FALSE)
head(kegg_ST_mouse)

barplot(kegg_ST_mouse, showCategory = 10)
dotplot(kegg_ST_mouse, showCategory = 10)
emapplot(kegg_ST_mouse, showCategory = 10)
goplot(kegg_ST_mouse)

#Wikipathways
wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")

# supported organisms (and case-dependent spelling)...
listOrganisms()

wp2gene <- readPathwayGMT(wp.hs.gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

cat("\n\nWhich column contains my new Entrez IDs?\n")
head(up.genes.entrez_ST_mouse)
keytypes(org.Hs.eg.db)
ewp.up <- clusterProfiler::enricher(
  up.genes.entrez_ST_mouse[[2]],
  universe = bkgd.genes.entrez_ST_mouse[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)


ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)
ewp.up_sim <- pairwise_termsim(ewp.up)

barplot(ewp.up_sim, showCategory = 20)
dotplot(ewp.up_sim, showCategory = 20)
emapplot(ewp.up_sim, showCategory = 20)
goplot(ewp.up_sim)

## convert gene ID to Symbol
ewp.upx <- setReadable(ewp.up, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(ewp.upx)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(ewp.upx, categorySize="pvalue")
p3 <- cnetplot(ewp.upx, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
cowplot::plot_grid(p3, ncol=1, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

p4 <- heatplot(ewp.upx)
p5 <- emapplot(ewp.up)
cowplot::plot_grid(p4, ncol=1, labels=LETTERS[1:2])


ewp.dn <- clusterProfiler::enricher(
  dn.genes.entrez_ST_mouse[[2]],
  universe = bkgd.genes.entrez_ST_mouse[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)


ewp.dn <- DOSE::setReadable(ewp.dn, org.Hs.eg.db, keyType = "ENTREZID")

dotplot(ewp.dn, showCategory = 20)

counts.metaST_mouse_allclust$fcsign <- sign(counts.metaST_mouse_allclust$avg_log2FC)
counts.metaST_mouse_allclust$logfdr <- -log10(counts.metaST_mouse_allclust$p_val + 1)
counts.metaST_mouse_allclust$sig <- counts.metaST_mouse_allclust$logfdr/counts.metaST_mouse_allclust$fcsign
sig.counts.metaST_mouse_allclust.entrez<-merge(counts.metaST_mouse_allclust, bkgd.genes.entrez_sc, by.x = "gene", by.y = "SYMBOL")
gsea.sig.counts.metaST_mouse_allclust <- sig.counts.metaST_mouse_allclust.entrez[,10]
names(gsea.sig.counts.metaST_mouse_allclust) <- as.character(sig.counts.metaST_mouse_allclust.entrez[,11])
gsea.sig.counts.metaST_mouse_allclust <- sort(gsea.sig.counts.metaST_mouse_allclust,decreasing = TRUE)
head(gsea.sig.counts.metaST_mouse_allclust)

gwp.counts.metaST_mouse_allclust <- clusterProfiler::GSEA(
  gsea.sig.counts.metaST_mouse_allclust,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.73, #p.adjust cutoff
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

gwp.sig.counts.metaST_mouse_allclust.df = data.frame(ID=gwp.counts.metaST_mouse_allclust$ID,
                                                     Description=gwp.counts.metaST_mouse_allclust$Description,
                                                     enrichmentScore=gwp.counts.metaST_mouse_allclust$enrichmentScore,
                                                     NES=gwp.counts.metaST_mouse_allclust$NES,
                                                     pvalue=gwp.counts.metaST_mouse_allclust$pvalue,
                                                     p.adjust=gwp.counts.metaST_mouse_allclust$p.adjust,
                                                     rank=gwp.counts.metaST_mouse_allclust$rank,
                                                     leading_edge=gwp.counts.metaST_mouse_allclust$leading_edge
)


gwp.sig.counts.metaST_mouse_allclust.df[which(gwp.sig.counts.metaST_mouse_allclust.df$NES > 1),] #pathways enriched for upregulated cardiac genes
gwp.sig.counts.metaST_mouse_allclust.df[which(gwp.sig.counts.metaST_mouse_allclust.df$NES < -1),] #pathways enriched for downregulated cardiac genes
