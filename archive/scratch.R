### Scratch for new WHAM! revisions 3.28
#library(R.utils)
#library(SummarizedExperiment)
#sourceDirectory("ALDEx2", verbose=T)

# source("ALDEx2/AllGenerics.R")
# source("ALDEx2/AllClasses.R")
# source("ALDEx2/progress.R")
# source("ALDEx2/clr_effect.r")
# source("ALDEx2/clr_function.r")
# source("ALDEx2/clr_ttest.r")
# source("ALDEX2/aldex.r")
# source("ALDEX2/stats.fast.R")
# source("ALDEX2/rdirichlet.R")
# source("ALDEX2/iqlr_features.r")
# source("ALDEX2/clr_glm.r")
# source("ALDEX2/clr_corr.r")
# source("ALDEX2/plot.aldex.r")

#####

setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28/")

library(ggplot2)
### figure out how to  get taxa levels from EBI with semicolons

buggah = fread("ERP103961_taxonomy_abundances_LSU_v4.0.tsv", header=TRUE, sep='\t')

bugg_fix <- gsub(";k__", "_k__", buggah$Taxa)
bug_levels <- strsplit(bugg_fix, ";")

tax_lev <- c()
for (i in 1:length(bug_levels)){
  bug_len <- length(bug_levels[[i]])
  if (bug_len == 4){
    tax_lev <- c(tax_lev, i)
  }
}

buggah_new <- buggah[tax_lev,]

specs <- buggah_new

### why wont it work ###
spec_col = (ncol(specs))-1
spec_row = nrow(specs)

spec_RelExp = data.frame(specs[,-1])
spec_RelExp2 = t(spec_RelExp)
df_spec_RelExp = data.frame(spec_RelExp2)

spec_datas = c()
for (i in 1:spec_row) {
  spec_data=as.vector(t(df_spec_RelExp[i]))
  spec_datas = c(spec_datas,spec_data)
}

group_titles2 <- rep(c(rep("g1", 300), rep("g2", 416)), 16)

spec_sample_num = paste(1:(spec_col))
spec_samp=strtoi(spec_sample_num)
spec_isamp = rep(spec_samp,spec_row)

spec_spec <- rep(specs$Taxa, each = spec_col)

spec_g_data = data.frame(spec_datas,
                         spec_isamp, group_titles2, spec_spec)
colnames(spec_g_data) = c("Relative_Abundance", "Sample_num", "Group", "Taxa")

uniq_spec_num = length(unique(spec_g_data$Taxa))

spec_g_data$Taxa <- reorder(spec_g_data$Taxa, -spec_g_data$Relative_Abundance)

print(spec_g_data)

spec_g_data2 <- aggregate(spec_g_data$Relative_Abundance, list(spec_g_data$Taxa), sum)

summer <- sum(spec_g_data2[,2])
spec_g_data2$prop <- spec_g_data2[,2]/summer
spec_g_data3 <- subset(spec_g_data2, prop < 1)#input$upper_limit)
spec_g_data4 <- subset(spec_g_data3, prop > 0)#input$lower_limit)

cc = quantile(spec_g_data2$prop, 0.5)

keep_taxa <- as.character(unlist(spec_g_data4[,1]))
spec_g_data_filt <- subset(spec_g_data, Taxa %in% keep_taxa)

exp_plot = ggplot(spec_g_data_filt, aes(x = Group,
                                        y = Relative_Abundance,
                                        fill = Taxa)) +
  stat_summary(fun.y = "mean", geom = "bar", position = "fill")+
  scale_fill_manual(values = randomColor(uniq_spec_num)) +
  ylab("Relative Abundance") + theme(legend.position = 'none')

pp = ggplotly(exp_plot)


## biobakery
buggah = fread("gene_families_countEst_go_names_small.tsv", header=TRUE, sep='\t')


bug_levels <- strsplit(buggah$Taxa, "\\.")
head(bug_levels)

bug_keep <- gsub("\\..*", "", buggah$Taxa)
tail(bug_keep)
for (i in 1:length(bug_levels)){
  bug_len <- bug_levels[[i]]
  bug_keep <- c(bug_keep, bug_len[1])
}
buggah$Taxa <- bug_keep

buggah_new <- buggah[tax_lev,]





#
#
#
#

#
#
#
#### differential abundance play

####### raw count estimation
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.7")
library(data.table)
#library(metagenomeSeq)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)

go_relab <- fread("gene_families_normRELAB_go_names_small.tsv",
                  header=T, sep='\t')

go_df <- fread("gene_families_countEst_go_names_small.tsv", header=T, sep='\t')

go_est <- data.frame(Gene_Family = go_df$Gene_Family, go_df[,4:ncol(go_df)])

go_est_slim <- aggregate(go_est[,-1], list(go_est$Gene_Family), sum)

go_counts <- go_est_slim[,-1]
rownames(go_counts) <- go_est_slim[,1]

####get the right groups

arm = c("SRR532024_Abundance-RPKs", "SRR532015_Abundance-RPKs", "SRR532006_Abundance-RPKs",
        "SRR532040_Abundance-RPKs", "SRR532504_Abundance-RPKs", "SRR532507_Abundance-RPKs",
        "SRR638753_Abundance-RPKs", "SRR640357_Abundance-RPKs", "SRR640340_Abundance-RPKs",
        "SRR640452_Abundance-RPKs", "SRR640499_Abundance-RPKs", "SRR545546_Abundance-RPKs")
vag = c("SRR062353_Abundance-RPKs", "SRR062357_Abundance-RPKs", "SRR062301_Abundance-RPKs",
        "SRR062276_Abundance-RPKs", "SRR1804686_Abundance-RPKs", "SRR1804628_Abundance-RPKs",
        "SRR514191_Abundance-RPKs", "SRR514180_Abundance-RPKs", "SRR513168_Abundance-RPKs",
        "SRR514231_Abundance-RPKs", "SRR513448_Abundance-RPKs")
sal = c("SRR062435_Abundance-RPKs", "SRR062441_Abundance-RPKs", "SRR062389_Abundance-RPKs",
        "SRR062413_Abundance-RPKs", "SRR062402_Abundance-RPKs", "SRR062396_Abundance-RPKs",
        "SRR346673_Abundance-RPKs", "SRR346681_Abundance-RPKs", "SRR062371_Abundance-RPKs",
        "SRR062372_Abundance-RPKs", "SRR062462_Abundance-RPKs", "SRR062415_Abundance-RPKs")
sto = c("SRR528423_Abundance-RPKs", "SRR528353_Abundance-RPKs", "SRR528300_Abundance-RPKs",
        "SRR528261_Abundance-RPKs", "SRR528183_Abundance-RPKs", "SRR528155_Abundance-RPKs",
        "SRR532178_Abundance-RPKs", "SRR532183_Abundance-RPKs", "SRR532190_Abundance-RPKs",
        "SRR532191_Abundance-RPKs", "SRR533152_Abundance-RPKs", "SRR533153_Abundance-RPKs")

orderer <- c(arm, vag, sal, sto)

sitesArm <- c(rep(1, 12), rep(0, 35))
sitesVag <- c(rep(0, 12), rep(1, 11), rep(0, 24))
sitesSal <- c(rep(0, 23), rep(1, 12), rep(0, 12))
sitesSto <- c(rep(0, 35), rep(1, 12))

hmp_meta <- data.frame(sample_id = orderer, 
                       sample = 1:47, 
                       BodySiteArm = sitesArm,BodySiteSal = sitesSal,
                       BodySiteSto = sitesSto,BodySiteVag = sitesVag)
rownames(hmp_meta) <- hmp_meta$sample_id
hmp_meta <- hmp_meta[colnames(go_counts),]
hmp_ADF <- AnnotatedDataFrame(data = hmp_meta)


hmp_wham <- newMRexperiment(counts = go_counts, phenoData = hmp_ADF)

normed_hmp <-  cumNorm(hmp_wham, p = 0.75)
hmp_sample_data <-  pData(normed_hmp)


full_results <- data.frame(logFC=NA, se = NA, pvalues = NA, adjPvalues = NA, gene = NA)
bodysites <- c("BodySiteArm", "BodySiteVag", "BodySiteSal", "BodySiteSto")
for (i in bodysites){
  this_one = as.formula(paste0("~", i))
  mod <-  model.matrix(this_one, data = hmp_sample_data)
  
  results_hmp <-  fitFeatureModel(normed_hmp, mod)
  logFC_hmp <- MRcoefs(results_hmp, number = nrow(normed_hmp))
  logFC_hmp$gene <- rownames(logFC_hmp)
  full_results <- rbind(full_results, logFC_hmp)
}

go_de <- subset(full_results, logFC > 5 | logFC < -5)
go_de_sig <- subset(go_de, adjPvalues < 0.05)
keepers <- unique(go_de_sig$gene)

##add taxa for now
#go_est$bug <- gsub(".*\\|", "", go_df$`Gene Family`)

go_show <- subset(go_est, `Gene Family` %in% keepers)

go_show_slim <- aggregate(go_show[,-1], list(go_show$`Gene Family`), sum)

go_show_nums <- go_show_slim[,-1]
rownames(go_show_nums) = go_show_slim[,1]

go_show_nums2 <- log10(go_show_nums + 1)

pdf("hmp_heat.pdf", height = 12, width = 10)
heatty = heatmap.3(go_show_nums2, scale = 'row', col = viridis(100),
                   margins = c(12,20))
dev.off()



##### differential taxa play
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)
#library(metagenomeSeq)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)

ebi_df <- fread("ERP103961_taxonomy_abundances_LSU_v4.0.tsv", header=T, sep='\t')

ebi_counts <- ebi_df[,-1]
rownames(ebi_counts) <- ebi_df$Taxa



####get the right groups

### this might be the hard part... need spec mat reorder
group1 <- colnames(ebi_counts)[1:50]
group2 <- colnames(ebi_counts)[101:150]

orderer <- c(group1, group2)

sitesG1 <- c(rep(1, 50), rep(0, 50))
sitesG2 <- c(rep(0, 50), rep(1, 50))

groups <- c("group1", "group2", "group3")
group_nums <- c(20, 40, 40)
meta <- data.frame(sample_id = colnames(ebi_counts)[1:100],
                   #sample = 1:ncol(ebi_counts))
                   sample = 1:100)

bin <- c()
for (i in 1:length(groups)){
  use <- rep(groups[i], group_nums[i])
  bin <- c(bin, use)
}
meta$bins <- bin

for (i in 1:length(groups)){
  new_bin <- as.integer(meta$bins == groups[i])
  meta <- cbind(meta, new_bin)
}

colnames(meta) <- c("sample_id", "sample", "bins", groups)
rownames(meta) <- meta$sample_id


ebi_counts <- data.frame(ebi_counts)
keep_samps <- as.character(unlist(meta$sample_id))

ebi_counts2 <- ebi_counts[,keep_samps]
rownames(ebi_counts2) <- ebi_df$Taxa

meta <- meta[colnames(ebi_counts2),]
hmp_ADF <- AnnotatedDataFrame(data = meta)


hmp_wham <- newMRexperiment(counts = ebi_counts2, phenoData = hmp_ADF)

normed_hmp <-  cumNorm(hmp_wham, p = 0.75)
hmp_sample_data <-  pData(normed_hmp)


full_results <- data.frame(logFC=NA, se = NA, pvalues = NA, adjPvalues = NA, gene = NA)
bodysites <- c("group1", "group2", "group3")
for (i in bodysites){
  this_one = as.formula(paste0("~", i))
  mod <-  model.matrix(this_one, data = hmp_sample_data)
  
  results_hmp <-  fitFeatureModel(normed_hmp, mod)
  logFC_hmp <- MRcoefs(results_hmp, number = nrow(normed_hmp))
  logFC_hmp$gene <- rownames(logFC_hmp)
  full_results <- rbind(full_results, logFC_hmp)
}

go_de <- subset(full_results, logFC > 0.5 | logFC < -0.5)
go_de_sig <- subset(go_de, adjPvalues < 0.6)
keepers <- unique(go_de_sig$gene)

##add taxa for now
#go_est$bug <- gsub(".*\\|", "", go_df$`Gene Family`)

go_show <- subset(ebi_df[,1:101], Taxa %in% keepers)

go_show_slim <- aggregate(go_show[,-1], list(go_show$Taxa), sum)

go_show_nums <- go_show_slim[,-1]
rownames(go_show_nums) = go_show_slim[,1]

go_show_nums2 <- log10(go_show_nums + 1)

colorss <- meta$bins

cc = length(unique(meta$bins))
cols <- randomColor(3)

for (i in 1:cc){
  colorss = gsub(groups[i], cols[i], colorss)
}

col_col <- t(matrix(colorss))
col_col <- rbind(col_col, col_col)

#pdf("hmp_heat.pdf", height = 12, width = 10)
heatty = heatmap.3(t(t(go_show_nums2)), scale = 'row', col = viridis(100),
          margins = c(5,25), Colv = F, labCol = "", 
          ColSideColors = t(col_col), side.height.fraction = 0.5)
#dev.off()


x <- sweep(go_show_nums2, 1L, rowMeans(go_show_nums2, na.rm = T), check.margin = FALSE)
sx <- apply(x, 1L, sd, na.rm = T)
x <- sweep(x, 1L, sx, "/", check.margin = FALSE)

Z_scored_CPM <- unlist(x)

names <- rep(colnames(go_show_nums2), nrow(go_show_nums2))
taxa <- rep(rownames(go_show_nums2), each = ncol(go_show_nums2))
groups <- rep(bin, nrow(go_show_nums2))

plotter <- data.frame(Z_scored_CPM, names, taxa, groups)
###
##

gnums <- paste0("samp", 1:100)
gnames <- bin
expr <- rep(1,100)
gplot <- data.frame(gnums, gnames, expr)
gplot$gnums <- factor(gplot$gnums, levels = as.character(unlist(gplot$gnums)))

gEmpty <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_blank() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank())

gTop <- ggplot(gplot, aes(gnums, expr, fill = gnames)) +
  geom_tile()+
  scale_fill_manual(values = c("red", "purple", "blue")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_blank())


gMain <- ggplot(plotter, aes(x=names, y = taxa, fill = Z_scored_CPM)) +
  geom_tile(aes(fill = Z_scored_CPM)) +
  ylab("") +
  scale_fill_viridis(discrete=F) +
  #scale_fill_gradient2(low="purple", mid = "green4", high = "yellow", midpoint = 1.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
library(gridExtra)
grid.arrange(gTop, gMain,
             ncol = 1, nrow = 2, widths = c(1), heights = c(1, 10))



##

taxa_df <- fread("ERP103961_taxonomy_abundances_LSU_v4.0.tsv", header=T, sep='\t')
feature_df <- fread("ERP103961_small.tsv", header=T, sep='\t')

tax <- taxa_df[,-1]

feat <- feature_df[,-1]
feat <- feat[,-1]

length(intersect(colnames(feat), colnames(tax)))
setdiff(colnames(feat), colnames(tax))




#
#
#
#
#
#
#### diff abundance by ANOVA with post hoc

setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.7")
library(data.table)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)

go_relab <- fread("gene_families_normRELAB_go_names_small.tsv",
                  header=T, sep='\t')

go_nums <- go_relab[,2:48]

go_relab$Gene_Family <- gsub("\\|.*", "", go_relab$`Gene Family`)
go_relab$Gene_Family <- gsub(".*:", "", go_relab$Gene_Family)
go_relab$Gene_Family <- gsub(" \\[", "[", go_relab$Gene_Family)

go_relab$Acc <- gsub(".*GO:", "", go_relab$`Gene Family`)
go_relab$Acc <- gsub(":.*", "", go_relab$Acc)
go_relab$Acc <- gsub("\\|.*", "", go_relab$Acc)

go_relab$Species <- gsub(".*\\|", "", go_relab$`Gene Family`)

new_go_df <- data.frame(Acc = go_relab$Acc, Gene_Family = go_relab$Gene_Family,
                        Species = go_relab$Species, go_nums)
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
#write.table(new_go_df, "Relab_go_names_wham.tsv", quote = F, sep ='\t', row.names = F)

###

### diff abundance by ANOVA
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)


go_df <- fread("Relab_go_names_wham.tsv", header=T, sep='\t')

go_nums <- go_df[,4:50]

## spec da

spec_df <- aggregate(go_nums, list(go_df$Taxa), sum)

spec_nums <- spec_df[,2:48]
rownames(spec_nums) <- spec_df[,1]

arm = c("SRR532024_Abundance.RPKs", "SRR532015_Abundance.RPKs", "SRR532006_Abundance.RPKs",
        "SRR532040_Abundance.RPKs", "SRR532504_Abundance.RPKs", "SRR532507_Abundance.RPKs",
        "SRR638753_Abundance.RPKs", "SRR640357_Abundance.RPKs", "SRR640340_Abundance.RPKs",
        "SRR640452_Abundance.RPKs", "SRR640499_Abundance.RPKs", "SRR545546_Abundance.RPKs")
vag = c("SRR062353_Abundance.RPKs", "SRR062357_Abundance.RPKs", "SRR062301_Abundance.RPKs",
        "SRR062276_Abundance.RPKs", "SRR1804686_Abundance.RPKs", "SRR1804628_Abundance.RPKs",
        "SRR514191_Abundance.RPKs", "SRR514180_Abundance.RPKs", "SRR513168_Abundance.RPKs",
        "SRR514231_Abundance.RPKs", "SRR513448_Abundance.RPKs")
sal = c("SRR062435_Abundance.RPKs", "SRR062441_Abundance.RPKs", "SRR062389_Abundance.RPKs",
        "SRR062413_Abundance.RPKs", "SRR062402_Abundance.RPKs", "SRR062396_Abundance.RPKs",
        "SRR346673_Abundance.RPKs", "SRR346681_Abundance.RPKs", "SRR062371_Abundance.RPKs",
        "SRR062372_Abundance.RPKs", "SRR062462_Abundance.RPKs", "SRR062415_Abundance.RPKs")
sto = c("SRR528423_Abundance.RPKs", "SRR528353_Abundance.RPKs", "SRR528300_Abundance.RPKs",
        "SRR528261_Abundance.RPKs", "SRR528183_Abundance.RPKs", "SRR528155_Abundance.RPKs",
        "SRR532178_Abundance.RPKs", "SRR532183_Abundance.RPKs", "SRR532190_Abundance.RPKs",
        "SRR532191_Abundance.RPKs", "SRR533152_Abundance.RPKs", "SRR533153_Abundance.RPKs")

orderer <- c(arm, vag, sal, sto)


spec_new <- spec_nums[,orderer]
library(matrixStats)
spec_var <- rowVars(as.matrix(spec_new))
which(spec_var > 0.001)

library(compositions)
spec_trans <- clr(spec_new+1)
#spec_trans <- spec_new
#spec_var <- rowVars(as.matrix(spec_trans))
#keeps <- which(spec_var > 0.1)

cc <- data.frame(spec_trans)
boxplot.matrix(spec_trans)

#cc2 <- cc[keeps,]

#RA <- unlist(spec_new)

x <- sweep(spec_new, 1L, rowMeans(spec_new, na.rm = T), check.margin = FALSE)
sx <- apply(x, 1L, sd, na.rm = T)
x <- sweep(x, 1L, sx, "/", check.margin = FALSE)

heatmap.3(x, dendorgram = 'none', Colv = F, Rowv = F, margin=c(10,3))

RA <- unlist(x)
RA <- unlist(spec_new)
RA_clr <- unlist(cc)
sample_id <- rep(colnames(spec_new), each = 250)
group <- c(rep("arm", 12), rep("vag", 11), rep('sal', 12), rep("sto", 12))
Group <- rep(group, each = 250)
Taxa <- rep(rownames(spec_new), 47)

stat_df <- data.frame(Sample = sample_id, Group = Group, RA = RA,
                      RA_clr = RA_clr, Taxa = Taxa)

head(RA)
head(sample_id)
head(Taxa)
table(Group)
#a_test <- anova(lm(RA_clr ~ Group + Taxa, stat_df))
uniq_bugs <- as.character(unlist(unique(stat_df$Taxa)))
stat_results <- data.frame(Taxa = NA, F_val = NA, P_val = NA)
hoc_results <- data.frame(Taxa = NA, Int = NA, p_adj = NA)
i=194
for (i in 1:length(uniq_bugs)){
  bug <- uniq_bugs[i]
  test <- subset(stat_df, Taxa == bug)
  
  a_test <- anova(lm(RA ~ Group, test))
  result <- data.frame(Taxa = uniq_bugs[i],
                       F_val = a_test$`F value`[1],
                       P_val = a_test$`Pr(>F)`[1])
  stat_results <- rbind(stat_results, result)
  if (result$P_val < 0.05){
    h_test <- aov(RA ~ Group, test)
    tk <- TukeyHSD(h_test, "Group")
    result <- data.frame(Taxa = rep(uniq_bugs[i], ncol(combn(length(unique(test$Group)), 2))),
                         Int = rownames(tk$Group),
                         p_adj = tk$Group[,4])
    hoc_results <- rbind(hoc_results, result)
  }
}
stat_results$p_adj <- p.adjust(stat_results$P_val, method = "BH")
keepers <- subset(stat_results, p_adj < 0.05)$Taxa

hoc2 <- subset(hoc_results, p_adj < 0.05)

keepers <- unique(hoc2$Taxa)

## post hoc?
hoc_df <- subset(stat_df, Taxa %in% keepers)
uniq_bugs <- as.character(unlist(unique(hoc_df$Taxa)))
hoc_results <- data.frame(Taxa = NA, Int = NA, p_adj = NA)

for (i in 1:length(uniq_bugs)){
  bug <- uniq_bugs[i]
  test <- subset(stat_df, Taxa == bug)
  
  h_test <- aov(RA ~ Group, test)
  tk <- TukeyHSD(h_test, "Group")
  result <- data.frame(Taxa = rep(uniq_bugs[i], ncol(combn(length(unique(hoc_df$Group)), 2))),
                       Int = rownames(tk$Group),
                       p_adj = tk$Group[,4])
  hoc_results <- rbind(hoc_results, result)
}


#full hoc
hoc_good <- subset(hoc_results, p_adj < 0.05)
length(unique(hoc_good$Taxa))

##


go_show <- subset(spec_new, rownames(spec_new) %in% keepers)


go_show_nums <- go_show[,-1]
rownames(go_show_nums) = go_show[,1]

go_show_nums2 <- log10(go_show_nums + 1)

#pdf("hmp_heat.pdf", height = 12, width = 10)
heatty = heatmap.3(go_show_nums2, scale = 'row', col = viridis(100),
                   margins = c(12,20), Colv = F)
#dev.off()




#

####### reactive post hoc

## post hoc?
hoc_df <- subset(stat_df, Taxa %in% keepers)
uniq_bugs <- as.character(unlist(unique(hoc_df$Taxa)))
hoc_results <- data.frame(Taxa = NA, Int = NA, p_adj = NA)
i=22
bug <- uniq_bugs[i]
test <- subset(stat_df, Taxa == bug)

h_test <- aov(RA ~ Group, test)
tk <- TukeyHSD(h_test, "Group")

resm <- matrix(NA, 4, 4)
resm[lower.tri(resm) ] <-round(tk$Group[,4], 4)
resm <- t(resm)
resm[lower.tri(resm) ] <-round(tk$Group[,4], 4)
rownames(resm) <- levels(test$Group)
colnames(resm) <- levels(test$Group)
resm

if(any(resm<0.05)){
  print('yay')
}

heatmap.3(resm,
          #cellnote = lab_mat,
          notecol="black",
          main=bug,
          breaks = seq(0, 0.05, by = 0.0005),
          col=c(colorRampPalette(c("red", "white"))(n=100)),
          dendrogram = 'none',
          Rowv=F,
          Colv=F,
          margins=c(10,10),
          cexRow=1.2,
          cexCol=1.2#,
          #na.color="gray60"
)

#
#
#
#

#
#
########
#
#
#
#
#
#
#

##### fix stupid shit

setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)

taxa = fread("ERP103961_taxonomy_abundances_LSU_v4.0.tsv", header=TRUE, sep='\t')
feature_df <- fread("ERP103961_small.tsv", header=T, sep='\t')

tax <- taxa[,-1]

feat <- feature_df[,-1]
feat <- feat[,-1]

length(intersect(colnames(feat), colnames(tax)))
setdiff(colnames(feat), colnames(tax))



## fix super kingdom shit

taxa = fread("ERP103961_taxonomy_abundances_LSU_v4.0.tsv", header=TRUE, sep='\t')

bugs <- taxa$Taxa
bugg_fix <- gsub(";k__", "_k__", bugs)

bug_levels <- strsplit(bugg_fix, ";")

i=200
for (i in 1:length(bug_levels)){
  bugg <- bug_levels[i]
  keep <- length(bugg[[1]])
  bugg[[1]][7]
}





#
#
#
#
#
#
######### differential taxa usign DESeq2

##### differential taxa play
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)
#library(DESeq2)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)

ebi_df <- fread("ERP103961_taxonomy_abundances_LSU_v4.0.tsv", header=T, sep='\t')

ebi_counts <- ebi_df[,-1]
rownames(ebi_counts) <- ebi_df$Taxa



####get the right groups

### this might be the hard part... need spec mat reorder
group1 <- colnames(ebi_counts)[1:100]
group2 <- colnames(ebi_counts)[101:200]
group3 <- colnames(ebi_counts)[201:300]


groups <- c("group1", "group2", "group3")
group_nums <- c(100, 100, 100)
meta <- data.frame(sample_id = colnames(ebi_counts)[1:300],
                   #sample = 1:ncol(ebi_counts))
                   sample = 1:300)

bin <- c()
for (i in 1:length(groups)){
  use <- rep(groups[i], group_nums[i])
  bin <- c(bin, use)
}
meta$bins <- bin

#for (i in 1:length(groups)){
#  new_bin <- as.integer(meta$bins == groups[i])
#  meta <- cbind(meta, new_bin)
#}

#colnames(meta) <- c("sample_id", "sample", "bins", groups)
rownames(meta) <- meta$sample_id


ebi_counts <- data.frame(ebi_counts)
keep_samps <- as.character(unlist(meta$sample_id))

ebi_counts2 <- ebi_counts[,keep_samps]
rownames(ebi_counts2) <- ebi_df$Taxa

meta <- meta[colnames(ebi_counts2),]

## deseq2

countData <- ebi_counts2
condition <- factor(meta$bins)

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~condition)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
diagdds = estimateSizeFactors(dds, geoMeans = geoMeans)

diagdds = DESeq(diagdds, fitType="local")


## loop contrasts
taxa_select <- c("sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Propionibacteriales;f__Propionibacteriaceae;g__Propionibacterium;s__Propionibacterium_acnes")
contrast_results <- data.frame(Int = NULL, P_adj = NULL)
cond_mat <- combn(levels(condition), 2)
i=1
for (i in 1:ncol(cond_mat)){
  use <- cond_mat[,i]
  contrast <- c("condition", use[1], use[2])
  res<-results(diagdds, contrast = contrast)
  cc <- data.frame(res)
  cc2 <- subset(cc, rownames(cc) %in% taxa_select)
  
  interaction <- paste0(use[1], "_", use[2])
  P_adj <- cc2$padj
  results <- data.frame(Int = interaction, P_adj = P_adj)
  contrast_results <- rbind(contrast_results, results)
}

group_num <- length(levels(condition))
group_num <- 2
contrast_results$P_adj2 <- c(0.03)

resm <- matrix(NA, group_num, group_num)
resm[lower.tri(resm) ] <-round(contrast_results$P_adj2, 4)
resm <- t(resm)
resm[lower.tri(resm) ] <-round(contrast_results$P_adj2, 4)
rownames(resm) <- levels(condition)[1:2]
colnames(resm) <- levels(condition)[1:2]
resm

any(resm<0.05, na.rm = T)
breaker <- seq(0, 0.05, by = 0.001)

col_breaks = c(seq(0,0.1,length=5),  
               seq(0.11,0.2,length=5),              
               seq(0.21,0.3,length=5))

coler <- c(colorRampPalette(c("grey90", "white"))(n=10))
heatmap.3(data.matrix(resm), dendrogram = 'none', Colv = F, Rowv = F,
          breaks = breaker, 
          #col = coler, 
          key = F)

#### done with contrasts
res<-results(diagdds, contrast = c("condition", "group1", "group2"))
res<-res[order(res$padj),]
head(res)
cc <- data.frame(res)


### regular
res<-results(diagdds)
res<-res[order(res$padj),]
head(res)
cc1 <- data.frame(res)


go_de <- subset(cc, log2FoldChange > 1.5 | log2FoldChange < -1.5)
go_de_sig <- subset(go_de, padj < 0.05)
keepers <- unique(rownames(go_de_sig))

##add taxa for now
#go_est$bug <- gsub(".*\\|", "", go_df$`Gene Family`)

go_show <- subset(ebi_df[,1:301], Taxa %in% keepers)

go_show_slim <- aggregate(go_show[,-1], list(go_show$Taxa), sum)

go_show_nums <- go_show_slim[,-1]
rownames(go_show_nums) = go_show_slim[,1]

go_show_nums2 <- log10(go_show_nums + 1)

colorss <- meta$bins

cc = length(unique(meta$bins))
cols <- randomColor(3)

for (i in 1:cc){
  colorss = gsub(groups[i], cols[i], colorss)
}

col_col <- t(matrix(colorss))
col_col <- rbind(col_col, col_col)

#pdf("hmp_heat.pdf", height = 12, width = 10)
heatty = heatmap.3(t(t(go_show_nums2)), scale = 'row', col = viridis(100),
                   margins = c(5,25), Colv = F, labCol = "", 
                   ColSideColors = t(col_col), side.height.fraction = 0.5)
#dev.off()





#

#
#
#
#
#
#
#
#
#
#
######## practice heat with bar and DA by anova

### diff abundance by ANOVA
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)


go_df <- fread("Relab_go_names_wham.tsv", header=T, sep='\t')

go_nums <- go_df[,4:50]

## spec da

feat_df <- aggregate(go_nums, list(go_df$Feature), sum)

feat_nums <- feat_df[,2:48]
rownames(feat_nums) <- feat_df[,1]

arm = c("SRR532024_Abundance.RPKs", "SRR532015_Abundance.RPKs", "SRR532006_Abundance.RPKs",
        "SRR532040_Abundance.RPKs", "SRR532504_Abundance.RPKs", "SRR532507_Abundance.RPKs",
        "SRR638753_Abundance.RPKs", "SRR640357_Abundance.RPKs", "SRR640340_Abundance.RPKs",
        "SRR640452_Abundance.RPKs", "SRR640499_Abundance.RPKs", "SRR545546_Abundance.RPKs")
vag = c("SRR062353_Abundance.RPKs", "SRR062357_Abundance.RPKs", "SRR062301_Abundance.RPKs",
        "SRR062276_Abundance.RPKs", "SRR1804686_Abundance.RPKs", "SRR1804628_Abundance.RPKs",
        "SRR514191_Abundance.RPKs", "SRR514180_Abundance.RPKs", "SRR513168_Abundance.RPKs",
        "SRR514231_Abundance.RPKs", "SRR513448_Abundance.RPKs")
sal = c("SRR062435_Abundance.RPKs", "SRR062441_Abundance.RPKs", "SRR062389_Abundance.RPKs",
        "SRR062413_Abundance.RPKs", "SRR062402_Abundance.RPKs", "SRR062396_Abundance.RPKs",
        "SRR346673_Abundance.RPKs", "SRR346681_Abundance.RPKs", "SRR062371_Abundance.RPKs",
        "SRR062372_Abundance.RPKs", "SRR062462_Abundance.RPKs", "SRR062415_Abundance.RPKs")
sto = c("SRR528423_Abundance.RPKs", "SRR528353_Abundance.RPKs", "SRR528300_Abundance.RPKs",
        "SRR528261_Abundance.RPKs", "SRR528183_Abundance.RPKs", "SRR528155_Abundance.RPKs",
        "SRR532178_Abundance.RPKs", "SRR532183_Abundance.RPKs", "SRR532190_Abundance.RPKs",
        "SRR532191_Abundance.RPKs", "SRR533152_Abundance.RPKs", "SRR533153_Abundance.RPKs")

orderer <- c(arm, vag, sal, sto)


feat_new <- feat_nums[,orderer]

library(matrixStats)
library(compositions)
feat_trans <- clr(feat_new)
feat_var <- rowVars(as.matrix(feat_trans))
keeps <- which(feat_var > 0.1)

cc <- data.frame(feat_trans)
boxplot.matrix(feat_trans)

#cc2 <- cc[keeps,]

RA <- unlist(feat_new)
RA_clr <- unlist(cc)
sample_id <- rep(colnames(feat_nums), 3706)
group <- c(rep("arm", 12), rep("vag", 11), rep('sal', 12), rep("sto", 12))
Group <- rep(group, 3706)
Feature <- rep(rownames(cc), each = 47)

stat_df <- data.frame(Sample = sample_id, Group = Group, RA = RA,
                      RA_clr = RA_clr, Feature = Feature)


#a_test <- anova(lm(RA_clr ~ Group + Taxa, stat_df))
uniq_feat <- as.character(unlist(unique(stat_df$Feature)))
stat_results <- data.frame(Feature = NA, F_val = NA, P_val = NA)
for (i in 1:length(uniq_feat)){
  feat <- uniq_feat[i]
  test <- subset(stat_df, Feature == feat)
  
  a_test <- anova(lm(RA ~ Group, test))
  result <- data.frame(Feature = uniq_feat[i],
                       F_val = a_test$`F value`[1],
                       P_val = a_test$`Pr(>F)`[1])
  stat_results <- rbind(stat_results, result)
}

keepers <- subset(stat_results, P_val < 0.05)$Feature

go_show <- subset(feat_df, Group.1 %in% keepers)


go_show_nums <- go_show[,-1]
rownames(go_show_nums) = go_show[,1]

go_show_nums2 <- log10(go_show_nums + 1)

#pdf("hmp_heat.pdf", height = 12, width = 10)
heatty = heatmap.3(go_show_nums2, scale = 'row', col = viridis(100),
                   margins = c(12,20), Colv = F) 
                   #Colv = T,
                   #labCol = "", labRow = "")


#dev.off()
gg_nums <- go_show_nums2

gg_nums_feat <- gg_nums[heatty$rowInd,]
gg_nums_samp <- gg_nums_feat[,heatty$colInd]

gg_names <- factor(colnames(gg_nums_samp), levels = colnames(gg_nums_samp))
gg_feat <- factor(rownames(gg_nums_samp), levels = rownames(gg_nums_samp))

x <- sweep(gg_nums_samp, 1L, rowMeans(gg_nums_samp, na.rm = T), check.margin = FALSE)
sx <- apply(x, 1L, sd, na.rm = T)
x <- sweep(x, 1L, sx, "/", check.margin = FALSE)

## for taxa plot later
gene_order= hclust(dist(x))
plot(gene_order)
###

Z_scored_CPM <- unlist(x)
head(Z_scored_CPM)


names <- rep(gg_names, each = nrow(gg_nums))
feat <- rep(gg_feat, ncol(gg_nums))
#groups <- rep(bin, nrow(go_show_nums2))

plotter <- data.frame(Z_scored_CPM, names, feat)
###
##
high = viridis(n = 100)[100]
low = viridis(n = 100)[1]
mid = viridis(n = 100)[50]

pal = viridis(n=100)

gMain <- ggplot(plotter, aes(x=names, y = feat, fill = Z_scored_CPM)) +
  geom_tile(aes(fill = Z_scored_CPM)) +
  ylab("") + xlab("")+
  #scale_fill_continuous() +
  #scale_fill_viridis(discrete=F, limits = c(-6,6), breaks = seq(-6,6,1)) +
  scale_fill_gradientn(colors = pal,
                       limits = c(-6,6), breaks = seq(-6,6,1)) +
  #scale_fill_gradient2(low = 'blue', high = 'red', mid = "white",
  #                     limits = c(-6,6), breaks = seq(-6,6,1)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


### make gmain in plotly!!!!
library(plotly)

ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  tickcolor = 'white',
  showgrid = FALSE
)
p <- plot_ly(
  x = plotter$names, y = plotter$feat,
  z = plotter$Z_scored_CPM, type = "heatmap",
  colorbar = list(x = -0.2, 
                  xanchor = 'left',
                  tickmode='array',
                  tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
                  ticktext = c("low", "high"))) %>%
  layout(yaxis = ax, xaxis = ax)
p


##taxa plot###
go_da_specs <- subset(go_df, Feature %in% keepers)


taxa_nums <- data.matrix(go_da_specs[,4:50])
rownames(taxa_nums) = paste0(go_da_specs$Feature, ":", go_da_specs$Taxa)

data2 <- log2(taxa_nums+1)

log_abundance <- c(data2)
genes2 <- rep(go_da_specs$Feature, ncol(taxa_nums))
bugs2 <- rep(go_da_specs$Taxa, ncol(taxa_nums))
samples2 <- rep(colnames(taxa_nums), each = nrow(taxa_nums))
#groups <- c(rep("SF", 17), rep("NYU", 6))
#groups2 <- rep(groups, each = nrow(taxa_nums))
plotter2 <- data.frame(genes2, bugs2, samples2, log_abundance)

#samp_order= hclust(dist(t(taxa_nums)))
#gene_order= hclust(dist(taxa_nums))
length(unique(plotter2$genes2))
#plotter$samples <- factor(plotter$samples, levels = samp_order$labels)
plotter2$genes2 <- factor(plotter2$genes2, levels = gene_order$labels)
#plotter$groups <- rep(c(rep("SF", 16), c(rep("NYU", 5))), each = 26)

plotter2$genera <- gsub("\\.s.*", "", plotter2$bugs2)

g_small <- aggregate(plotter2$log_abundance, list(plotter2$genera), mean, na.rm=T)
g_small <- g_small[order(g_small$x, decreasing=T),]
g_high <- g_small[1:10,]
plotter2 <- subset(plotter2, genera %in% g_high$Group.1)
plotter2$genera <- factor(plotter2$genera, levels = g_high$Group.1)



#pdf("figures/DEgo_term_enrichment_by_taxa.pdf", height = 10, width = 13)
gTaxa <- ggplot(data = plotter2, aes(x=genes2, y=log_abundance, fill=genera)) +
  stat_summary(aes(genes2), fun.y = 'mean', geom='bar', 
               position = 'fill') +
  ylab("Relative Abundance") +
  xlab("GO Term") +
  guides(fill=guide_legend(ncol=5, title = "Genera")) +
  scale_fill_manual(values = randomColor(length(unique(plotter2$genera)))) +
  #coord_flip() +
  theme(panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_line(color = "white"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = 'top')




# library(gridExtra)
# grid.arrange(gTaxa, gMain, ncol = 1)


#main_plotly <- ggplotly(gMain)


taxa_plotly <- ggplotly(gTaxa)
taxa_plotly$x$layout$showlegend = FALSE
#taxa_plotly$x$data[[10]]

####

p2 <- subplot(style(taxa_plotly, showlegend = FALSE),
              p, nrows = 2, margin = c(0,0,-0.01,0),
              heights = c(0.2, 0.8))
p2

p2$width = 1400
p2$height = 1000

export(p2, "file.png")


##### make fake color bar to add to the plotly object!
greenSeriesData <- matrix(c(1,1,NA,NA), nrow = 2)
redSeriesData <- matrix(c(NA,NA,1,1), nrow = 2)

greenColor <- data.frame(x = c(1,2), y = c("#63a803", "#63a803"))
colnames(greenColor) <- NULL

redColor <- data.frame(x = c(1,2), y = c("#a80b03", "#a80b03"))
colnames(redColor) <- NULL


colorer <- plot_ly(
  type = "heatmap"
) %>% add_trace(
  z = greenSeriesData,
  colorscale = greenColor, showscale = F
) %>% add_trace(
  z = redSeriesData,
  colorscale = redColor, showscale = F
) %>%
  layout(yaxis = ax, xaxis = ax)


### make real fake color bar

names <- c(rep("arm", 12), rep("vag", 11), rep("sal", 12), rep("sto", 12))
uniq_names <- unique(names)

meta <- data.frame(names)

for (i in 1:length(uniq_names)){
  new_bin <- as.integer(meta$names == uniq_names[i])
  meta <- cbind(meta, new_bin)
}

series_mat <- meta[,-1]
colnames(series_mat) <- uniq_names
#series_mat <- t(series_mat)
series_mat[series_mat==0] <- NA

cols_cols <- c("red", "blue", "orange", "purple")

g1 = t(as.matrix(series_mat[,1]))
g1col <- data.frame(x = c(0,1), y = c(cols_cols[1], cols_cols[1]))
colnames(g1col) <- NULL

my_text <- gsub(1, "arm", g1)

colorer <- plot_ly(
  type = "heatmap"
) %>% add_trace(
  x=1:47,
  z = g1,
  colorscale = g1col, showscale = F,
  hoverinfo = 'all'
) #%>%
  #layout(yaxis = ax, xaxis = ax)

pdf("check_it.pdf")
print(colorer)
dev.off()

  
for (i in 2:length(uniq_names)){
  g = t(as.matrix(series_mat[,i]))
  gcol <- data.frame(x = c(0,1), y = c(cols_cols[i], cols_cols[i]))
  colnames(gcol) <- NULL
  
  colorer <- colorer %>% add_trace(
    z = g,
    colorscale = gcol, showscale = F,
    hoverinfo = 'all'
  )
}

colorer
colorer$x$attrs$`16856574d717`$hoverinfo
####
p3 <- subplot(colorer,
              p, nrows = 2, margin = c(0,0,-0.01,0),
              heights = c(0.1, 0.9), shareX = T)
p3



### 
cc <- g_legend(gTaxa)
plot(cc)

class(p3)
class(cc)

class(gTaxa)



#
#
#
#

#
#
#
#
#
#####
#
#
####### taxa plot that responds to clicks

setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)


go_df <- fread("Relab_go_names_wham.tsv", header=T, sep='\t')

full_select <- data.frame(subset(go_df, Feature %in% "[BP] sulfate assimilation"))

go_nums <- full_select[,4:50]

sumin <- rev(1:47)
go_reorder <- full_select[,sumin]

RA <- unlist(go_nums)
Sample_id <- rep(colnames(go_nums), each = nrow(go_nums))
Taxa <- rep(full_select$Taxa, ncol(go_nums))
Group <- rep(c(rep("arm", 12), rep("vag", 11), rep('sal', 12), rep("sto", 12)), each = nrow(go_nums))

plotter <- data.frame(Sample_id, Taxa, Group, RA)

bug_sorter <- plotter[order(plotter$RA, decreasing = T),]
bug_sort <- as.character(unique(bug_sorter$Taxa))
plotter$Taxa <- factor(plotter$Taxa, levels = bug_sort)

tax_sel <- ggplot(plotter, aes(x=Group, y = RA, fill = Taxa)) + 
  stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
  scale_fill_manual(values = randomColor(length(unique(plotter$Taxa)))) +
  theme(legend.position='none',
        panel.background = element_blank())

taxly <- ggplotly(tax_sel, tooltip = c("x", 'fill'))
taxly



#
#
#
#
#
#
#
#
#
#
#
#
#### colorer legend for all
Group <- c('arm', 'leg', 'thing')
coll <- c('red', 'blue', 'green')

plotplot <- data.frame(Group)

gg = ggplot(plotplot, aes(1:3, 1:3, color = Group)) +
  geom_point() +
  guides(color = guide_legend(direction = "horizontal",
                              override.aes = list(shape = 15, size=9)))+
  scale_color_manual(values = coll) + 
  theme(legend.key = element_rect(fill = "transparent"),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 42))


g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(gg) 
plot(legend)
library(gridExtra)
grid.arrange(legend)




#
#
#
#

#
#
#
#
#
#
#
#
#
### pick nice colors

library(randomcoloR)

cols <- randomColor(10, hue = 'red')

pp <- data.frame(x=1:50, y = 1:50, col = rep(factor(c(letters[1:10])), each = 5))

ggplot(pp, aes(x=x, y=y, color = col)) +
  geom_point() +
  scale_color_manual(values=cols_keep2)

cols[6]

cols_keep <- c("#5beabd", "#a0080f", "#6d89e8", "#482e91", "#91dd68",
               "#fcdb23", "#fc9207", "#01871e", "#e87ac3", "#f70c1c")


cols_keep2 <- c("#f70c1c", "#6d89e8", "#91dd68", "#482e91", "#fc9207",
                "#fcdb23", "#e87ac3", "#5beabd", "#01871e", "#a0080f")




#
#
#

#
#
#
#
#
#
#


x <- data.frame(c(0.03, NA),
                 c(NA, 0.03))
x_mat <- data.matrix(x)
x2 <- cbind(x_mat[,1], x_mat[,1])

heatmap.2(x2, Colv=F, Rowv = F, na.rm = T)


library(gridExtra)
d <- head(iris, 3)
g <- tableGrob(x_mat)
g$heights
grid.newpage()
grid.arrange(g)







### percentile play

cc <- runif(200, 0,100)
quantile(cc, 0.5)

max(cc)




title("Your title", line= -14.5)

c1 <- c(0.03, 1)
c2 <- c(4, 0.03)
resm <- rbind(c1, c2)

resm[1,2] = 0.03 + 0.0000001




#
#
#

#
######## 

x <- replicate(10, rnorm(20)) 
heatmap.3(x)

cc <- hclust(dist(x))
plot(cc)

cc$order


cc <- data.frame(clr(x))
unlist(cc)

heatmap.3(data.matrix(cc))



#
#
#
#
#
#
#
#
#
#
#

#
######### deseq2 concerns

##### differential taxa play
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)
library(metagenomeSeq)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)

ebi_df <- fread("ERP104179_IPR_abundances_v4.1.tsv", header=T, sep='\t')

ebi_counts <- data.frame(ebi_df[,-1])
ebi_counts <- data.frame(ebi_counts[,-1])
rownames(ebi_counts) <- ebi_df$Feature

ebi_counts[is.na(ebi_counts)]

g_maker <- hclust(dist(t(ebi_counts)))
plot(g_maker)

ord = colnames(ebi_counts)[g_maker$order]

ebi_order <- data.frame(ebi_counts)[,ord]

ebi_relab <- sweep(ebi_order, 2, colSums(ebi_order), '/')
rownames(ebi_relab) <- rownames(ebi_counts)



####get the right groups






#
#
#
#
#
#
#
## colorer with one group
series_mat <- rep(1, 20)

names(series_mat) <- rep("Group1", 20)
#series_mat <- t(series_mat)
series_mat[series_mat==0] <- NA

#cols_cols <- randomColor(length(uniq_names))

g1 = t(as.matrix(series_mat))
g1col <- data.frame(x = c(0,1), y = c('red', 'red'))
colnames(g1col) <- NULL

ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  tickcolor = 'white',
  showgrid = FALSE
)

colorer <- plot_ly(
  type = "heatmap"
) %>% add_trace(
  x=1:20,
  z = g1,
  colorscale = g1col, showscale = F,
  hoverinfo = 'all'
) %>%
  layout(yaxis = ax, xaxis = ax)
}
colorer





#####




#
#
#
#
#### all zero matrix play

mat <- rbind(c(NA, 0, 0),
             c(0,NA,0),
             c(0,0,NA))
if (sum(mat, na.rm = T) == 0){
  mat[1,2] = 1e-7
}
heatmap.3(mat)




#
#

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##
#
#
#
#
#
#
#
#
#
#
#
#
#
######
#
#
#
#
###### ALDEX tests and estimating HUMANN2 relab for count data!!!!
setwd("/Users/jcooperdevlin/Desktop/Kelly/wham_revisions3.28")
library(data.table)
source("/Users/jcooperdevlin/Desktop/Itai/ST_app/heater3/R/heatmap.3.R")
library(viridis)
library(ggplot2)
library(randomcoloR)


go_df <- fread("Relab_go_names_wham.tsv", header=T, sep='\t')
go_nums <- go_df[,4:50]

total_reads_df <- read.table("raw_reads.csv", header=T, sep='\t')[,-1]
total_reads_df$Sample <- gsub(".fastq", "_Abundance.RPKs", total_reads_df$Sample)

rownames(total_reads_df) = total_reads_df$Sample
total_reads_df <- total_reads_df[colnames(go_nums),]
total_reads_df <- t(total_reads_df$Total_Reads)

go_count_est <- sweep(go_nums, 2, total_reads_df, "*")
go_count_est <- round(go_count_est, 0)

go_count_est2 <- cbind(go_df[,1:3], go_count_est)
write.table(go_count_est2, "HMP_count_est_go_names_wham.tsv", row.names=F, sep='\t',
            quote=F)

### try aldex stats
feat_df <- aggregate(go_count_est, list(go_df$Taxa), sum)

feat_nums <- feat_df[,2:48]
rownames(feat_nums) <- feat_df[,1]

arm = c("SRR532024_Abundance.RPKs", "SRR532015_Abundance.RPKs", "SRR532006_Abundance.RPKs",
        "SRR532040_Abundance.RPKs", "SRR532504_Abundance.RPKs", "SRR532507_Abundance.RPKs",
        "SRR638753_Abundance.RPKs", "SRR640357_Abundance.RPKs", "SRR640340_Abundance.RPKs",
        "SRR640452_Abundance.RPKs", "SRR640499_Abundance.RPKs", "SRR545546_Abundance.RPKs")
vag = c("SRR062353_Abundance.RPKs", "SRR062357_Abundance.RPKs", "SRR062301_Abundance.RPKs",
        "SRR062276_Abundance.RPKs", "SRR1804686_Abundance.RPKs", "SRR1804628_Abundance.RPKs",
        "SRR514191_Abundance.RPKs", "SRR514180_Abundance.RPKs", "SRR513168_Abundance.RPKs",
        "SRR514231_Abundance.RPKs", "SRR513448_Abundance.RPKs")
sal = c("SRR062435_Abundance.RPKs", "SRR062441_Abundance.RPKs", "SRR062389_Abundance.RPKs",
        "SRR062413_Abundance.RPKs", "SRR062402_Abundance.RPKs", "SRR062396_Abundance.RPKs",
        "SRR346673_Abundance.RPKs", "SRR346681_Abundance.RPKs", "SRR062371_Abundance.RPKs",
        "SRR062372_Abundance.RPKs", "SRR062462_Abundance.RPKs", "SRR062415_Abundance.RPKs")
sto = c("SRR528423_Abundance.RPKs", "SRR528353_Abundance.RPKs", "SRR528300_Abundance.RPKs",
        "SRR528261_Abundance.RPKs", "SRR528183_Abundance.RPKs", "SRR528155_Abundance.RPKs",
        "SRR532178_Abundance.RPKs", "SRR532183_Abundance.RPKs", "SRR532190_Abundance.RPKs",
        "SRR532191_Abundance.RPKs", "SRR533152_Abundance.RPKs", "SRR533153_Abundance.RPKs")

orderer <- c(arm, vag, sal, sto)


feat_new <- feat_nums[,orderer]

arm_vag <- feat_new[,24:47]

library(ALDEx2)

conds <- c(rep("arm", 12), rep("vag", 11), rep("sal", 12), rep("sto", 12))
#x <- aldex.clr(feat_new, conds, mc.samples=16, verbose=TRUE)

#x.glm <- aldex.glm(x, conds)

combos <- combn(unique(conds), 2)
stat_results = data.frame(we.ep = NA, we.eBH = NA, wi.ep = NA, wi.eBH = NA,
                          rab.all = NA, rab.win.a = NA, rab.win.b = NA,
                          diff.btw = NA, diff.win = NA, effect = NA, overlap = NA,
                          Int = NA, Feature = NA)
for (i in 1:ncol(combos)){
  ab <- combos[,i]
  a <- ab[1]
  b <- ab[2]
  a_sel = which(conds == a)
  b_sel = which(conds == b)
  
  feat_subset <- cbind(feat_new[,a_sel], feat_new[,b_sel])
  cond_subset <- c(conds[a_sel], conds[b_sel])
  
  x <- aldex.clr(feat_subset2, cond_subset, mc.samples=16, verbose=TRUE)
  #x.glm <- aldex.glm(x, cond_subset)
  x.tt <- aldex.ttest(x, cond_subset, paired.test = F)
  x.effect <- aldex.effect(x, cond_subset, include.sample.summary=FALSE, verbose=TRUE)
  x.all <- data.frame(x.tt, x.effect, stringsAsFactors=FALSE)
  x.all$Int <- rep(paste0(a,"_", b), nrow(x.all))
  x.all$Feature <- rownames(x.all)
  names(x.all) = names(stat_results)
  stat_results <- rbind(stat_results, x.all)
}


aldex.plot(x.all, type="MA", test="welch")

aldex.plot(x.all, type="MW", test="welch")



##
#
#######
haves <- c("arm-vag", "arm-sal", "arm-sto", "sal-sto")
total <- c("arm-vag", "arm-sal", "arm-sto", "vag-sal", "vag-sto", "sal-sto")

setdiff(total, haves)



#
#
#
#
#### 1-D aldex?
library(ALDEx2)
library(SummarizedExperiment)
library(GenomicRanges)
library(IRanges)
library(Biobase)
library(BiocGenerics)
library(multtest)


wilcox.test(t.input[2,1:25], t.input[2,26:50])


names <- paste0("g_", 1:50)
nums3 <- abs(round(c(rnorm(25), rnorm(25,100)), 0))
nums <- rbind(nums,nums2, nums3)

h <- data.frame(nums)
h <- t(h)

names(h) <- names
conds <- c(rep('a', 25), rep('b', 25))

hh <- data.frame(h)
names(hh) <- names
#hh2 <- rbind(hh, hh)

x <- aldex.clr(hh, conds = conds, mc.samples = 16)
#x.tt <- try(aldex.ttest(x, conds, paired.test = T))
aldex.glm(x, conds)
x.effect = aldex.effect(x, conditions = conds)

######ttets debug?


clr <- x
conditions <- conds
paired.test = F
hist.plot = F
aldex.ttest2 <- function (clr, conditions, paired.test = FALSE, hist.plot = FALSE) 
{
  smpl.ids <- getSampleIDs(clr)
  feature.names <- getFeatureNames(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)
  conditions <- as.factor(conditions)
  levels <- levels(conditions)
  if (length(conditions) != numConditions(clr)) {
    stop(paste("mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions):", 
               length(conditions), "len(names(clr)):", numConditions(clr)))
  }
  if (length(levels) != 2) 
    stop("only two condition levels are currently supported")
  levels <- vector("list", length(levels))
  names(levels) <- levels(conditions)
  sets <- names(levels)
  setAsBinary <- as.numeric(conditions == sets[1])
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  wi.p.matrix <- as.data.frame(matrix(1, nrow = feature.number, 
                                      ncol = mc.instances))
  wi.BH.matrix <- wi.p.matrix
  we.p.matrix <- wi.p.matrix
  we.BH.matrix <- wi.p.matrix
  print("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for (mc.i in 1:mc.instances) {
    numTicks <- progress(mc.i, mc.instances, numTicks)
    t.input <- sapply(mc.all, function(y) {
      y[, mc.i]
    })
    for (j in 1:nrow(t.input)){
      wi.p.matrix[j, mc.i] <- wilcox.test(t.input[j,1:length(setA)], t.input[j,c(length(setA)+1):c(length(setA)+length(setB))], 
                                         paired = paired.test)$p.value
      
      wi.BH.matrix[j, mc.i] <- p.adjust(wi.p.matrix[j, mc.i], 
                                       method = "BH")
      
      #we.p.matrix[, mc.i] <- t.fast(t.input, setAsBinary, paired.test)
      
      we.p.matrix[j, mc.i] <- t.test(t.input[j,1:length(setA)], t.input[j,c(length(setA)+1):c(length(setA)+length(setB))], 
                                    paired = paired.test)$p.value
      
      we.BH.matrix[j, mc.i] <- p.adjust(we.p.matrix[j, mc.i], 
                                       method = "BH")
    }
  }
  if (hist.plot == TRUE) {
    par(mfrow = c(2, 2))
    hist(we.p.matrix[, 1], breaks = 99, main = "Welch's P values Instance 1")
    hist(wi.p.matrix[, 1], breaks = 99, main = "Wilcoxon P values Instance 1")
    hist(we.BH.matrix[, 1], breaks = 99, main = "Welch's BH values Instance 1")
    hist(wi.BH.matrix[, 1], breaks = 99, main = "Wilcoxon BH values Instance 1")
    par(mfrow = c(1, 1))
  }
  we.ep <- rowMeans(we.p.matrix)
  we.eBH <- rowMeans(we.BH.matrix)
  wi.ep <- rowMeans(wi.p.matrix)
  wi.eBH <- rowMeans(wi.BH.matrix)
  z <- data.frame(we.ep, we.eBH, wi.ep, wi.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}




#
#
#
#

#
#
#
#
#

#
#
#####


#

#
#
#
#
#
### samp dfs
Acc <- paste0("0000", 1:5)
Features <- paste0("Feature", 1:5)
Taxa <- c(paste0("Bug", 1:3), paste0("Bug", 1:2))
num_mat <- matrix(NA,5,5)
for (i in 1:ncol(num_mat)){
  num_mat[,i] <- sample(1000,5)
}
colnames(num_mat) = paste0("Sample", 1:5)
sample_df <- cbind(Acc, Features, Taxa, num_mat)


##ebi
Acc <- paste0("0000", 1:5)
Features <- paste0("Feature", 1:5)
#Taxa <- c(paste0("Bug", 1:3), paste0("Bug", 1:2))
num_mat <- matrix(NA,5,5)
for (i in 1:ncol(num_mat)){
  num_mat[,i] <- sample(1000,5)
}
colnames(num_mat) = paste0("Sample", 1:5)
sample_df <- cbind(Acc, Features, num_mat)





#

#
#
#
#
#
#
#
#
#
#
#
#
#### practice colorkey manips
resm <- matrix(NA, 3,3)
resm[2,1] = 0.03
resm[1,3] = 0.1
#resm[1,1] = 0
#resm[3,3] = 1

resm_lab <- resm
resm_lab[resm_lab < 0.0001] <- "<0.0001"


if (any(resm<0.05, na.rm = T)){
  breaker = seq(0, 0.05, by = 0.0005)
  coler = c(colorRampPalette(c("red", "white"))(n=100))
} else {
  breaker <- seq(0, 1, by = 0.01)
  coler <- c(colorRampPalette(c("white"))(n=100))
}

breaker <- seq(0, 1, by = 0.005)
coler <- c(c(colorRampPalette(c("red", "white"))(n=11)),
           c(colorRampPalette(c("white"))(n=189)))

select_feat <- "Stuff"


par(cex.main=0.8)
heatmap.2(data.matrix(resm),
          cellnote = resm_lab,
          notecol="black",
          key.title = NULL,
          key = T,
          sepcolor="black",
          breaks = breaker,
          col = coler,
          dendrogram = 'none',
          Rowv=F,
          Colv=F,
          margins=c(5,5),
          cexRow=1.2,
          cexCol=1.2#,
)
title(select_feat, line= -2.5)


#####




img_names <- read.table("/Users/jcooperdevlin/Desktop/img_names.txt",
                        header=F, sep='\t')
keeps <- read.table("/Users/jcooperdevlin/Desktop/keeper_table_full4.21.txt",
                    header=F, sep='\t')

img_names2 <- gsub("_Abundance.RPKs.tsv", "", img_names$V1)

add_names <- setdiff(img_names2, keeps$V1)

###

rownames(keeps) <- paste0(keeps$V1, "_Abundance.RPKs.tsv")
keeps2 <- keeps[img_names$V1,]

keeps3 <- subset(keeps2, !is.na(V1))

keeps3$filename <- rownames(keeps3)


####
keeps4 <- data.frame(V1 = NA, V2=NA, V3 = NA, V4 = NA, filename = NA)
for (i in 1:length(add_names)){
  curr_name <- add_names[i]
  if (curr_name != "feature_index.txt"){
  base <- gsub("_.*", "", curr_name)
  item <- subset(keeps3, V1 == base)
  item$filename <- curr_name
  keeps4 <- rbind(keeps4, item)
  }
}




keeps5 <- rbind(keeps3, keeps4[-1,])


