# Path settings
plink_bin <- file.path("plink", "plink")
stringname3 <- file.path("reference", "g1000_eur")
block_file <- file.path("reference", "all.bed")
data_dir <- file.path("data")
temp_dir <- file.path("temp")
output_dir <- file.path("output")
beta_se_dir <- file.path(output_dir, "beta.se")
gam3seo_res_dir <- file.path(output_dir, "MERLIN.res")

library(MERLIN)
library(readr)
library(parallel)
library(MR.LDP)
library(mr.raps)
library(Rcpp)
library(PDSCE)
library(data.table)
library(MendelianRandomization)

bim <- read.table(paste0(stringname3, ".bim")) 

phe.list <- c("MDD","MDDR","ALCD","OCD","BD","SCZ","ASD","PTSD")

exp <- "Testosterone" # Testosterone/Oestradiol
expgwas <- read_delim(file.path(data_dir, paste0(exp,".GWAS.txt")), col_names=TRUE)
expgwas <- data.frame(expgwas)
expgwis <- read_delim(file.path(data_dir, paste0(exp,".GWIS.txt")), col_names=TRUE) 
expgwis <- data.frame(expgwis)
expgwis$BETA <- expgwis$BETA/2
expgwis$SE <- expgwis$SE/2

process_data <- function(data) {
tryCatch({
  
# phe <- phe.list[3]
phe <- data[1]

### For IVW, Raps, Egger
expgwas.LD.cmd <- paste(plink_bin, " --bfile ", stringname3, 
                        " --clump-p1 0.00000005 --clump-r2 0.01 --clump-kb 1024 --maf 0.05 --clump ", 
                        file.path(data_dir, paste0(exp,".GWAS.txt")), 
                        " --clump-snp-field SNP --clump-field P --out ", 
                        file.path(temp_dir, paste0("Tes1.", exp,".",phe)), sep="")
system(expgwas.LD.cmd)
expgwas.LD <- fread(paste0(file.path(temp_dir, paste0("Tes1.",exp,".",phe)), ".clumped"))
snp.expgwas.LD <- expgwas.LD$SNP
length(snp.expgwas.LD)
file.remove(paste0(file.path(temp_dir, paste0("Tes1.",exp,".",phe)), ".clumped"))
file.remove(paste0(file.path(temp_dir, paste0("Tes1.",exp,".",phe)), ".log"))

outgwas <- read_delim(file.path(data_dir, paste0(phe,".GWAS.txt")), col_names=TRUE) 
outgwas <- data.frame(outgwas)

intersect.snp1 <- intersect(snp.expgwas.LD, outgwas$SNP) 
print(paste0(exp,"-",phe," cutoff1 causal snps:",length(intersect.snp1)))

expgwas.use <- expgwas[match(intersect.snp1, expgwas$SNP), ]
expgwas.use <- expgwas.use[complete.cases(expgwas.use), ]
outgwas.use <- outgwas[match(intersect.snp1, outgwas$SNP), ]
outgwas.use <- outgwas.use[complete.cases(outgwas.use), ]

bh14ld <- expgwas.use$BETA
s124ld <- expgwas.use$SE
bh24ld <- outgwas.use$BETA
s224ld <- outgwas.use$SE

res.ivw <- TwoSampleMR::mr_ivw(bh14ld, bh24ld, s124ld, s224ld)
res.raps <- mr.raps.simple(bh14ld, bh24ld, s124ld, s224ld)
res.egger <- mr_egger(mr_input(bx = bh14ld, bxse = bh24ld, by = s124ld, byse = s224ld))

### For MR.LDP, MRGEI_dropgam3seo
expgwas.LD.cmd <- paste(plink_bin, " --bfile ", stringname3, 
                        " --clump-p1 0.00000005 --clump-r2 0.3 --clump-kb 1024 --maf 0.05 --clump ", 
                        file.path(data_dir, paste0(exp,".GWAS.txt")), 
                        " --clump-snp-field SNP --clump-field P --out ", 
                        file.path(temp_dir, paste0("Tes2.", exp,".",phe)), sep="")
system(expgwas.LD.cmd)
expgwas.LD <- fread(paste0(file.path(temp_dir, paste0("Tes2.",exp,".",phe)), ".clumped"))
snp.expgwas.LD <- expgwas.LD$SNP
length(snp.expgwas.LD)
file.remove(paste0(file.path(temp_dir, paste0("Tes2.",exp,".",phe)), ".clumped"))
file.remove(paste0(file.path(temp_dir, paste0("Tes2.",exp,".",phe)), ".log"))

intersect.snp2 <- intersect(snp.expgwas.LD, outgwas$SNP) 
print(paste0(exp,"-",phe," cutoff2 causal snps:",length(intersect.snp2)))

avbIndex <- match(intersect.snp2, bim$V2)
avbIndex <- as.matrix(avbIndex[order(avbIndex)])
intersect.snp2 <- bim[avbIndex, ]$V2

expgwas.order <- expgwas[match(intersect.snp2, expgwas$SNP), ]
outgwas.order <- outgwas[match(intersect.snp2, outgwas$SNP), ]

bh14ld <- expgwas.order$BETA
s124ld <- expgwas.order$SE
bh24ld <- outgwas.order$BETA
s224ld <- outgwas.order$SE

lam <- 0.1; coreNum <- 1;
bp <- expgwas.order$BP
chr <- expgwas.order$CHR
idx4panel <- matrix(numeric(0), nrow = 0, ncol = 1)

Rblockres <- Cal_block_Rmatrix(bp, chr, avbIndex-1, idx4panel, block_file, stringname3, coreNum, lam);
R <- Rblockres$R; diag(R) <- 1 

# LDP
p <- length(bh14ld); gamma <- rep(0.01,p)
alpha <- rep(0.01,p); sgga2 <- 0.01
sgal2 <- 0.01; maxIter <- 10000; 
diagnostics <- FALSE
beta0 <- 0.1; epsStopLogLik <- 1e-6

SimMRLDP_Hb <- MRLDP_SimPXvb(bh14ld, bh24ld, s124ld, s224ld, gamma, alpha,  beta0, sgga2, sgal2, R,
                             0, epsStopLogLik, maxIter, model = 2)
SimMRLDP_H0 <- MRLDP_SimPXvb(bh14ld, bh24ld, s124ld, s224ld, gamma, alpha,  beta0, sgga2, sgal2, R,
                             1, epsStopLogLik, maxIter, model = 2)
tstatLDP <- 2*(SimMRLDP_Hb$tstat - SimMRLDP_H0$tstat)
beta_hatLDP <- SimMRLDP_Hb$beta0;
se_hatLDP <- abs(beta_hatLDP/sqrt(tstatLDP));
pvalLDP <- pchisq(tstatLDP, 1, lower.tail = F)

### For MERLIN
outgwis <- read_delim(file.path(data_dir, paste0(phe,".GWIS.txt")), col_names=TRUE) 
outgwis <- data.frame(outgwis)
outgwis$BETA <- outgwis$BETA/2
outgwis$SE <- outgwis$SE/2

expgwas.LD.cmd <- paste(plink_bin, " --bfile ", stringname3, 
                        " --clump-p1 0.001 --clump-r2 0.3 --clump-kb 1024 --maf 0.05 --clump ", 
                        file.path(data_dir, paste0(exp,".GWAS.txt")), 
                        " --clump-snp-field SNP --clump-field P --out ", 
                        file.path(temp_dir, paste0("Tes3.", exp,".",phe)), sep="")
system(expgwas.LD.cmd)
expgwas.LD <- fread(paste0(file.path(temp_dir, paste0("Tes3.",exp,".",phe)), ".clumped"))
snp.expgwas.LD <- expgwas.LD$SNP
length(snp.expgwas.LD)
file.remove(paste0(file.path(temp_dir, paste0("Tes3.",exp,".",phe)), ".clumped"))
file.remove(paste0(file.path(temp_dir, paste0("Tes3.",exp,".",phe)), ".log"))

expgwis.snpcut <- expgwis[expgwis$P <= 1*10^(-3),"SNP"]

intersect.snp3 <- intersect(snp.expgwas.LD, expgwis.snpcut) 
intersect.snp3 <- intersect(intersect.snp3, outgwas$SNP) 
intersect.snp3 <- intersect(intersect.snp3, outgwis$SNP) 
print(paste0(exp,"-",phe," cutoff3 causal snps:",length(intersect.snp3)))

avbIndex <- match(intersect.snp3, bim$V2)
avbIndex <- as.matrix(avbIndex[order(avbIndex)])
intersect.snp3 <- bim[avbIndex, ]$V2

expgwas.order <- expgwas[match(intersect.snp3, expgwas$SNP), ]
expgwis.order <- expgwis[match(intersect.snp3, expgwis$SNP), ]
outgwas.order <- outgwas[match(intersect.snp3, outgwas$SNP), ]
outgwis.order <- outgwis[match(intersect.snp3, outgwis$SNP), ]

bh14ld <- expgwas.order$BETA
s124ld <- expgwas.order$SE
bh24ld <- outgwas.order$BETA
s224ld <- outgwas.order$SE
bh14ld_2 <- expgwis.order$BETA 
s124ld_2 <- expgwis.order$SE
bh24ld_2 <- outgwis.order$BETA
s224ld_2 <- outgwis.order$SE

lam <- 0.1; coreNum <- 1;
bp <- expgwas.order$BP
chr <- expgwas.order$CHR
idx4panel <- matrix(numeric(0), nrow = 0, ncol = 1)

Rblockres <- Cal_block_Rmatrix(bp, chr, avbIndex-1, idx4panel, block_file, stringname3, coreNum, lam);
R <- Rblockres$R; diag(R) <- 1 

save(bh14ld, bh24ld, bh14ld_2, bh24ld_2,
     s124ld, s224ld, s124ld_2, s224ld_2, R, intersect.snp3,
     file = file.path(beta_se_dir, paste0("inputGEI.",exp,".",phe,".RData")))

rho1 <- 0; rho2 <- 0
res2 <- MRGEI_Gam3seo(bh14ld, bh14ld_2, bh24ld, bh24ld_2, s124ld, s124ld_2, s224ld, s224ld_2, R, rho1, rho2)
str(res2)

save(res2, file = file.path(gam3seo_res_dir, paste0("res.",exp,".",phe,".RData")))

# write results
results <- c(exp, phe, 
             length(intersect.snp1),
             res.ivw$b, res.ivw$se, res.ivw$pval, 
             res.egger$Estimate, res.egger$StdError.Est, res.egger$Pvalue.Est,
             mean(res.raps$beta.hat), res.raps$beta.se, res.raps$beta.p.value,
             length(intersect.snp2),
             beta_hatLDP, se_hatLDP, pvalLDP,
             length(intersect.snp3),
             mean(res2$Beta1res), res2$Beta1.se, res2$Beta1.pval,
             mean(res2$Beta4res), res2$Beta4.se, res2$Beta4.pval)
results.df <- data.frame(t(results))
results.df

write.csv(results.df,file.path(output_dir, paste0(exp,".",phe,".results.csv")), row.names = FALSE)

}, error = function(e) {
  cat("Runtime error in this round：", e$message, "\n")
})
  
}

# set parallel work
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

batch_size <- 1

start_index <- (task_id - 1) * batch_size + 1
end_index <- min(task_id * batch_size, length(phe.list))

current_batch <- phe.list[start_index:end_index]

mclapply(current_batch, process_data, mc.cores = 1)

cat("Batch", task_id, "processed\n")


