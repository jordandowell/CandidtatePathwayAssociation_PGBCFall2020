# introduction to candidate pathway association 
#install necessary packages

library(qqman)



#ensure what directory you are in 
getwd()
# starting off recoding  our tped, & tfam files into bed files via plink

#give your computer permission to run external programs
# this may require you to do some googling on your part if this does not work.
#kewords: command line grant access folder binary
system("chmod -R 755 ../CandidtatePathwayAssociation_PGBCFall2020") # grant permissions (to avoid access denied errors when tyring to run the external binary)



#change NEWSNPs to whatever file you are interested 
#this needs to be the precursor characters for the .tped and .tfam files
getwd()
SNPset<-"Global_SNPlist"
#name the phenotype you are using 
Phenotype<-"CarotenoidContent"
#the following is using the plink software to recode the files into more usable forms
system(paste("Software/plink --tfile",paste0("data/", SNPset),"--covar",paste0("data/", SNPset,"_COVARIATES.txt") , "--make-bed", "--allow-no-sex --allow-extra-chr", "--out ",paste0("data/", SNPset)))





#create directories to store your output
dir.create(paste0(SNPset,"_Plots/"))
dir.create(paste0(SNPset,"_Plots/LMM/"))
dir.create(paste0(SNPset,"_Plots/BSLMM/"))
dir.create(paste0(SNPset,"_Tables/"))
dir.create(paste0(SNPset,"_Tables/LMM/"))
dir.create(paste0(SNPset,"_Tables/BSLMM/"))


#using GEMMA to run a mixed linear model 
system(paste("Software/gemma -bfile ",paste0("data/",SNPset)," -k ",paste0("data/",SNPset,".cXX.txt"), "-c ",paste0("data/",SNPset,".PCA_EV"), "-lmm 1 -outdir ",paste0(SNPset,"_Tables/LMM/"), "-o ",paste0(SNPset,"_",Phenotype,"_LMM")))


#Using GEMMA to run a Bayesian sparse linear mixed model( this make take like 10-15 minutes) 
system(paste("Software/gemma -bslmm 1 -bfile ",paste0("data/",SNPset)," -w 250000"," -rpace 100 -wpace 1000"," -k ",paste0("data/",SNPset,".cXX.txt")," -outdir", paste0(SNPset,"_Tables/BSLMM/"), "-o " ,paste0(SNPset,"_",Phenotype,"_BSLMM")))

# both analyses are run so lets compare results 



# lets read in the LMM outputs to make a manhattan plot 
snips<-read.table(paste0(SNPset,"_Tables/LMM/",SNPset,"_",Phenotype,"_LMM",".assoc.txt"),sep="",header = T)


snp.map<-read.table(paste0("data/",SNPset,".map"), sep="",header=F, col.names = c("chr","rs","nothing","ps"))

head(snp.map)

#identify which column has the pvalue

pvalue<-snips$p_wald

ASSOC_PVALUES<- merge(snp.map,snips,by=c("chr","rs","ps"))
#we need to convert the chr to numbers 

ASSOC_PVALUES$chr_num<-as.numeric(as.factor(ASSOC_PVALUES$chr))

#lets make a qq plot....notice there are fewer high pvalue that expected.
qq(ASSOC_PVALUES$p_wald, main = "Q-Q plot of GWAS p-values")

# manhattan plot can tell up similar information here the pvalue are so far below the genome wide line(bonferroni correction )
manhattan(ASSOC_PVALUES,chr = "chr_num", bp = "ps", p = "p_wald", snp = "rs",logp = T,
          suggestiveline = -log10(0.05),genomewideline=-log10(0.05/(687)), chrlabs = as.character(unique(ASSOC_PVALUES$chr_num)))



# lets check we may have a strong polygenic architecture on our hands 
#modified from Victor Soria-Carrasco


# Load hyperparameter file for p-value
# ==============================================================================
hyp.params<-read.table(paste0(SNPset,"_Tables/BSLMM/",SNPset,"_",Phenotype,"_BSLMM",".hyp.txt"),header = T)

# ==============================================================================
# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))
# rho -> approximation to proportion of genetic variance explained by variants with major effects,
#rho close to zero trait is highly polygenic
rho<-c("rho",mean(hyp.params$rho),quantile(hyp.params$rho, probs=c(0.5,0.025,0.975)))


# ==============================================================================

# get table of hyperparameters
# ==============================================================================
hyp.params.table<-as.data.frame(rbind(pve,pge,pi,n.gamma,rho),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
View(hyp.params.table)
write.table(hyp.params.table,paste0(paste0(SNPset,"_Tables/BSLMM/",paste0(SNPset,"_",Phenotype,".Hyperparametertable.txt"))))
#

# plot traces and distributions of hyperparameters
# ==============================================================================
#pdf(file=paste0(paste0(SNPset,"_Plots/BSLMM/",paste0(SNPset,"_",Phenotype,".Hyperparameters.pdf"))), width=8.3,height=11.7)
#layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi 
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# N_gamma #number of large effect snps
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------

# rho
# ------------------------------------------------------------------------------
plot(hyp.params$rho, type="l", ylab="rho", main="rho - trace")
hist(hyp.params$rho, main="rho - posterior distribution", xlab="rho")
plot(density(hyp.params$rho), main="rho - posterior distribution", xlab="rho")
# ------------------------------------------------------------------------------
#h
# ------------------------------------------------------------------------------
plot(hyp.params$h, type="l", ylab="h", main="h - trace")
hist(hyp.params$h, main="h - posterior distribution", xlab="h")
plot(density(hyp.params$h), main="h - posterior distribution", xlab="h")





#dev.off()



# Load parameters output
# ==============================================================================
params<-read.table(paste0(SNPset,"_Tables/BSLMM/",SNPset,"_",Phenotype,"_BSLMM.param.txt"),header=T,sep="\t")
# ==============================================================================
# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta(effect all variants have) * gamma(effect of the SNP)) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
View(head(params.effects.sort, 10) )

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file=paste0(paste0(SNPset,"_Tables/BSLMM/",paste0(SNPset,"_",Phenotype,".top1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top01, file=paste0(paste0(SNPset,"_Tables/BSLMM/",paste0(SNPset,"_",Phenotype,".top0.1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top001, file=paste0(paste0(SNPset,"_Tables/BSLMM/",paste0(SNPset,"_",Phenotype,".top0.01eff.dsv"))), quote=F, row.names=F, sep="\t")


# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip01.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip10, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip10.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip25, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip25.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip50, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip50.dsv"))), quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
# ==============================================================================

#add the numeric chromosome 
params$chr_num<-as.numeric(as.factor(params$chr))

# manhattan plot can tell up similar information here the pvalue are so far below the genome wide line(bonferroni correction )
manhattan(params,chr = "chr_num", bp = "ps", p = "gamma", snp = "rs",logp = F,
          suggestiveline = F,genomewideline=F, chrlabs = as.character(unique(params$chr_num)))




