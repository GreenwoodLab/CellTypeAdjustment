#Running the different cell type adjustment methods on the simulated beta value matrix
# sim_beta is the matrix of simulated beta values (probes x samples).
# disease_status is the phenotype (0=control, 1=case)
# dmr contains the indices of CpGs that are differentially methylated with the phenotype

load("simulated_data.RData")

#Random sample of data to speed up demonstration, if desired
set.seed(24601)
nprobes = 50000
rs = sample(1:dim(sim_beta)[1], nprobes)
sim_beta = sim_beta[rs,]
dmr = dmr[dmr %in% rs]

###################################################
#Model not adjusting for cell type
library(limma)
X.unadj = model.matrix(~disease_status)
lm.unadj = eBayes(lmFit(sim_beta, X.unadj))
p.unadj = lm.unadj$p.value[,2]

###################################################
#Code for reference-based and reference-free methods based on Houseman tutorial 2014
#http://people.oregonstate.edu/~housemae/software/TutorialLondon2014/
###################################################
#Reference-based method
library(minfi)
library(FlowSorted.Blood.450k)
data(FlowSorted.Blood.450k.JaffeModelPars)
commondmr <- intersect(rownames(FlowSorted.Blood.450k.JaffeModelPars), rownames(sim_beta))

#Directly estimating cell type composition
cellcomps = minfi:::projectCellType(sim_beta[commondmr,], FlowSorted.Blood.450k.JaffeModelPars[commondmr,], lessThanOne=TRUE)

X.ref = model.matrix(~disease_status+cellcomps)
lm.ref = eBayes(lmFit(sim_beta, X.ref))
p.ref = lm.ref$p.value[,2]

###################################################
#Reference-free method
library(RefFreeEWAS)

#Estimated latent dimension
est_dim = EstDimRMT(t(scale(t(sim_beta - lm.unadj$coef %*% t(X.unadj)))), plot=FALSE)$dim
#Main function
mod = RefFreeEwasModel(sim_beta, cbind(1,disease_status), K=est_dim)
nboot = 100 #Number of bootstrap samples to obtain standard errors
boot = BootRefFreeEwasModel(mod, nboot)
#Extracting standard errors and calculating test statistic
se = apply(boot[,,"B",], 1:2, sd)
ts = abs(mod$Beta)/se

df = -diff(dim(model.matrix(~disease_status)))-apply(is.na(sim_beta),1,sum)
p.reffree = 2*pt(ts[,2], df=df, lower.tail=FALSE)

###################################################
#Surrogate variable analysis
library(sva)

#Null model matrix must be nested in the full model matrix
model_mat = model.matrix(~disease_status)
null_model_mat = model.matrix(~1, data.frame(disease_status))

#Main SVA function
SVA = sva(sim_beta, model_mat, null_model_mat, method="irw")

#Calculate p-values
update_model = cbind(model_mat, SVA$sv)
update_null_model = cbind(null_model_mat, SVA$sv)
p.sva = f.pvalue(sim_beta, update_model, update_null_model)

###################################################
#Independent Surrogate Variable Analysis
library(isva)

ISVA = DoISVA(sim_beta, matrix(disease_status))
p.isva = ISVA$spv

###################################################
#Deconfounding
library(deconf)

#Running deconfounding on subset to speed up demonstration
sampdec = sample(length(sim_beta[,1]), 1000)
tmpdec = deconfounding(sim_beta[sampdec,],2)

X.tmpdec = model.matrix(~disease_status+t(tmpdec$C$Matrix))
lm.tmpdec = eBayes(lmFit(sim_beta, X.tmpdec))
p.dec = lm.tmpdec$p.value[,2]

##################################################
#RUV
library(ruv)

#Choosing Jaffe CpG sites to be control probes for ruv
ruv_controls = which(rownames(sim_beta) %in% rownames(FlowSorted.Blood.450k.JaffeModelPars))
#Taking only the controls that do not correlate with the phenotype
subset_controls = ruv_controls[lm.unadj$p.value[ruv_controls,2]>0.01]
#Turning control probe vector into logical
ctl = rep(FALSE, dim(sim_beta)[1]); ctl[subset_controls] = TRUE

#getK estimates the latent dimension, but for this dataset it ends up being much too high
ruvK = getK(t(sim_beta), matrix(disease_status), ctl)$k
ruv_mod = RUV4(t(sim_beta), matrix(disease_status), ctl, k=ruvK)
p.ruv = ruv_mod$p

##################################################
#EWASher and CellCDecon need to be run externally... See scripts for these two methods

#Function to save data to text files for use in EWASher and CellCDecon 
#Function from http://research.microsoft.com/en-us/downloads/472fe637-7cb9-47d4-a0df-37118760ccd1/
makefiles_FLE = function(rSrcDir, data, phen, covar)
{        
  #write each input to file, then call RunFLEFromFiles
  write.table(data, file = "datafileTmp.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(phen, file = "phenfileTmp.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(covar, file = "covarfileTmp.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

#Need to manipulate beta matrix to what EWASher and CCD expect
tmp_beta = data.frame(rownames(sim_beta), sim_beta)
colnames(tmp_beta) = c("ID", colnames(sim_beta))

#Directory to output files to
path_to_dir = "/home/data1/homeldi/kevin.mcgregor/research/marie_hudson/r_scripts/mixtures/github_tutorial"

#Make the EWASher input files
makefiles_FLE(path_to_dir, tmp_beta, data.frame(colnames(tmp_beta)[-1], disease_status), covar=NULL) 

#Or save another version with patients' true disease (coded) from original 450K data before simulation
#as a covariate to compare performance of EWASher
makefiles_FLE(path_to_dir, tmp_beta, 
              data.frame(colnames(tmp_beta)[-1], disease_status),
              data.frame(1, true_disease))  




