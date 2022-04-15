# run & process sim

oldw <- getOption("warn")
options(warn = -1)
library(optparse) # for the cmdline args
library(seqgendiff)
library(sva)
library(tictoc) # for timing
options(warn = oldw)

# cmd line arguments are only passed to script if invoked by the Rscript cmd
option_list = list(make_option(c("--n_obs"),type="integer",default=10,
                               help="Number of Patients to sample (0 < && <= 98)",
                               metavar="integer"),
                   make_option(c("--n_resp"),type="integer",default=5,
                               help="number of responders (0 <= && <= n_obs)",
                               metavar = "integer"),
                   make_option(c("--seed"),type="integer",default=1234,
                               help="seed for the simulation",
                               metavar = "integer"),
                   make_option(c("--eff_sz_var"),type="double",default=0.8,
                               help="variance for effect size distribution",
                               metavar = "double"),
                   make_option(c("--prop_deg"),type="double",default=0.05,
                               help="proportion differentially expressed genes",
                               metavar = "double"),
                   make_option(c("--rda_file_header"),type="character",default="../P4_big_stuffs/data/test/",
                               help="file location relative to where the script is executing from",
                               metavar="character"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

n_obs = opt$n_obs
n_resp = opt$n_resp
seed = opt$seed
eff_sz_var = opt$eff_sz_var
prop_diff_exp_genes = opt$prop_deg

replace_obs = TRUE
Lib_size = 4e8
mean_lib_size = 5e7 # mean(kidney_sample_data$recount_qc.star.all_mapped_reads)
sim_file_loc = sprintf("%s",opt$rda_file_header)
sim_file_out = sprintf("%s/seed_%d.rda",sim_file_loc,seed)

load(file="./data/kidney_sample_data.rda") # kidney_sample_data
load(file="./data/X_ref.rda") # X_ref

set.seed(seed)


# Sample the patients
ix <- sample(ncol(X_ref), size = n_obs, replace = replace_obs)
pat.id <- kidney_sample_data$external_id[ix]
X_samp <- X_ref[,pat.id]

# Binomial thinning to get to target lib size
p_adj = Lib_size/(n_obs*mean_lib_size)
X_thinned <- thin_lib(X_samp, thinlog2 = rep(log(1/p_adj, base = 2), n_obs))

# Add the gene_exp_differences
X_np <- thin_2group(mat = X_thinned$mat,
                    prop_null = 1-prop_diff_exp_genes,
                    signal_fun = rnorm,
                    signal_params = list(mean =0, sd = eff_sz_var),
                    group_prop = n_resp/n_obs)


run_lm <- function(Y, D, Z = NULL, sandwich = FALSE, fixL_coef = FALSE){
  # Minimal argument checking
  stopifnot(is.matrix(Y))
  stopifnot(is.vector(D))
  stopifnot(length(D) == ncol(Y))
  
  L <- colSums(Y) 
  Y <- log2(Y + 0.5)
  n <- ncol(Y)
  
  # Set up the model matrix
  if(!fixL_coef){
    X <- cbind(rep(1, n), D, log2(L))
  }else{
    Y <- t(t(Y) - log2(L))
    X <- cbind(rep(1, n), D)
  }
  if(!is.null(Z)) X <- cbind(X, Z)
  
  # Compute the regression coefficients
  #XtXinv <- solve(t(X)%*%X)
  XtXinv <- solve(crossprod(X))
  #Px <- XtXinv%*%t(X)
  Px <- tcrossprod(XtXinv, X)
  #bhat <- Px%*%t(Y)
  bhat <- tcrossprod(Px, Y)
  #Ytilde <- t(bhat) %*% t(X)
  Ytilde <- tcrossprod(t(bhat), X)
  p <- ncol(X)
  l2fc <- bhat[2,]
  
  # Compute standard errors
  if(!sandwich){
    s2 <- rowSums((Y-Ytilde)^2)/(n-p)
    s <- sqrt(s2*XtXinv[2,2])
  }else{
    H <- X %*% Px
    Yresid2 <- (Y - Ytilde)^2
    omega <- t(t(Yresid2)/((1-diag(H))^2))
    s <- sapply(seq(nrow(Y)), function(i){
      O <- diag(omega[i,])
      V <- tcrossprod(crossprod(t(Px), O), Px)
      return(sqrt(V[2,2]))
    })
  }
  # Compute p-values
  p <- 2*pt(-abs(l2fc/s), df = n - p)
  res <- data.frame(l2fc=l2fc, se = s, p = p)
  return(res)
}

X <- cbind(X_np$design_obs, X_np$designmat, 
           log2(colSums(X_np$mat)))
Y <- log2(X_np$mat + 0.5)
sva_np <- sva(log2(X_np$mat + 0.5), mod =X, mod0 = X[,-2])

f2 <- run_lm(Y = X_np$mat, D = X_np$designmat[,1], Z = sva_np$sv)
truth = X_np$coefmat[,1]

dir.create(file.path(getwd(),sim_file_loc), showWarnings = FALSE,recursive = TRUE)
print(sim_file_out)
save(truth,f2,X_np,pat.id,file=sim_file_out)

BH_adj_pval <- p.adjust(f2$p, method = "BH")
ndeg <- sum(truth != 0)
genes_disc_05 <- sum(BH_adj_pval < 0.05)
tpr_05 <- sum(BH_adj_pval < 0.05 & truth != 0)/sum(truth != 0)
genes_disc_10 <- sum(BH_adj_pval < 0.10)
tpr_10 <- sum(BH_adj_pval < 0.10 & truth != 0)/sum(truth != 0)
fdr_05 <- 0
fdr_10 <- 0
if(genes_disc_05 > 0){
  fdr_05 <- sum(BH_adj_pval < 0.05 & truth == 0)/genes_disc_05
}
if(genes_disc_10 > 0){
  fdr_10 <- sum(BH_adj_pval < 0.10 & truth == 0)/genes_disc_10
}
lib_size <- colSums(X_np$mat)
avg_lib_size <- mean(lib_size)


report <- sprintf("NOBS: %d\nNRESP: %d\nseed: %s\nNDEG: %d\nGDISC05: %d\nTPR05: %f\nFDR05: %f\nGDISC10: %d\nTPR10: %f\nFDR10: %f\nAVGLIBSIZE: %f\nEFFVAR: %f\nPROPDEG: %f",n_obs,n_resp,seed,ndeg,genes_disc_05,tpr_05,fdr_05,genes_disc_10,tpr_10,fdr_10,avg_lib_size,eff_sz_var,prop_diff_exp_genes)
cat(report)
