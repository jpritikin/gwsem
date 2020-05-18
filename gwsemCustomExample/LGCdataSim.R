# This is a script to simulate data for the LGC GWAS model
# Load 
require(gwsem)

# Specify the basic model parameters

manifests<-c("time1","time2","time3","time4","snp")
latents<-c("intercept","slope")
path <- list(mxPath(from="intercept",to=c("time1","time2","time3","time4"), free=c(FALSE,FALSE,FALSE,FALSE), value=c(1.0,1.0,1.0,1.0) , arrows=1, label=c("intercept__time1","intercept__time2","intercept__time3","intercept__time4") ),
             mxPath(from="slope",to=c("time2","time3","time4"), free=c(FALSE,FALSE,FALSE), value=c(1.0,2.0,3.0) , arrows=1, label=c("slope__time2","slope__time3","slope__time4") ),
             mxPath(from="one",to=c("intercept","slope"), free=c(TRUE,TRUE), value=c(1.0,1.0) , arrows=1, label=c("const__intercept","const__slope") ),
             mxPath(from="snp",to=c("intercept","slope"), free=c(TRUE,TRUE), value=c(-.5,0.1) , arrows=1, label=c("snp__intercept","snp__slope") ),
             mxPath(from="intercept",to=c("intercept","slope"), free=c(TRUE,TRUE), value=c(1.0,0.1) , arrows=2, label=c("sigma_i","COV_intercept_slope") ),
             mxPath(from="slope",to=c("slope"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("sigma_s") ),
             mxPath(from="time1",to=c("time1"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="time2",to=c("time2"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="time3",to=c("time3"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="time4",to=c("time4"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="snp",to=c("snp"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_snp") ),
             mxPath(from="one",to=c("time1","time2","time3","time4","snp"), free=F, value=0, arrows=1))


# Begin the data simulation process 
model <- mxModel("lgc",  type="RAM", manifestVars = manifests, latentVars = latents, path)        # Build a shell of a model
sim1 <- mxGenerateData(model, 6000)                                                               # Simulate data for the basic model

for (ii in 1:5) {                                                                                 # Simulate 5 pcs
  sim1[[paste0('pc', ii)]] <- rnorm(6000)
}

head(sim1)

# Select a set of SNPs to make the GWAS more interesting

m1 <- buildItem(sim1, depVar = "time1")          # We are going to use this to to grab snps from the genetic file in their proper format.

# Choose K SNPs that we will make associated with the latent variable
set.seed(54321)
k <- 5                                                  # We are going to simulate 4 hits: 2 with the intecept and 2 with the slope
hits <- list()                                          # Make an object to put the "hits" in
hits[['lgm']] <- sample.int(2000, k)                    # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp1 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['lgm']][1])
snp2 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['lgm']][2])
snp3 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['lgm']][3])
snp4 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['lgm']][4])
snp5 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['lgm']][5])


# Add the effect of the SNPs to the latent variables
sim1[,1:4] <- (sim1[,1:4] + .24 * snp1$data$observed$snp + .27 * snp2$data$observed$snp + .35 * snp5$data$observed$snp)
sim1[,1:4] <- sim1[,1:4] + as.matrix(snp3$data$observed$snp) %*% .34 %*% matrix(c(0,1,2,3),1,4) + 
                           as.matrix(snp4$data$observed$snp) %*%  -.36 %*% matrix(c(0,1,2,3),1,4) +
                           as.matrix(snp5$data$observed$snp) %*%  -.42 %*% matrix(c(0,1,2,3),1,4)

write.table(sim1, "lgcData.txt", quote = F, row.names = F)



