# divDyn routine
require(divDyn)
require(here)
require(dplyr)

# stages can be associated to the data -stages will be the basis for the plots

data("stages")
head(stages)
tsplot(stages, boxes="sys", shading= "series")
tsplot(stages, boxes="system", shading="stage", xlim=59:81)
tsplot(stages, boxes="sys", shading="series",
       labels.args=list(col="red", font=3), shading.col=c("white", "wheat"))
tsplot(stages, boxes=c("sys"), shading="sys", boxes.col="systemCol")


# load data
load (here ("processed_data", "table_data_array_cynodontia.RData"))
coll_occ_taxa_perm_cret$mid <- rowMeans (coll_occ_taxa_perm_cret [,c ("min_ma","max_ma")])

# stages in our data
stages <- stages [which(stages$stage %in% coll_occ_taxa_perm_cret$Interval),]

tsplot(stages, shading="stage", boxes="sys",xlim=c(270,62))
ranges(coll_occ_taxa_perm_cret, 
       tax="unique_name", bin="mid",
       labs=T,
       labels.args=list(cex=0.6),  occs=T)


#  survivorship curves. The function survivors() calculates the proportions of survivors
# from every bin to all the remaining bins
surv_forward <- survivors(coll_occ_taxa_perm_cret, 
                  bin="bin_assignment", 
                  tax="unique_name",
                  method = "forward")

# backward
surv_backward <- survivors(coll_occ_taxa_perm_cret, 
                          bin="bin_assignment", 
                          tax="unique_name",
                          method = "backward")

# it starts in 1; delete empty cols and rows
surv_backward<-(surv_backward [which(rowSums(surv_backward,na.rm=T)>0),
                              which(colSums(surv_backward,na.rm=T)>0)])
surv_forward<-(surv_forward [which(rowSums(surv_forward,na.rm=T)>0),
                               which(colSums(surv_forward,na.rm=T)>0)])


# time scale plot
tsplot(stages, shading="stage", boxes="sys",
       xlim=c(260,60), ylab="proportion of survivors (red) and constituint (black) present",
       ylim=c(0.01,1),plot.args=list(log="y"))

for(i in 1:ncol(surv_backward)) lines(stages$mid, surv_backward[,i])
for(i in 1:ncol(surv_forward)) lines(stages$mid, surv_forward[,i],col="red")

# basic dataset statistics
samp <-sumstat(coll_occ_taxa_perm_cret, 
               tax="unique_name", 
               bin="bin_assignment", duplicates=T)
samp
samp <-binstat(coll_occ_taxa_perm_cret, 
               tax="unique_name", 
               bin="bin_assignment",
               indices=T)
samp

# time scale plot
oldPar <- par(mar=c(4,4,2,4))
tsplot(stages, shading="stage", boxes="sys",boxes.col="systemCol",
       xlim=c(266,64), ylab="Number of genus",
       ylim=c(0.01,600))
lines(stages$mid, rev(samp$richness[is.na(samp$richness)!=T]))
lines(stages$mid, rev(samp$chao1occ[is.na(samp$chao1occ)!=T]),col= "red")
lines(stages$mid, rev(samp$chao2[is.na(samp$chao2)!=T]),col= "blue")


# the collections (rescaled, other axis)
lines(stages$mid, rev(samp$occs[is.na(samp$occs)!=T]), col="green4",lwd=2)
axis(4, col="green4",
     col.ticks="green4",
     col.axis="green4",
     at=seq(0,600,100), labels=seq(0,600,100))
mtext(4, text="number of occurrences (records)", col="green4", line=2)
par(oldPar)

#plotting
tsplot(stages, shading="series", boxes="sys", #xlim=52:95,
       ylab="number of occurences", ylim=c(0,1000))
parts(coll_occ_taxa_perm_cret$mid, 
      coll_occ_taxa_perm_cret$lithology1.x) # must be character


# diversity dynamics
ddFirst<-divDyn(coll_occ_taxa_perm_cret, 
                bin="bin_assignment", 
                tax="unique_name", 
                noNAStart=F)

# plot
tsplot(stages, shading="stage", boxes="sys",boxes.col="systemCol",
       xlim=c(266,64), ylab="Number of genus",
       ylim=c(0.01,100))

# pull of the recent pattern 

lines(stages$mid, rev(ddFirst$tThrough),col="red",lwd=2)
lines(stages$mid, rev(ddFirst$divRT),col="red",lwd=2)
lines(stages$mid, rev(ddFirst$divBC),col="blue",lwd=2)
lines(stages$mid, rev(ddFirst$divSIB),col="green4",lwd=2)
lines(stages$mid, rev(ddFirst$divCSIB),col="orange",lwd=2)

# legend
legend("topleft", legend=c("RT", "BC", "SIB","CSIB"),
       col=c("red", "blue", "green","orange"), lwd=c(2,2,2,2), bg="white")

# origination and extinciton


# plot
tsplot(stages, shading="stage", boxes="sys",boxes.col="systemCol",
       xlim=c(266,64), ylab="Number of genus",
       ylim=c(0.01,3))

# proportional rates
lines(stages$mid, ddFirst$extPC,col="red",lwd=2)
lines(stages$mid, ddFirst$oriPC,col="blue",lwd=2)

# plot
tsplot(stages, shading="stage", boxes="sys",boxes.col="systemCol",
       xlim=c(266,64), ylab="Number of genus",
       ylim=c(0.01,5))

# per capita rates
lines(stages$mid, ddFirst$extPC,col="black",lwd=2)
#lines(stages$mid, ddFirst$oriPC[-c(1:6)],col="blue",lwd=2)
# gap filler
lines(stages$mid, ddFirst$extGF,col="blue",lwd=2)
#lines(stages$mid, ddFirst$oriGF[-c(1:6)],col="blue",lwd=2)
# second-third
lines(stages$mid, ddFirst$ext2f3,col="green",lwd=2)
#lines(stages$mid, ddFirst$ori2f3[-c(1:6)],col="blue",lwd=2)
# legend
legend("topright", legend=c("per capita", "gap-filler", "second-for-third"),
       col=c("black", "blue", "green"), lwd=c(2,2,2), bg="white")


# continuous time
ddIDbin <- divDyn(coll_occ_taxa_perm_cret, 
                  bin="mid", 
                  tax="unique_name", revtime=TRUE)

# basic plot
tsplot(stages, shading="stage", boxes="sys", #xlim=52:95,
       ylab="Diversity, range-through", ylim=c(0,100))
lines(ddIDbin$mid, (ddIDbin$divRT), col="black", lwd=2)
lines(stages$mid, rev(ddFirst$divRT), col="red", lwd=2)
legend("topleft", legend=c("unique mid entries", "stg stages"),
       col=c("black", "red"), lwd=c(2,2), bg="white")


# subsampling

## subsampled, stable
withoutFail <- subsample(coll_occ_taxa_perm_cret, bin="bin_assignment", 
                     tax="unique_name",
                     #coll="collection_no",
                     iter=50, q=0.5, duplicates=T,
                     useFailed=F,
                     type="sqs",
                     output="dist")
withFail <- subsample(coll_occ_taxa_perm_cret, bin="bin_assignment", 
                         tax="unique_name",
                         #coll="collection_no",
                         iter=50, 
                      q=0.5, duplicates=T,
                         useFailed=T,
                      type="sqs",
                      output="dist")

# basic plot
tsplot(stages, shading="stage", boxes="sys", #xlim=52:95,
       ylab="number of occurrences", ylim=c(0,100))

# uncertainty
shades(stages$mid, withoutFail$divCSIB,  col="gray25",res=c(0.05,0.25,0.75,0.95))
shades(stages$mid, withFail$divCSIB, col="blue",res=c(0.05,0.25,0.75,0.95))

# subsampled
lines(stages$mid, 
      rowMeans (withoutFail$divCSIB,na.rm=T), col="black", lwd=2)

# subsampled, with failed
lines(stages$mid, 
      rowMeans (withFail$divCSIB,na.rm=T), col="cyan", lwd=2)

legend("topleft", legend=c("SQS, q=0.5, without failed",
                           "SQS, q=0.5, with failed"), col=c("black", "blue"),
       lwd=c(2,2), bg="white")

# environmental affinity
coll_occ_taxa_perm_cret$temperature_bin <- ifelse (coll_occ_taxa_perm_cret$temperature < 15,"cold", "warm")

# majority rule method
affMajor <- affinity(coll_occ_taxa_perm_cret,
                     bin="bin_assignment", 
                     tax="unique_name",
                     method="majority", 
                     #coll="collection_no", 
                     env="temperature_bin")
table(affMajor)

# binomial
affBin1 <- affinity(coll_occ_taxa_perm_cret,
                    bin="bin_assignment", 
                    tax="unique_name",
                    method="binom", 
                    #coll="collection_no", 
                    env="temperature_bin", alpha=0.1)
table(affBin1)

