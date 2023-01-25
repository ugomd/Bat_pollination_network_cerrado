####################################################################################
#                         Diniz and Aguiar - Data for:                             #
# Beyond a 'perfect fit: the interplay between morphology and spatiotemporal over- #
#                     lap in a bat-flower network                                  #
#             Scientific Reports, Under review (v. 18.01.2022)                     #
#                                                                                  #     
# Authors: Ugo M. Diniz (Technische Universität München); Ludmilla M. S. Aguiar    #
# (Universidade de Brasília)                                                       #
#                                                                                  #
####################################################################################

####################################################################################
########################## Preparation: Packages ###################################
####################################################################################

### Installing/loading required packages
install.packages("bipartite")
install.packages("ggplot2")
install.packages("gdata")
install.packages("vegan")

require(bipartite)
require(gdata)
require(ggplot2)
require(vegan)


####################################################################################
###################### 1. Sampling Completeness analysis ###########################
####################################################################################

# Loading vectors conataning abundance data of bats and flowering plants (all plant species, and chiropterophilous only)
 
batspp<- c(24,5,1,79,24,62,17,36,43,22,4,37,13,1)
plantspp <- c(46,38,27,28,25,23,22,13,11,10,8,7,5,5,5,4,4,3,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
plantsspp_c <- c(46,38,27,28,25,22,13,11,10,8,7,5,4,4,3)

# Estimating maximum richness (Chao1)

estimateR(batspp, index =c("chao1"))    #bats
estimateR(plantspp, index =c("chao1"))   #all plants
estimateR(plantsspp_c, index =c("chao1"))  #chiropterophilous plants

# Rarefaction curves (loading matrices with bat/plant individuals found per sampling event and plotting curves)

sampbats <- read.xls("sampbats.xlsx", h=T)
sampplantsfull <- read.xls("sampplantsfull.xlsx", h=T)
sampplantschiro <- read.xls("sampplantschiro.xlsx",h=T)

# Bats
tiff("sampbats.tiff", width = 10, height = 15, units = "cm", res = 600)  
curve_bat<- specaccum(sampbats, method="rarefaction")
plot(curve_bat, xvar="individuals", ci.type = "poly", 
     ci.col="grey80", ci.lty = 0, ylim=c(0,18)) +
  abline(h=15, lty=1, col="black", lwd=1.5)
abline(h=(15+2.29), lty=3, col="grey50", lwd=2) # estimate + SD (Chao1)
abline(h=(15-2.29), lty=3, col="grey50", lwd=2)
dev.off()

# Plants (all)
tiff("sampplantsfull.tiff", width = 10, height = 15, units = "cm", res = 600)  
curve_plants_full<- specaccum(sampplantsfull, method="rarefaction")
plot(curve_plants_full, xvar="individuals", ci.type = "poly", 
     ci.col=rgb(0.7, 0.4, 0.7, 0.3), ci.lty = 0, ylim=c(0,100),xlim=c(0,300)) +
  abline(h=70, lty=1, col="black", lwd= 2)
abline(h=(70+25.59), lty=3, col="grey50", lwd=2) # estimate + SD (Chao1)
abline(h=(70-25.59), lty=3, col="grey50", lwd=2)
dev.off()

# Plants (Chiropterophilous)
tiff("sampplantschiro.tiff", width = 10, height = 15, units = "cm", res = 600)  
curve_plants_full<- specaccum(sampplantschiro, method="rarefaction")
plot(curve_plants_full, xvar="individuals", ci.type = "poly", 
     ci.col=rgb(1, 0, 0, 0.2),   ci.lty = 0, ylim=c(5,16)) +
  abline(h=15, lty=1, col="black", lwd=2)
dev.off()


####################################################################################
############## 2. Temporal and spatial patterns of bat species #####################
####################################################################################

# Loading temporal data of species' relative abundances (months are letter coded, A = October, B = November...)
temporal<-read.xls("temporal.xlsx",h=T)

tempnectar<-temporal[temporal$Guild=="nect",] # Selecting nectarivores only

tiff("tempnectar.tiff", width = 20, height = 12, units = "cm", res = 600)  
ggplot(data=tempnectar, aes(x=Month, y = Frequency, group=Species)) +  
  geom_line(aes(color = Species), lwd=.8) + geom_point(aes(color=Species, shape=Species), cex=3) + theme_classic()+
  scale_color_brewer(palette="Accent") + theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
                                               axis.text= element_text(color="black", size=15),
                                               legend.title = element_text(color="black", size=15))
dev.off()


tempfrug<- temporal[temporal$Guild=="frug",] # Selecting frugivores only

tiff("tempfrug.tiff", width = 20, height = 12, units = "cm", res = 600) 
ggplot(data=tempfrug, aes(x=Month, y = Frequency, group=Species)) +  
  geom_line(aes(color = Species), lwd=.8) + geom_point(aes(color=Species, shape=Species), cex=3) + theme_classic()+
  scale_color_brewer(palette="Paired") + theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
                                               axis.text= element_text(color="black", size=15),
                                               legend.title = element_text(color="black", size=15))
dev.off()


# Loading spatial data (habitats are letter-coded, A = savanna, B = edge, C = forest)

spatial <-read.xls("spatial.xlsx",h=T)

spatnectar<-spatial[spatial$Guild=="nect",] # Selecting nectarivores only

tiff("spatnectar.tiff", width = 20, height = 12, units = "cm", res = 600)  
ggplot(data=spatnectar, aes(x=Habitat, y = Frequency, group=Species)) +  
  geom_line(aes(color = Species), lwd=.8) + geom_point(aes(color=Species, shape=Species), cex=3) + theme_classic()+
  scale_color_brewer(palette="Accent") + theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
                                               axis.text= element_text(color="black", size=15),
                                               legend.title = element_text(color="black", size=15))
dev.off()


spatfrug <-spatial[spatial$Guild=="frug",] # Selecting frugivores only

tiff("spatfrug.tiff", width = 20, height = 12, units = "cm", res = 600)  
ggplot(data=spatfrug, aes(x=Habitat, y = Frequency, group=Species)) +  
  geom_line(aes(color = Species), lwd=.8) + geom_point(aes(color=Species, shape=Species), cex=3) + theme_classic()+
  scale_color_brewer(palette="Paired") + theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
                                               axis.text= element_text(color="black", size=15),
                                               legend.title = element_text(color="black", size=15))
dev.off()


####################################################################################
###################### 3. Calculating network-level indices ########################
####################################################################################


# Loading the interaction matrix and visualizing the network
network <- as.matrix(read.delim("network.txt", row.names=1))
dim(network)
plotweb(network)
visweb (network)

# Null model and network macrostructure (Qw, H2', and WNODA) 
null <- nullmodel(network, N=1000, method="vaznull")

Qw <- computeModules(network, method = "Beckett") 
Qw@likelihood # Modularity
plotModuleWeb(Qw, labsize=0.5) # module visualization

nullmod <- sapply(null, computeModules)
modunull <- sapply(nullmod, function(x) x@likelihood)
(modu@likelihood - mean(modunull))/sd(modunull) # F-value modularity
sum(modunull>(modu@likelihood)) / length(modunull) # p-value modularity


H2 <- networklevel(network, index="H2") 
H2 # Specialization

specnull <- unlist(sapply(null, networklevel, index="H2"))
(spec-mean(specnull))/sd(specnull) # F-value H2'
sum(specnull>nest)/length(specnull) # p-value H2'


# WNODA and compound topology

######################################################################################################
# OBS: WNDOA calculation and the compound topology analysis were taken from Pinheiro et al. (2022),  #  
# who kindly made their code public. The functions, the script and observations have been made       #
# entirely by them and adapted to our data. Please refer to the original publication for a thorough  #
# and detailed explanation on the functions.                                                         #
#                                                                                                    #
# Pinheiro, R. B., Felix, G. M., & Lewinsohn, T. M. (2022). Hierarchical compound topology           #
# uncovers complex structure of species interaction networks.                                        #
# Journal of Animal Ecology, 91(11), 2248-2260.                                                      #
#                                                                                                    #
# Data and codes available :  https://doi.org/10.5281/zenodo.7007952                                 #
######################################################################################################
  
# required packages
install.packages("rstudioapi")
library(rstudioapi)
source("RestNullModel.R")
source("PosteriorProb.R")
source("MyDiamond.R")

Part <- bipartite::module2constraints(Qw)
row.Part <- Part[1:nrow(network)]
col.Part <- Part[(nrow(network)+1):(nrow(network)+ncol(network))]
WNDOA <- unlist(bipartite::nest.smdm(x = network, 
                                       constraints = Part, 
                                       weighted = T, 
                                       decreasing = "abund"))


#### Functions ###
PosteriorProb <- function(M, R.partitions, C.partitions, Prior.Pij, Conditional.level){
  
  # Test of assumptions
  if (!is.matrix(M)){stop("M is not a matrix")}
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {stop("M is degenerated. There are rows and/or columns without interactions in the matrix. Remove them before proceding")}
  if (!is.numeric(R.partitions) | !is.numeric(C.partitions)) {stop("Partitions are not numeric")}
  if (length(R.partitions) != nrow(M) | length(C.partitions) != ncol(M)) {stop("Partitions and matrix dimensions have different sizes")}
  if (!(Conditional.level %in% c("matrix","modules","areas"))) {stop("Conditional.level should be 'matrix','modules' or 'areas'")}
  if (Prior.Pij != "degreeprob" & Prior.Pij != "equiprobable" & Prior.Pij != "degreeprob.byarea") {stop("Pij.probs should be 'equiprobable' or 'degreeprob' or 'degreeprob.byarea")}
  
  # M dimensions
  r <- dim(M)[1] # Number of rows
  c <- dim(M)[2] # Number of columns
  array()
  
  # Making an array with r rows, c columns, and 3 slices. This array represents the modular structure.
  # The first slice informs if a given cell M(rc) is within (1) or outside (0) a module.
  # The second slice informs to which module the species in the row (r) of a given cell M(rc) belongs.
  # The third slice informs to which module the species in the column (c) of a given cell M(rc) belongs .
  
  Matrix.mod <- array(0, dim = c(r, c, 3))
  for (rr in 1:r){
    for (cc in 1:c){
      Matrix.mod[rr,cc,1] <- ifelse(R.partitions[rr] == C.partitions[cc], 1,0)
      Matrix.mod[rr,cc,2] <- R.partitions[rr]
      Matrix.mod[rr,cc,3] <- C.partitions[cc]
    }
  }
  
  # Defining a priori Pij probabilities.
  if (Prior.Pij == "equiprobable"){
    Pi <- rep(1 / r, times = r) 
    Pj <- rep(1 / c, times = c)
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }else if (Prior.Pij == "degreeprob"){
    Pi <- rowSums(M) / sum(rowSums(M)) 
    Pj <- colSums(M) / sum(colSums(M))
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }else if(Prior.Pij == "degreeprob.byarea"){
    Prior.Pij.species <- M
    RMod <- sort(unique(R.partitions))
    CMod <- sort(unique(C.partitions))
    for (rr in RMod){
      for (cc in CMod){
        M.rr.cc <- matrix(M[R.partitions == rr,C.partitions == cc], sum(1*(R.partitions == rr)), sum(1*(C.partitions == cc)))
        Pi.rr.cc <- rowSums(M.rr.cc) / sum(rowSums(M.rr.cc)) 
        Pj.rr.cc <- colSums(M.rr.cc) / sum(colSums(M.rr.cc))
        Prior.Pij.species[R.partitions == rr, C.partitions == cc] <- tcrossprod(Pi.rr.cc, Pj.rr.cc)
      }
    }
  }
  
  # Defining conditional probabilities by area based on species degrees and connectance by area.
  
  if (Conditional.level == "matrix"){
    Post.Pij <- Prior.Pij.species
  }else {
    Prior.Pij.area <- matrix(NA,r,c)
    Cond.Pij.area <- matrix(NA,r,c)
    if (Conditional.level == "modules"){
      WMod.prior <- sum(Prior.Pij.species[Matrix.mod[,,1] == 1])
      OMod.prior <- sum(Prior.Pij.species[Matrix.mod[,,1] == 0])
      Prior.Pij.area[Matrix.mod[,,1] == 1] <- WMod.prior
      Prior.Pij.area[Matrix.mod[,,1] == 0] <- OMod.prior
      
      WMod.cond <- sum(M[Matrix.mod[,,1] == 1]) / sum(M)
      OMod.cond <- sum(M[Matrix.mod[,,1] == 0]) / sum(M)
      Cond.Pij.area[Matrix.mod[,,1] == 1] <- WMod.cond
      Cond.Pij.area[Matrix.mod[,,1] == 0] <- OMod.cond
    }else if (Conditional.level == "areas"){
      RMod <- sort(unique(R.partitions))
      CMod <- sort(unique(C.partitions))
      for (rr in RMod){
        for (cc in CMod){
          WArea.prior <- sum(Prior.Pij.species[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc])
          Prior.Pij.area[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc] <- WArea.prior
          
          WArea.cond <- sum(M[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc]) / sum(M)
          Cond.Pij.area[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc] <- WArea.cond
        }
      } 
    }
    
    # Adjusting the prior Pij prob by conditional probabilities. 
    
    Post.Pij <- Prior.Pij.species * (Cond.Pij.area / Prior.Pij.area)
  }
  
  return(Post.Pij = Post.Pij)
}



RestNullModel <- function(M, Pij.Prob, Numbernulls, Print.null = F, allow.degeneration = F, 
                          return.nonrm.species = T, connectance = T, byarea = F, R.partitions = NULL, C.partitions = NULL){
  
  ### Test of assumptions
  
  if (!is.matrix(M)){stop("M is not a matrix")}
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {stop("M is degenerated")}
  
  if (!is.matrix(Pij.Prob)){stop("Pij is not a matrix")}
  if (T %in% c(Pij.Prob < 0)){stop("Pij must contain only numbers >= 0")}
  
  if (nrow(M) != nrow(Pij.Prob) | ncol(M) != ncol(Pij.Prob)){stop("Dimensions of M and Pij.Prob should be identical")}
  
  if (byarea == T){
    if(is.null(C.partitions) | is.null(R.partitions)){stop("Partitions missing")}
    if (length(unique(c(length(R.partitions),nrow(M),nrow(Pij.Prob)))) != 1){stop("The number of elements of R.partition should be the same as the number of rows of M and Pij.prob")}
    if (length(unique(c(length(C.partitions),ncol(M),ncol(Pij.Prob)))) != 1){stop("The number of elements of C.partition should be the same as the number of column of M and Pij.prob")}
    if(!identical(sort(unique(R.partitions)), sort(unique(C.partitions)))){stop("The number and labels of modules in R.partition and C.partition must be the same")}
  }  
  
  if (Numbernulls <= 0 | !is.numeric(Numbernulls)) {stop("Numbernulls should be a number > 0")}
  if (!is.logical(connectance)){stop("connectance should be logical (T or F)")}
  if (!is.logical(allow.degeneration)){stop("allow.degeneration should be logical (T or F)")}
  if (!is.logical(return.nonrm.species)){stop("return.nonrm.species should be logical (T or F)")}
  if (!is.logical(byarea)){stop("byarea should be logical (T or F)")}
  
  ### M dimensions
  r <- dim(M)[1] # Number of rows
  c <- dim(M)[2] # Number of collums
  
  ### Constructing a array with r rows, c columns and 2 slices. This array represents the matrix area structure
  if (byarea == T){
    Matrix.area <- array(0, dim = c(r, c, 2))
    for (rr in 1:r){
      for (cc in 1:c){
        Matrix.area[rr,cc,1] <- R.partitions[rr]
        Matrix.area[rr,cc,2] <- C.partitions[cc]
      }
    }
    
  }else if (byarea == F){
    ## Assigning all rows and columns to the same partition in order to run the code bellow
    Matrix.area <- array(1, dim = c(r, c, 2))
    R.partitions <- rep(1, nrow(M))
    C.partitions <- rep(1, ncol(M))
    
  }
  
  ### Null model simulation
  
  NullMatrices <- list() # list where the null matrices will be storage
  length(NullMatrices) <- Numbernulls #assigning the number of null matrices to be saved in NullMatrices 
  
  ## Drawing interaction in each null matrix 
  for (nn in 1:Numbernulls){
    
    R.part <- sort(unique(as.vector(Matrix.area[,,1])))
    C.part <- sort(unique(as.vector(Matrix.area[,,2])))
    finalmat <- matrix(NA, r, c)
    
    for (R.p in R.part){
      for (C.p in C.part){
        
        M.a <- as.matrix(M[R.partitions == R.p, C.partitions == C.p])
        Pij.a <- Pij.Prob[R.partitions == R.p, C.partitions == C.p]
        
        r.a <- dim(M.a)[1]
        c.a <- dim(M.a)[2]
        
        P.a <- P1.a <- Pij.a
        finalmat.a <- matrix(0, r.a, c.a)
        
        if(allow.degeneration == F & R.p == C.p){
          
          ## Ensuring that the dimensions of the null matrix will be the same of the original matrix
          
          D.int.finalmat.a <- 0 # The number of rows + columns occupied of the null matrix 
          while (D.int.finalmat.a < sum(dim(M.a))) { # While the dimensions of the null matrix was smaller then the original matrix, keep going
            sel <- sample(1:length(M.a), 1, prob = P.a) # Sample an cell of M.a with probability P.a
            selc <- floor((sel - 1)/(dim(M.a)[1])) + 1 # Recovering column and 
            selr <- ((sel - 1)%%dim(M.a)[1]) + 1 # row of the cell sampled
            if (sum(finalmat.a[, selc]) == 0 | sum(finalmat.a[selr,]) == 0) { # Checking if row or column of the sampled cell is empty
              finalmat.a[sel] <- 1 
              P.a[sel] <- 0
            }
            D.int.finalmat.a <- sum(rowSums(finalmat.a) > 0) + sum(colSums(finalmat.a) > 0) # Setting the new number of dimensions occupied
          }
          # When the number of occupied dimensions of the null matrix was the same as the original matrix, continue
        }
        
        conn.remain <- sum(M.a > 0) - sum(finalmat.a > 0) # The number of cells remaining to be occupied to mantain the original connectance
        
        if (conn.remain > 0) {
          if(connectance == T){
            if (length(which(finalmat.a == 0)) == 1) {
              add <- which(finalmat.a == 0)
            } else {
              add <- sample(which(finalmat.a == 0), conn.remain, 
                            prob = P1.a[finalmat.a == 0], replace = F)
            }
          }else {
            add <- sample(1:length(finalmat.a), conn.remain, 
                          prob = P1.a, replace = T)
          }
          for (add1 in add){
            finalmat.a[add1] <- finalmat.a[add1] + 1
          }
        }
        
        ### Checking if there are still interactions to be drawn. If applicable, draw.
        int.remain <- (sum(M.a) - sum(finalmat.a))
        if (int.remain > 0) {
          if (length(which(finalmat.a > 0)) == 1) {
            add <- rep(which(finalmat.a > 0),int.remain)
          }else{
            add <- sample(which(finalmat.a > 0), int.remain, prob = P1.a[which(finalmat.a >0)], replace = T)
          }
          finalmat.a[as.numeric(names(table(add)))] <- finalmat.a[as.numeric(names(table(add)))] + (table(add)) 
        }
        
        finalmat[R.partitions == R.p, C.partitions == C.p] <- finalmat.a
      }
    }
    
    # Saving outputs
    R2keep <- which(rowSums(finalmat) != 0)
    C2keep <- which(colSums(finalmat) != 0)
    finalmat2 <- finalmat[R2keep,C2keep]
    if (return.nonrm.species == T){
      NullMatrices[[nn]] = list(NullMatrix = finalmat2, R.Kept = R2keep, C.Kept = C2keep)
    }else if(return.nonrm.species == F){
      NullMatrices[[nn]] = finalmat2
    }
    if (Print.null == T){print(nn)}
  }
  return(NullMatrices = NullMatrices)
}




#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = as.matrix(network), 
                     R.partitions = row.Part, #Input the modular structured recovered 
                     C.partitions = col.Part, #Input the modular structured recovered 
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Generate randomized networks with the null model of your choice, considering the interaction probabilities calculated before. 

permutations <- 1000

nulls.com <- RestNullModel(M = as.matrix(network), 
                           Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                           Numbernulls = permutations, #This step may take long, so start experimenting with low values
                           Print.null = T, 
                           allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                           return.nonrm.species = F, 
                           connectance = T, byarea = T, 
                           R.partitions = row.Part, 
                           C.partitions = col.Part)

#Calculate the nestedness within and between modules
rest.nest <- nest.smdm(network, constraints = Part, 
                       weighted = T, 
                       decreasing = "abund", 
                       sort = T)

unlist(rest.nest)

null.com <- sapply(nulls.com, 
                   function(x) bipartite::nest.smdm(x = x,
                                                    constraints = Part, 
                                                    weighted = T, 
                                                    decreasing = "abund"))
WNODA.null.com <- unlist(null.com[3,])
WNODAsm.null.com <- unlist(null.com[8,])
WNODAdm.null.com <- unlist(null.com[9,])


#Estimate the p-values

#Nestedness in the entire network
praw.WNODA <- sum(WNODA.null.com>obs.com[3]) / length(WNODA.null.com)


p.WNODA <- ifelse(praw.WNODA > 0.5, 1- praw.WNODA, praw.WNODA)    # P-value
p.WNODA


((rest.nest$WNp = ODAmatrix)-mean(WNODA.null.com))/sd(WNODA.null.com)

#Nestedness within the modules
praw.WNODAsm <- sum(WNODAsm.null.com>obs.com[8]) / length(WNODAsm.null.com)
p.WNODAsm <- ifelse(praw.WNODAsm > 0.5, 1- praw.WNODAsm, praw.WNODAsm)    # P-value
p.WNODAsm

(nest_n-mean(nestnull_n))/sd(nestnull_n)
((rest.nest$WNODA_SM_matrix)-mean(WNODA.null.com))/sd(WNODA.null.com)

#Nestedness between the modules
praw.WNODAdm <- sum(WNODAdm.null.com>obs.com[9]) / length(WNODAdm.null.com)
p.WNODAdm <- ifelse(praw.WNODAdm > 0.5, 1- praw.WNODAdm, praw.WNODAdm)    # P-value
p.WNODAdm


((rest.nest$WNODA_DM_matrix)-mean(WNODA.null.com))/sd(WNODA.null.com)

par(mfrow = c(1,1))



####################################################################################
################# 4. Morphological analysis and density plots ######################
####################################################################################

# loading morphological data of bats
morph_bats <- read.xls("morph_bats.xlsx", h=T)


# testing for differences between modules
aov_rcr <- aov(morph_bats$rcr~morphdist$module) #Rostrum-to-cranium ratio
summary(aov_rcr) 
TukeyHSD(aov_rcr)

aov_bci <- aov(morphdist$bci~morphdist$module) #body condition index
summary(aov_bci)
TukeyHSD(aov_bci)


# Plotting denstity graphs

nect1 <- morph_bats[morph_bats$module=="nect1",] #separation by modules
nect2 <- morph_bats[morph_bats$module=="nect2",]
frug1 <- morph_bats[morph_bats$module=="frug1",]
frug2 <- morph_bats[morph_bats$module=="frug2",]


# First, RCR

tiff("densitynect1.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect1, aes(x=rcr)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,0.4) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densitynect2.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect2, aes(x=rcr)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,0.4) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug1.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug1, aes(x=rcr)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,0.4) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug2.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug2, aes(x=rcr)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,0.4) + theme_classic() +
    theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))

dev.off()

# Second, BCI

tiff("densitynect1bmi.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect1, aes(x=bci)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,1.8) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densitynect2bmi.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect2, aes(x=bci)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,1.8) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug1bmi.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug1, aes(x=bci)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,1.8) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug2bmi.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug2, aes(x=bci)) +
  geom_density(alpha=.8, color = "lightgoldenrod2", fill="lightgoldenrod2") + 
  xlim(0.1,1.8) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()


# Now, plants

morphplants<-read.xls("morph_plants.xlsx", h=T)

aov_length <- aov(morphplants$length~morphplants$module) #floral tube length
summary(aov_length)
TukeyHSD(aov_length)


aov_girth <- aov(morphplants$girth~morphplants$module) #floral girth/tube diameter
summary(aov_girth)
TukeyHSD(aov_girth)

nect1p <- morphplants[morphplants$module=="nect1",] #separation by module
nect2p <- morphplants[morphplants$module=="nect2",]
frug1p <- morphplants[morphplants$module=="frug1",]
frug2p <- morphplants[morphplants$module=="frug2",]


# First, length

tiff("densitynect1p.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect1p, aes(x=length)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-15,50) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densitynect2p.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect2p, aes(x=length)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-15,50) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug1p.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug1p, aes(x=length)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-15,50) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug2p.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug2p, aes(x=length)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-15,50) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

# Second, diameter

tiff("densitynect1pg.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect1p, aes(x=girth)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-5,40) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densitynect2pg.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(nect2p, aes(x=girth)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-5,40) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug1pg.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug1p, aes(x=girth)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-5,40) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()

tiff("densityfrug2pg.tiff", width = 20, height = 7, units = "cm", res = 600)        
ggplot(frug2p, aes(x=girth)) +
  geom_density(alpha=.8, color = "darkolivegreen3", fill="darkolivegreen3") + 
  xlim(-5,40) + theme_classic() +
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        axis.text= element_text(color="black", size=15),
        legend.position = c(0.85, 0.82),
        legend.title = element_text(face="bold", size = 17),
        legend.text = element_text( size = 15))
dev.off()
  

####################################################################################
################### 5. Likelihood analysis (microstructure) ########################
####################################################################################

# Loading a separate folder with the probability matrices (needs to be set up according to where it was saved in the pc)
setwd("C:/Users/.../Analyses/likelihood models/full")

# The folder contains .txt files (used for analyses) and Excel files (can be used for a more easy multiplication process)

fO <- as.matrix(read.delim("fO.txt")) # full (all bats) observed interaction matrix
fM <- as.matrix(read.delim("fM.txt")) # full morphology interaction matrix
fA <- as.matrix(read.delim("fA.txt")) # full abundance interaction matrix
fS <- as.matrix(read.delim("fS.txt")) # full spatial interaction matrix
fF <- as.matrix(read.delim("fF.txt")) # full temporal interaction matrix
fN <- as.matrix(read.delim("fN.txt")) # full null interaction matrix 
fMF <- as.matrix(read.delim("fMF.txt")) # Multiplications from here onwards...
fMS <- as.matrix(read.delim("fMS.txt")) 
fAF <- as.matrix(read.delim("fAF.txt")) 
fAS <- as.matrix(read.delim("fAS.txt")) 
fSF <- as.matrix(read.delim("fSF.txt")) 
fMFS <- as.matrix(read.delim("fMSF.txt")) 
fAFS <- as.matrix(read.delim("fAFS.txt")) 


############################################################################################
# OBS: The function mlik implemented here was created by/taken from Vázquez et al. (2016), #
# please refer to the original reference for a more thorough cover of the function.        #
#                                                                                          #
# P. Vázquez, Diego; Chacoff, Natacha P.; Cagnolo, Luciano (2016): Evaluating multiple     #
# determinants of the structure of plant-animal mutualistic networks. Wiley. Collection.   #
# https://doi.org/10.6084/m9.figshare.c.3301145.v1                                         #
############################################################################################

# Setting the function to calculate likelihood 
mlik<-function(o,m.p,par){
  lik<--dmultinom(o,prob=m.p,log=T) #Negative log likelihood
  aic=2*lik+2*par
  res=data.frame(lik,aic)
  return(res)
}

# Calculating the likelihood of probability matrices explaining the observed interaction matrix (O), 
# and the respective AIC values

mlik(fO, fO, 28) # Observed (the benchmark)
mlik(fO, fM, 28)
mlik(fO, fA, 28)
mlik(fO, fS, 28)
mlik(fO, fF,28)
mlik(fO, fN,28)
mlik(fO, fMF,56)
mlik(fO, fMS,56)
mlik(fO, fAF,56)
mlik(fO, fAS,56)
mlik(fO, fSF,56)
mlik(fO, fMFS,84)
mlik(fO, fAFS,84)


# plotting results

faic <- read.xls("faic.xlsx",h=T) #sheet containing AIC values from section above

tiff("likelihood_full.tiff", width = 10, height = 15, units = "cm", res = 600)  
ggplot(faic, aes(x=Order, y=delta)) +
  geom_col(col="lightgoldenrod2",fill="lightgoldenrod2",width=0.5) + coord_flip() +theme_classic()+
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        panel.border = element_rect(colour = "black", fill=NA, size=.5), 
        axis.text= element_text(color="black", size=12))
dev.off()



# Now, we repeat with the network containing specialized nectarivores only

setwd("C:/Users/.../Analysis/likelihood models/nectar")

nO <- as.matrix(read.delim("nO.txt"))
nM <- as.matrix(read.delim("nM.txt"))
nA <- as.matrix(read.delim("nA.txt")) 
nS <- as.matrix(read.delim("nS.txt")) 
nF <- as.matrix(read.delim("nF.txt")) 
nN <- as.matrix(read.delim("nN.txt")) 
nMF <- as.matrix(read.delim("nMF.txt")) 
nMS <- as.matrix(read.delim("nMS.txt")) 
nAF <- as.matrix(read.delim("nAF.txt")) 
nAS <- as.matrix(read.delim("nAS.txt")) 
nSF <- as.matrix(read.delim("nSF.txt")) 
nMFS <- as.matrix(read.delim("nMFS.txt")) 
nAFS <- as.matrix(read.delim("nAFS.txt")) 


mlik(nO, nO, 28)
mlik(nO, nM, 28)
mlik(nO, nA, 28)
mlik(nO, nS, 28)
mlik(nO, nF,28)
mlik(nO, nN,28)
mlik(nO, nMF,56)
mlik(nO, nMS,56)
mlik(nO, nAF,56)
mlik(nO, nAS,56)
mlik(nO, nSF,56)
mlik(fO, nMFS,84)
mlik(nO, nAFS,84)

# Plotting
naic <- read.xls("naic.xlsx",h=T)

tiff("likelihood_nectar.tiff", width = 10, height = 15, units = "cm", res = 600)  
ggplot(naic, aes(x=Order, y=delta)) +
  geom_col(col="lightgoldenrod2",fill="lightgoldenrod2",width=0.5) + coord_flip() +theme_classic()+
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        panel.border = element_rect(colour = "black", fill=NA, size=.5), 
        axis.text= element_text(color="black", size=12))
dev.off()


# Finally, repeating the analysis with the network containing bats from other guilds only

setwd("C:/Users/.../Analysis/likelihood models/other")

oO <- as.matrix(read.delim("oO.txt"))
oM <- as.matrix(read.delim("oM.txt"))
oA <- as.matrix(read.delim("oA.txt")) 
oS <- as.matrix(read.delim("oS.txt")) 
oF <- as.matrix(read.delim("oF.txt")) 
oN <- as.matrix(read.delim("oN.txt")) 
oMF <- as.matrix(read.delim("oMF.txt")) 
oMS <- as.matrix(read.delim("oMS.txt")) 
oAF <- as.matrix(read.delim("oAF.txt")) 
oAS <- as.matrix(read.delim("oAS.txt")) 
oSF <- as.matrix(read.delim("oSF.txt")) 
oMFS <- as.matrix(read.delim("oMSF.txt")) 
oAFS <- as.matrix(read.delim("oAFS.txt")) 

mlik(oO, oO, 28)
mlik(oO, oM, 28)
mlik(oO, oA, 28)
mlik(oO, oS, 28)
mlik(oO, oF,28)
mlik(oO, oN,28)
mlik(oO, oMF,56)
mlik(oO, oMS,56)
mlik(oO, oAF,56)
mlik(oO, oAS,56)
mlik(oO, oSF,56)
mlik(oO, oMFS,84)
mlik(oO, oAFS,84)


oaic <- read.xls("oaic.xlsx",h=T)

tiff("likelihood_others.tiff", width = 10, height = 15, units = "cm", res = 600)  
ggplot(oaic, aes(x=Order, y=delta)) +
  geom_col(col="lightgoldenrod2",fill="lightgoldenrod2",width=0.5) + coord_flip() +theme_classic()+
  theme(axis.title.x = element_text(color="black", face ="bold", size =17),
        axis.title.y = element_text(color="black", face ="bold", size =17), 
        panel.border = element_rect(colour = "black", fill=NA, size=.5), 
        axis.text= element_text(color="black", size=12))
abli
dev.off()

