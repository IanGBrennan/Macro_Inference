#####
# Script for macro-comparison of shifts in multivariate traits
## inferred from different phylogenetic methods (ML, BI, MSC, SCC)

###################################################################################
## Idea and Introduction:
# we want to test macroevolutionary analyses across trees from different methods
# and data subsets to determine how robust the results and subsequent inferences
# are. We'll initially try two things (1) Disparity through time 'DTT' and (2) l1ou
# to test for shifts in morphological traits. The input trait data will stay the
# same for all analyses, but data subsets and methods will change.

# Method 1 is 'Concatenation': We can either use RAxML or BEAST2. BEAST2 is preferrable
# because it will output a distribution of plausible ultrametric trees of varied topologies
# and branch lengths all at once. Alternatively, we can run RAxML for 100 tree
# searches, which will give us some branch length/phylogenetic variation, but we
# will need to transform the tree to make it ultrametric ('chronos', lambda = ...).

# Method 2 is 'Short-Cut Coalescent': These are methods (ASTRAL, ASTRID, MP-EST) 
# that use as inputs, gene trees estimated from another program (I will use either
# RAxML or *BEAST2). These methods can be sensitive to the input gene trees, so we
# can test this by using gene trees built in RAxML and in BEAST2. If we run a
# bootstrapped ASTRAL analysis, we will get an output of 102 trees (100 bs,
# 1 greedy algorithim best fit, and 1 consensus), perfect for our analyses!

# Method 3 is 'Full Coalescent': The only proper full coalescent method that we'll
# use is implemented in *BEAST2. 

library(geiger)
library(phytools)
library(ggbiplot);library(ggplot2);library(ggtree); library(ggridges)
library(l1ou)
library(Rmisc)
library(wesanderson); library(ggthemes)

# A list of the final sets of trees/data we're using in this study 
#############################################################
## (all in "/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/")
# or ("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/"):
  # Marsupials
      # RAxML concatenated = "Marsupials.RAxML.Concat.trees"
          # phylo regressed traits = "Marsupials.RAxML.Concat.PhyloResiduals.rds"
      # ASTRAL w/RAxML gene trees = "Marsupials.ASTRAL.RAxML.trees"
          # phylo regressed traits = "Marsupials.Astral.RAxML.PhyloResiduals.rds"
      # starBEAST = 
      # ASTRAL w/starBEAST gene trees =
  # Elapids
      # RAxML concatenated = "T222_concatenated_Ingroup.trimmed.trees"
          # phylo regressed traits = "Elapids.Concat.RAxML.PhyloResiduals.rds"
      # ASTRAL w/RAxML gene trees = "Elapidae.RAxML.ASTRAL.bs.SCALED.100.TRIMMED.trees"
          # phylo regressed traits = "Elapids.Astral.RAxML.PhyloResiduals.rds"
      # starBEAST = "Elapidae.*BEAST.by25s.TRIMMED.trees"
          # phylo regressed traits = "Elapids.BEAST.PhyloResiduals.rds"
      # ASTRAL w/starBEAST gene trees = "Elapidae.ASTRAL.*BEAST.TRIMMED.trees"
          # phylo regressed traits = "Elapids.Astral.BEAST.PhyloResiduals.rds"
      # phylogenetically regressed traits (for l1ou) = 
  # Protea
      # RAxML concatenated = "Protea.RAxML.Concat.tre"
      # ASTRAL w/RAxML gene trees = "Protea.ASTRAL.RAxML.SCALED.tre"
      # starBEAST = "Protea.starBEAST.trimmed.trees"
      # ASTRAL w/starBEAST gene trees = "Protea.ASTRAL.starBEAST.TRIMMED.SCALED.trees"
      # phylogenetically regressed traits (for l1ou) = 
  # Cichlids


# Just a staging area for TREEs:
#################################
trees = read.nexus("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/Marsupials.ASTRAL.RAxML.trees") #pick out our set of posterior trees
#################################


# Quickly standardize the data to body size
#######################################################################
#### I've decided to use a phylogenetic regression instead of a standard linear regression, 
##### the steps for running it for each tree are looped below
datum <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.RAW.csv", header=T, row.names=1)
datum <- log(datum) # quickly log the raw values (if not already done!)
total.data <- NULL # create a frame for holding all our residuals
against <- datum$M.avgBL; names(against) <- rownames(datum) # phyl.resid needs names to determine proper order
for (b in 1:length(trees)) {
  resid.data <- NULL # temporary object
  resids <- phyl.resid(trees[[b]], against, datum[,c("Brain.size", "M.avgWT")]) # regress brain size and weight against body length
  resid.data <- cbind(resid.data, resids$resid); resid.data <- cbind(resid.data, against); 
  colnames(resid.data) <- c("BrainSize", "Weight", "BodyLength")
  total.data[[b]] <- resid.data
}
saveRDS(total.data, "/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.ASTRAL.RAxML.PhyloResiduals.rds")

# You'll want to go through these steps each time you start again
trait.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.residuals.data.csv")
#trait.data <- subset(trait.data, trait.data$In_analyses == "Yes")
#trait.data <- trait.data[complete.cases(trait.data),] # trim any data columns that aren't complete
trait.data <- trait.data[,c("Name_in_Duchene", "logLength", "brain_res", "weight_res")]
trait.data <- trait.data[,c("Name_in_Duchene", "logLength", "nonlogbrain_res", "nonlogweight_res")]
rownames(trait.data) <- trait.data$Name_in_Duchene
# or
trait.data <- readRDS("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.PhyloResiduals.rds")


# Read in the trees we'll use:
astral.raxml <- read.tree ("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/Marsupials_ASTRAL_RAxML.trees")
raxml.concat <- read.tree ("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/")
starbeast    <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/")
astralbeast  <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/")
# if you've already cleaned up your trees assign them to the 'trees' object
trees = astral.raxml

# Add on the Thylacine next to Myrmecobius at a reasonable depth
depths <- seq(0.3,0.9, 0.01)
out.trees <- NULL
for (k in 1:length(trees)){
  at.depth <- sample(depths, 1)
  out.trees[[k]] <- bind.tip(trees[[k]], where=which(trees[[k]]$tip.label=="Myrmecobiidae_Myrmecobius_fasciatus"),
                         tip.label="Thylacinidae_Thylacinus_cynocephalus", 
                         position=at.depth*trees[[k]]$edge.length[which(trees[[k]]$edge[,2]==which(trees[[k]]$tip.label=="Myrmecobiidae_Myrmecobius_fasciatus"))])
}



# then we can clean up your tree(s) to match your data; unnecessary if using 'TRIMMED' trees
###############################################################
keep <- rownames(trait.data)
all <- astral.raxml[[1]]$tip.label # change to appropriate tree
drop <- setdiff(all, keep)
new.trees <- lapply(astral.raxml, drop.tip, tip=drop) # change to appropriate tree
class(new.trees) <- "multiPhylo"

name.check(new.trees[[1]], trait.data)
trees = new.trees

# now we can get into the analyses, start with DTT in geiger
##############################################################
# either use a raw trait
trait <- trait.data$Total_L # specify your trait of interest, make sure it has names attached!
names(trait) <- rownames(trans.data) # make sure it has names attached!
test <- dtt(trees[[1]], trait, plot=T, nsim=1000)
plot(test$dtt ~ test$times)
lines(test$times, test$dtt)
# or a PCA value, which is reducing the dimensionality of multiple traits
## if you want to do a PCA value (which you likely do, follow the steps below)
pca.test <- trait.data[complete.cases(trait.data),] # make sure to remove incomplete taxa
test.data <- pca.test[,2:length(trait.data)] # designate the variables of interest (columns)
#test.data <- log(test.data) # log transform data
#ln.test[ln.test=="-Inf"] <- 0
species <- pca.test[,1] # make note of the species names
#genus <- pca.test[,7]
#parity <- pca.test[,8]
trait.pca <- prcomp(test.data) # perform your PCA
plot(trait.pca, type="l") # visualize how many axes to keep
trait.pca # determine how much influence each axis holds
summary(trait.pca) # and the cumulative value of each
#(ggbiplot(protea.pca, obs.scale=1, var.scale=1, ellipse=T, circle=T))
#loadings <- as.data.frame(protea.pca$rotation)
axes <- predict(trait.pca, newdata = test.data) # pull out each taxon's value for each PC axis
trait.data[,"PC1"] <- axes[,1] # save the PC loadings as a new variable for the DTT
(ggbiplot(trait.pca, obs.scale=1, var.scale=1, ellipse=T, circle=T))

### now you can do your DTT on PC1!
trait <- trait.data[,"PC1"] # specify your trait of interest, make sure it has names attached!
names(trait) <- rownames(axes) # make sure it has names attached!
test <- dtt(trees[[1]], trait, plot=T, nsim=1000)
plot(test$dtt ~ test$times)
lines(test$times, test$dtt)



# otherwise, just jump right in and do this analysis for all the trees!
par(new=F)
mdi.estimates <- NULL
timing.of.disparity <- NULL
for (i in 1:length(trees)) {
  cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
  input <- dtt(trees[[i]], trait, plot=F, nsim=1000) # run the DTT analysis
  
  # create a set of objects to hold the outputs
  mdi.info <- as.data.frame(t(c(i, input$MDI)))
  disparity <- as.data.frame(input$dtt)
  timing <- as.data.frame(input$times)
  sim.means <- as.data.frame(rowMeans(input$sim)) # keep this for plotting inside this function
  sim.bounds <- as.data.frame(t(apply(input$sim, 1, CI))) # apply the CI function to each row (time interval) of the simulations
  colnames(sim.bounds) <- c("simulated.upper95", "simulated.mean", "simulated.lower95")
  
  # create an object to hold the combined timing/disparity data, and feed it into a bigger frame
  disp.timing <- NULL
  disp.timing <- cbind(timing, disparity, sim.bounds) # could use sim.means if you prefer (less uncertainty)
  timing.of.disparity <- rbind.data.frame(timing.of.disparity, disp.timing)
  
  # plot only the empirical trends (not the simulated)
  plot(input$dtt ~ input$times,
       xlim=c(0,1), ylim=c(0,2), col="red")
  lines(input$times, input$dtt); par(new=T)
  # plot only the simulated mean trends
  plot(sim.means[,1] ~ input$times,
       xlim=c(0,1), ylim=c(0,2), col="blue")
  lines(input$times, sim.means[,1]); par(new=T)
  
  mdi.estimates <- rbind(mdi.estimates, mdi.info)
}
colnames(timing.of.disparity) <- c("timing", "disparity", "simulated.upper95", "simulated.mean", "simulated.lower95")
timing.of.disparity <- timing.of.disparity[order(timing.of.disparity$timing),]

## Plot the timing of shifts as a density distribution (with mean estimate)
colnames(mdi.estimates) <- c("tree.num", "MDI")
concat <- (ggplot(mdi.estimates, aes(x=MDI)) 
           + geom_density(fill="green")
           + scale_x_continuous(limits = c(0, 1))
           + theme_classic())
#+ geom_vline(aes(xintercept=mean(timing, na.rm=T)),
#             color="red", linetype="dashed", size=1)
#+ scale_x_reverse(limits=c(20,0))) # use this to reverse then define the limits of the x axis

#### if you want to add another method to the MDI estimates plot:
total.mdi <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/All_Methods.MDI.estimates.csv", header=T)
total.mdi[,"X"] <- NULL
mdi.estimates[,"model"] <- "starBEAST"
total.mdi <- rbind(total.mdi, mdi.estimates)
write.csv(total.mdi, file="/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/All_Methods.MDI.estimates.csv")
total.mdi$model <- factor(total.mdi$model, levels=c("starBEAST", "raxml.concatenated", "astral.raxml", "astral.starbeast")) # this re-orders the models in the legend
(ggplot(total.mdi, aes(x = MDI, y = model, fill=model)) 
  + geom_density_ridges(scale=7)
  + scale_fill_manual(values=wes_palette("Zissou", 3, "discrete"))
  + theme_ridges())

#### Make a summary of the Disparity estimates within a sliding window (mean, lowerCI, upperCI)
tof <- timing.of.disparity
sim.nums <- seq(0,1,0.01)
emp.mean.estimates <- NULL; sim.mean.estimates <- NULL
for (t in sim.nums) {
  time.min <- t
  time.max <- t+0.1 # adjust the window width if you're getting NaN values (no observations within a time period)
  
  timed.chunk <- subset(tof, tof$timing>time.min & tof$timing<time.max)
  emp.estimates <- as.data.frame(t(CI(timed.chunk$disparity)))
  #simu.estimates <- as.data.frame(t(CI(timed.chunk$simulated.mean))) # this takes the CIs from the mean of simulations
  sim.estimates <- as.data.frame(t(apply(timed.chunk[,3:5], 2, mean))) # this takes the mean of the CIs from the simulations
  timing <- paste(time.min, "to", time.max); colnames(timing)
  emp.output.frame <- NULL; sim.output.frame <- NULL
  emp.output.frame <- as.data.frame(c(timing, emp.estimates))
  sim.output.frame <- as.data.frame(c(timing, sim.estimates))
  colnames(emp.output.frame) <- c("time.window", "empirical.upper95", "empirical.mean", "empirical.lower95")
  colnames(sim.output.frame) <- c("time.window", "simulated.upper95", "simulated.mean", "simulated.lower95")
  #colnames(output.frame) <- NULL
  emp.mean.estimates <- rbind(emp.mean.estimates, emp.output.frame)
  sim.mean.estimates <- rbind(sim.mean.estimates, sim.output.frame)
}
emp.mean.estimates[,"timing"] <- sim.nums # add a column with the times
emp.mean.estimates <- emp.mean.estimates[-c(101),] # drop the last row which is superfluous
sim.mean.estimates[,"timing"] <- sim.nums # add a column with the times
sim.mean.estimates <- sim.mean.estimates[-c(101),] # drop the last row which is superfluous
total.estimates <- cbind(emp.mean.estimates, c(sim.mean.estimates[,2:4]))

#external.total <- NULL
#read.csv(external.total) # read in the results from other tree sets
total.estimates[,"model"] <- "starBEAST" # change this according to the trees you're using
external.total <- rbind(external.total, total.estimates)
write.csv(external.total, file="/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/All_Methods.DTT.trends.csv")

#### Now let's plot the trend in disparity from the posterior of trees as a confidence ribbon
rib <- (ggplot(data=total.estimates)
        + geom_ribbon(aes(x=timing, ymin=empirical.lower95, ymax=empirical.upper95, fill="empirical"))
        + geom_ribbon(aes(x=timing, ymin=simulated.lower95, ymax=simulated.upper95, fill="simulated"))
        + theme_classic())
#+ geom_smooth()
#+ scale_x_reverse())
#+ geom_smooth(method="auto", aes(x=x, y=all.emp.ci.vj$y.meanCI), se=T))

#### Now combine the MDI plot and the tree into a single figure
test <- (ggtree(trees[[1]]) 
         + geom_tiplab(size=4))

multiplot(rib, test, rib, test, ncol=2)


mardata <- cbind(mardata, log(mardata[,2:4])) # log transform the raw measurements
colnames(mardata) <- c("Name_in_Duchene", "Brain.size", "M.avgBL", "M.avgWT",
                       "logBrain", "logLength", "logWeight")

trait.data
logbrain_res <- log(trait.data[,"nonlogbrain_res"])

########################################################
# Now we can start looking at shifts in size and shape
## using 'l1ou' we'll estimate morphological shifts
## then attempt to identify instances of convergence
########################################################

#### Check to make sure the tips match the data labels
name.check(trees[[23]], trait.data[[23]]);

#### Adjust the data and tree to fit (order matters in l1ou!)
data <- adjust_data(trees[[16]], trait.data[[16]]) # exclude your PC1 values from above!


#### Estimate the number and position of shifts a priori 
shift.fit <- estimate_shift_configuration(data$tree, data$Y, nCores=8, quietly=F, criterion="pBIC") # if you want to do only a single trait, 'data$Y[,x]'
plot(shift.fit, cex=0.8)
# if you get 'figure margins' error, do: par(mar=c(1,1,1,1))

#### Investigate convergence among the supported shift set
fit.conv <- estimate_convergent_regimes(shift.fit, nCores=8, criterion="pBIC")
plot(fit.conv, cex=0.5)


## This quick function pulls out the descendant tips from and edge number
#########################################################################
getDescendants.edges<-function(tree,edge,curr=NULL){
  names <- NULL
  if(is.null(curr)) curr<-vector()
  node.below <- tree$edge[edge,2]
  if(node.below <= Ntip(tree)) {
    input <- tree$tip.label[[node.below]]
    names <- append(names, input)
  }
  else {
    daughters<-tree$edge[which(tree$edge[,1]==node.below),2]
    curr<-c(curr,daughters)
    z<-which(daughters<=length(tree$tip))
    if(length(z)==2) for(i in 1:length(z)) {
      input <- tree$tip.label[[curr[[i]]]]
      names <- append(names, input)
    }
    if(length(z)==1) {
      target <- daughters[[z]]
      input <- tree$tip.label[[target]]
      names <- append(names, input)
    }
    w<-which(daughters>=length(tree$tip))
    if(length(w)>0) for(i in 1:length(w)) 
      curr<-getDescendants(tree,daughters[w[1]],curr)
    curr<-unique(curr)
    curr<-subset(curr, curr<=Ntip(tree))
    for (q in 1:length(curr)) {
      input <- tree$tip.label[[curr[[q]]]]
      names <- append(names, input)
    }
  }
  names <- unique(names)
  return(names)
}
#########################################################################

## Let's try building a loop to make sense of the shifts
########################################################
# *note, edge indices are rearranged in the 'adjust_data' function,
# so to get the proper edge index, call from data$tree
no.shifts <- NULL
shift.positions.by.tree <- list()
shift.positions.list <- NULL
l1ou.res <- NULL

# before you hit enter, make sure to change the names of the output files below!

for (i in 1:length(trees)) {
  cat("iteration", i, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  
  #### Adjust the data and tree to fit (order matters in l1ou!)
  data <- adjust_data(trees[[i]], trait.data[[i]]) # specify the data columns of interest
  
  #### Estimate the number and position of shifts a priori 
  shift.fit <- estimate_shift_configuration(data$tree, data$Y, nCores=8, quietly=F, criterion="pBIC") # if you want to do only a single trait, 'data$Y[,x]'
  l1ou.res[[i]] <- shift.fit
  #plot(shift.fit, cex=0.8)
  # if you get 'figure margins' error, do: par(mar=c(1,1,1,1))
  
  #### Create a data frame to hold the tree # and the # of shifts inferred
  shifts.frame <- NULL
  shifts.frame <- as.data.frame(t(c(i, shift.fit$nShifts)))
  no.shifts <- rbind(no.shifts, shifts.frame)
  
  # have to pull the shift positions (edges)
  shift.edges <- shift.fit$shift.configuration
  # match it to tips
  all.shifted.tips <- NULL
  if (length(shift.edges) == 0) {
    all.shifted.tips <- "no shifts"
  } else for (t in 1:length(shift.edges)) {
              names <- getDescendants.edges(data$tree, shift.edges[[t]])
              all.shifted.tips <- append(all.shifted.tips, names)
  }
  shift.positions.list <- append(shift.positions.list, all.shifted.tips)
  shift.positions.by.tree[[i]] <- all.shifted.tips
  
  saveRDS(no.shifts,               file="/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/l1ou_Output/Marsupials.ASTRAL.RAxML.num.shifts.RDS")
  saveRDS(shift.positions.list,    file="/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/l1ou_Output/Marsupials.ASTRAL.RAxML.list.shift.positions.RDS")
  saveRDS(shift.positions.by.tree, file="/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/l1ou_Output/Marsupials.ASTRAL.RAxML.shift.positions.by.tree.RDS")
  saveRDS(l1ou.res,                file="/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/l1ou_Output/Marsupials.ASTRAL.RAxML.Results.RDS")
}
colnames(no.shifts) <- c("tree.no", "n.shifts")

# add.lengths to appropriate edges of the tree
res.counts <- table(shift.positions.list) # make a table of the shift frequencies
shifted.tips <- as.data.frame(res.counts) # turn it into a data frame
# shifted.tips <- shifted.tips[1,] # remove any labelled "no shifts"
all.node.numbers <- as.data.frame(trees[[1]]$tip.label) # get all the nodes of the tree and the numbers, the tree must match the one you want to plot!
all.node.numbers[,"tip.no"] <- rownames(all.node.numbers); colnames(all.node.numbers) <- c("tip.name", "tip.no") # make a column that shows the tip number
target.numbers <-  all.node.numbers[all.node.numbers$tip.name %in% shifted.tips$shift.positions.list,] # subset all the tips, so show just the node numbers of the shifted tips
target.numbers[,2] <- as.numeric(target.numbers[,2]) # make it numeric
target.numbers <- target.numbers[order(match(target.numbers$tip.name, shifted.tips$shift.positions.list)),] # match the new frame to the output (shift) frame
target.numbers[,"shift.freq"] <- cbind(shifted.tips$Freq) # add on the frequencies of shifts

# now we need to designate the tree we want to plot, and adjust its branch lengths to match the shifts
chosen.shift <- l1ou.res[[1]] # designate which shift set you want to plot the adjustments on 
tree <- chosen.shift$tree
#alt.data <- adjust_data(trees[[60]], protea.data[,c(2:6)]) # it has to be adjust to fit the 'post order'
#tree <- alt.data$tree # designate your target tree, I've put in the tree from the last 'shift.fit' object or do "tree <- trees[[10]]"
#tree <- shift.fit$tree # designate your target tree, I've put in the tree from the last 'shift.fit' object
#shift.fit$shift.configuration <- l1ou.res[[90]]$shift.configuration # bring along the shift config from the tree you want to plot!

for (i in 1:length(target.numbers[,2])){
  rownames(tree$edge) <- c(1:length(tree$edge[,1])) # give the tree edge frame rownames
  target.almost <- subset(tree$edge, tree$edge[,2]==target.numbers[,2][i]) # pull out the ancestor and descendant nodes of the target edge
  interim.target <- subset(target.numbers, target.numbers[,2]==target.numbers[,2][i]) # subset the data frame to just the current tip of interest (descendant node)
  target.edge <- as.numeric(rownames(target.almost)) # get the number of the target edge
  tree$edge.length[[target.edge]] <- tree$edge.length[[target.edge]]+(0.01*interim.target[,3]) # add the desired length to the branch, per shift (here, 0.01)
}
chosen.shift$tree <- tree # set the tree in the 'shift.fit' object to our rescaled tree
#shift.fit$Y <- alt.data$Y
#plot(tree) # have a look to see if it worked
# if you get 'figure margins' error, do: par(mar=c(1,1,1,1))
plot(chosen.shift, cex=1)


#####################
## I've now saved the l1ou output as it's own object "l1ou.res", so we should be able to 
### call on it to plot any tree we want without having to go back. We'll try it next time
####################

shifted.tips$method <- "ASTRAL.starBEAST"

shifted.tips$colors <- cut(shifted.tips$Freq,
                           breaks = c(0, 50, 100, 200),
                           labels = c("noise", "real", "reallyreal"))
ggplot(shifted.tips, aes(x=shift.positions.list, y=Freq)) +
  geom_bar(stat="identity")


# if we want to join the shift.positions lists together, read them in and combined
load("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/l1ou_Outputs_and_Figures/Protea.starBEAST.list.shift.positions.RData")
sbeast.shifts <- shift.positions.list
sbeast.list <- as.data.frame(count(sbeast.shifts))
sbeast.list$method <- "starBEAST"
colnames(sbeast.list) <- c("shift.positions.list","Freq", "method")
# then we can export the document
all.shifted.tips <- rbind.data.frame(all.shifted.tips, sbeast.list)
write.csv(all.shifted.tips, file="/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/ShiftedTips_by_Method.csv")

# and read it back in if necessary, so that we can plot it
all.shifts.list <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/ShiftedTips_by_Method.csv")
ggplot(all.shifts.list, aes(x=shift.positions.list, y=Freq, fill=method)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.5), width=2) +
  #theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual( values=wes_palette("Moonrise3")) +
  coord_flip() +
  #scale_color_ptol() +
  theme_few()




###################################################################################
# One more thing we might want to show is the branch length difference among trees
## we can try to do this by plotting pairwise distances to show inherent bias
###################################################################################
# start by rescaling the trees to height 1
as <- lapply(astralbeast, rescale, "depth", 1); sb <- lapply(starbeast, rescale, "depth", 1)
class(as) <- "multiPhylo"; class(sb) <- "multiPhylo"
# now a loop to compare tree1/method1 to tree1/method2 via pairwise distances
total.dff <- NULL # make an empty object to store all distances (ntips! x ntrees)
method1 <- starbeast; method2 <- astralbeast
for (k in 1:length(starbeast)) {
  dff <- NULL
  in1 <- method1[[k]] # designate tree.n/method.n
  in2 <- method2[[k]] # designate tree.n/method.n+1
  
  inboth <- intersect(in1$tip.label, in2$tip.label) # check all the tips that match between trees
  
  # Get the pairwise distance matrices
  pw.in1 <- cophenetic.phylo(in1)
  pw.in2 <- cophenetic.phylo(in2)
  
  # now: compare pairwise distances from these 2 trees
  #   the complication is that they have different sets of taxa
  
  # Get unique combinations of this set:
  ucomb <-  combn(inboth, m = 2)
  
  # make vectors to hold results
  dist_1 <- rep(NA, ncol(ucomb))
  dist_2 <- rep(NA, ncol(ucomb))
  
  dff <- data.frame(species1 = ucomb[1,], species2 = ucomb[2,] , dist_1, dist_2, stringsAsFactors=F)
  
  # fill in the blanks....
  
  for (ii in 1:nrow(dff)){
    dff$dist_1[ii] <- pw.in1[ dff$species1[ii], dff$species2[ii] ]
    dff$dist_2[ii] <- pw.in2[ dff$species1[ii], dff$species2[ii] ]  
  }
  total.dff <- rbind.data.frame(total.dff, dff)
}

fit <- lm(dist_1 ~ dist_2, data=total.dff) # change this according to the parameter you simulated
plot.fit <- (ggplotRegression(fit))
plot.fit + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

#######################################################################################
# Interlude: the base 'plot' and 'abline' functions are alright, but we want to 
## make it (1) prettier, and (2) include the information from our linear regression
### into the plot, so that we know what our results were. Use custom 'ggplotRegression'
### if you want to change the saturation use 'alpha'
ggplotRegression <- function (fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(alpha=0.25, color="red") + # change to 0.25 and "red" for time plots
    stat_smooth(method = "lm", col = "black") + # change to "black" for time plots
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
#######################################################################################












wd = "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/"

# Run a loop that writes each one to an individual file in a named folder
dir.create(np(paste(addslash(wd), "Protea.RAxML.GeneTrees", sep=""))) #create folder for the trees
for (i in 1:length(trees)){
  name <- paste(wd,"/Protea.RAxML.gene.",i,".tre", sep="")
  write.tree(trees[[i]], file=name)
} #you should now have a folder with 100 tree separate tree files




butt <- load(file="Protea.ASTRAL.RAxML.list.shift.positions.RData")
View(butt)





# Quickly standardize the data to body size
#######################################################################
# regress (log transformed) Tail Length against SVL, extract residuals
length.weight <- lm(logWeight ~ logLength, data=datum)
length.weight <- lm(M.avgWT ~ M.avgBL, data=datum)
#plot(length.weight)
nonlogweight_res <- resid(length.weight)
mardata <- cbind(mardata, nonlogweight_res)

# regress (log transformed) Head Length against SVL, extract residuals
length.brain <- lm(logBrain ~ logLength, data=mardata)
length.brain <- lm(Brain.size ~ M.avgBL, data=datum)
#plot(length.brain)
nonlogbrain_res <- resid(length.brain)
mardata <- cbind(mardata, nonlogbrain_res)

datum <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.RAW.csv", header=T, row.names=1)
datum <- log(datum) # quickly log the raw values (if not already done!)
total.data <- NULL # create a frame for holding all our residuals
against <- datum$M.avgBL; names(against) <- rownames(datum) # phyl.resid needs names to determine proper order
for (b in 1:length(trees)) {
  resid.data <- NULL # temporary object
  resids <- phyl.resid(trees[[b]], datum$M.avgBL, datum[,c("Brain.size", "M.avgWT")]) # regress brain size and weight against body length
  resid.data <- cbind(resid.data, datum$M.avgBL); resid.data <- cbind(resid.data, resids$resid)
  colnames(resid.data) <- c("BodyLength", "BrainSize", "Weight")
  total.data[[b]] <- resid.data
}
saveRDS(total.data, "/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.RAxML.Concat.PhyloResiduals.rds")

butt <- readRDS("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.Astral.RAxML.PhyloResiduals.rds")
trees = read.tree("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/Marsupials.ASTRAL.RAxML.trees") #pick out our set of posterior trees
test <- phyl.resid(trees[[1]], against, datum[,c("Brain.size", "M.avgWT")])
test
butt[[1]]
