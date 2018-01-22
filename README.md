# Multivariate Macroevolutionary Inference Across Phylogenetic Tree-Building Methods

# Idea and Introduction:  
We want to test macroevolutionary analyses across trees from different phylogenetic reconstruction methods and data subsets to determine how robust the results and subsequent inferences are. We'll initially utilize two methods (1) Disparity through time 'DTT' and (2) l1ou to testing for shifts in morphological traits. The input trait data will stay the same for all analyses, but data subsets and methods will change.

## Method 1 is 'Concatenation':  
Concatenation forces a single topology for all loci included. It means we can use all the data (good) and it's fast (great), but it lacks fundamental biological reality when using tens or hundreds of loci. To get our concatenated trees we can either use RAxML or BEAST2. BEAST2 is preferrable because it will output a distribution of plausible ultrametric trees of varied topologies and branch lengths all at once. Alternatively, we can run RAxML for 100 tree searches, which will give us some branch length/phylogenetic variation, but we will need to transform the tree to make it ultrametric ('chronos', lambda = ...). Currently I'm trying the latter, as it should be the most different to alternate methods.

## Method 2 is 'Short-Cut Coalescent':  
To maximize the amount of data we can use (all!) and computational time (speedy!), we can use shortcut coalescent analyses (SCA). These methods (ASTRAL, ASTRID, MP-EST) use gene trees estimated from another program as inputs. SCA methods can be sensitive to the accuracy/quality of input gene trees, so we can test this by using gene trees built in RAxML and in starBEAST2. If we run a bootstrapped ASTRAL analysis, we will get an output of 102 trees (100 bs, 1 greedy algorithim best fit, and 1 consensus), perfect for our analyses!

## Method 3 is 'Full Coalescent':  
The only proper full coalescent method that we'll use is implemented in starBEAST2. Because it is computationally expensive, we will run analyses in batches of 10-25 loci, limiting loci to only those including all taxa (no missing data). Then we can use batches of MCC trees or something more clever to use for our analyses.
