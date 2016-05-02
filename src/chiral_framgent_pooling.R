options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/fragments')
library(rcdk)

pool_size = 5

# function to quickly just look at a structure in R interactive mode
plotmol = function(molecule,width=500,height=500) {
  par(mar=c(0,0,0,0)) # set margins to zero since this isn't a real plot
  temp1 = view.image.2d(molecule,width,height) # get Java representation into an image matrix. set number of pixels you want horiz and vertical
  plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='') # create an empty plot
  rasterImage(temp1,1,1,10,10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
}

# create a .-delimited list of SMILES for copying and pasting into ChemDraw
sepsmiles = function(smiles) {
  return (paste(smiles,collapse='.'))
}

# read in data
frags = read.table('data/cbip_potential_fragments.tsv',header=T,sep='\t',quote='',comment.char='',allowEscapes=FALSE)
colnames(frags) = tolower(gsub("[^A-Za-z0-9_]","_",colnames(frags)))
frags = frags[frags$expected_mass <= 300,]

# borrowed some code originally employed here: https://github.com/ericminikel/prp_knockdown_screens/blob/master/src/analysis.R
fps = list(nrow(frags))
for (i in 1:nrow(frags)) {
  mol = parse.smiles(frags$smiles[i])[[1]]
  fp = get.fingerprint(mol, type='circular', fp.mode = 'bit', depth=6, size=1024, verbose=FALSE) # circular==ECFP6
  fps[i] = fp
}

# create distance matrix and hierarchical clustering, which we never end up using
dist_matrix = 1 - fp.sim.matrix(fplist=fps, method='tanimoto')
clustering = hclust(as.dist(dist_matrix))
plot(clustering)

set.seed(2) # 2 is the best I've looked at

# first, create random clusters
frags$rand = runif(n=nrow(frags),min=0,max=1)
frags$random_ordering = rank(frags$rand)
frags$random_cluster = frags$random_ordering %/% pool_size
frags$reshuffled_cluster = frags$random_cluster

# then try to optimize the clusters by re-shuffling the fragments that ended up in a group with something similar to themself
for (iteration in 1:10) {
  frags$mindist = 0
  for (i in 1:nrow(frags)) {
    # minimum distance of this fragment to any other in its cluster, not including itself
    frags$mindist[i] = min(dist_matrix[i,frags$reshuffled_cluster==frags$reshuffled_cluster[i] & 1:nrow(frags)!=i])
  }
  # lowest minimum distance of any fragment on this round
  print(paste("iteration ",iteration,": ",min(frags$mindist)))
  
  # quality checks and visualizations
  # hist(dist_matrix)
  # hist(frags$mindist)
  # check that all clusters should have either 0 or >=2 "close" fragments
  #   frags$random_cluster[frags$mindist < .5]
  #   frags[frags$random_cluster==3,]
  
  # re-shuffle all the fragments that are <0.7 similar to their nearest fragment in their cluster
  to_reshuffle = frags$mindist < .7
  frags$reshuffled_cluster[to_reshuffle] = sample(frags$reshuffled_cluster[to_reshuffle],size=sum(to_reshuffle),replace=FALSE)
}

# view the new distribution
hist(frags$mindist)

# select clusters that have (relatively) close similarity within them
frags$reshuffled_cluster[frags$mindist < .75]

# view such a cluster
sepsmiles(frags$smiles[frags$reshuffled_cluster==61])

# write out results
output = frags[with(frags, order(reshuffled_cluster)),]
plate_layout = data.frame(well_name=paste(rep(LETTERS[1:8],each=12),rep(1:12,8),sep=''),well_number=0:95)
output$well_name = plate_layout$well_name[match(output$reshuffled_cluster,plate_layout$well_number)]
output = output[,c("well_name","broad_id","smiles","expected_mass__desalted_","chemist_or_library")]
write.table(output,'chiral_fragment_pools.tsv',sep='\t',row.names=F,col.names=T,quote=F)
