require(ape)
require(rstan)
#require(cmdstanr)

#error regarding mde(x): convert data frame to matrix

model.path <- paste(commandArgs(trailingOnly = TRUE)[1],'.stan',sep='')

curr.seed <- commandArgs(trailingOnly = TRUE)[2]

fams <- gsub('.posterior.tree','',dir('posteriorTrees'))

set.seed(curr.seed)

tree.index <- sample(1000,1)

families <- c()
isolates <- c()
smallfamilies <- c()
for (fam in fams) {
  trees <- read.tree(paste('posteriorTrees/',fam,'.posterior.tree',sep=''))
  if (class(trees) == 'multiPhylo') {
    if (length(trees[[1]]$tip.label)==2) {
      smallfamilies <- c(smallfamilies,fam)
    }
    else {
      families <- c(families,fam)
    }
  }
  else {
    isolates <- c(isolates,fam)
  }
}

data.df <- read.csv('charMtx.csv')
all.states <- c()
all.b.lens <- c()
all.parent <- c()
all.child <- c()
all.envir <- c()
roots <- c()
fam.ID <- c()
is.tip <- c()
max.node <- 0
max.b.len <- 0
mrca.matrix <- matrix(0,nrow=1391,ncol=1391)
tree.ind <- 1
tree.sizes <- c()
tree.sizes.nodes <- c()
node.heights <- c()
for (fam in families) {
  trees <- read.tree(paste('posteriorTrees/',fam,'.posterior.tree',sep=''))
  tree <- trees[[tree.index]]
  data.df.sub <- data.df[data.df$longname%in%tree$tip.label,]
  rownames(data.df.sub) = data.df.sub$longname
  data.df.sub <- data.df.sub[,3:ncol(data.df.sub)]
  tree <- reorder.phylo(tree,'pruningwise')
  states <- data.df.sub[tree$tip.label,]
  states <- cbind(states,1-states)
  node.states <- as.data.frame(matrix(1,nrow=tree$Nnode,ncol=ncol(states)))
  colnames(node.states) <- colnames(states)
  states <- rbind(states,node.states)
  states <- states[,c(rbind(c(1:(ncol(states)/2)),c(1:(ncol(states)/2))+(ncol(states)/2)))]
  parent <- tree$edge[,1] + max.node
  child <- tree$edge[,2] + max.node
  b.lens <- tree$edge.length
  all.states <- rbind(all.states,states)
  all.b.lens <- c(all.b.lens,b.lens)
  all.parent <- c(all.parent,parent)
  all.child <- c(all.child,child)
  max.node <- max(parent)
  roots <- c(roots,length(parent) + max.b.len)
  max.b.len <- max.b.len + length(parent)
  fam.ID <- c(fam.ID,rep(fam,length(b.lens)))
  is.tip <- c(is.tip,c(rep(1,length(tree$tip.label)),rep(0,tree$Nnode)))
  tree.length <- nrow(states)
  mrca.matrix[tree.ind:(tree.ind+tree.length-1),tree.ind:(tree.ind+tree.length-1)] <- mrca(tree,full=T) + tree.ind - 1
  tree.ind <- tree.ind + tree.length
  tree.sizes <- c(tree.sizes,length(b.lens))
  tree.sizes.nodes <- c(tree.sizes.nodes,nrow(states))
  node.heights.tree <- max(node.depth.edgelength(tree))-node.depth.edgelength(tree)
  node.heights <- c(node.heights,node.heights.tree)
}

#node.heights <- node.heights[all.parent]

#get branch lengths leading to specific node
child.branch.lengths <- rep(0,nrow(mrca.matrix))
for (i in 1:length(all.child)) {
  child.branch.lengths[all.child[i]] <- all.b.lens[i]
}

#ancestral branches of each possible mrca; sum(child.branch.lengths*ancestors) gives the distance from root to some mrca node
ancestors <- matrix(0,nrow(mrca.matrix),ncol(mrca.matrix))
for (i in 1:nrow(mrca.matrix)) {
  for (j in unique(mrca.matrix[i,])) {
    ancestors[i,j] = 1
  }
}

lens.to.root <- ancestors*matrix(child.branch.lengths,nrow=nrow(ancestors),ncol=length(child.branch.lengths),byrow=T)
lens.to.root <- lens.to.root[all.child,all.child]

languages <- read.csv('asjp/languages.csv')
languages$longname <- paste(languages$classification_wals,languages$Name,sep='.')
languages$longname <- gsub('-','_',languages$longname)
lonlat <- merge(data.df,languages,by='longname')[,c('longname','Longitude','Latitude')]
#rownames(lonlat) <- lonlat$longname
lon <- lonlat$Longitude
names(lon) <- lonlat$longname
lon <- lon[rownames(all.states)]
names(lon) <- rownames(all.states)

lat <- lonlat$Latitude
names(lat) <- lonlat$longname
lat <- lat[rownames(all.states)]

geo.missing <- ifelse(is.na(lat),1,0)

D <- ncol(all.states)/2
B <- length(all.b.lens)
T <- length(roots)
M = sum(geo.missing) - T

N <- length(geo.missing)
P <- N - (M + T)

mrca11 <- matrix(0,nrow=M,ncol=M)
mrca22 <- matrix(0,nrow=P,ncol=P)
mrca12 <- matrix(0,nrow=M,ncol=P)

n.present <- c()
n.absent <- c()

#don't get rid of NAs yet

i = 0
j = 1
k = 1
for (t in 1:length(tree.sizes.nodes)) {
  S = tree.sizes.nodes[t]
  root.ind <-i+((S+1)/2)+1
  missing.inds <- i + which(is.na(lon[(i+1):(i+S)]))
  missing.inds <- missing.inds[which(missing.inds!=root.ind)]
  present.inds <- i + which(!is.na(lon[(i+1):(i+S)]))
  X = length(missing.inds)
  Y = length(present.inds)
  Sigma11 = matrix(nrow=X,ncol=X)
  for (x in 1:X) {
    for (y in 1:X) {
      Sigma11[x,y] <- mrca.matrix[missing.inds[x],missing.inds[y]]
    }
  }
  Sigma22 = matrix(nrow=Y,ncol=Y)
  for (x in 1:Y) {
    for (y in 1:Y) {
      Sigma22[x,y] <- mrca.matrix[present.inds[x],present.inds[y]]
    }
  }
  Sigma12 = matrix(nrow=X,ncol=Y)
  for (x in 1:X) {
    for (y in 1:Y) {
      Sigma12[x,y] <- mrca.matrix[missing.inds[x],present.inds[y]]
    }
  }
  mrca11[j:(j+X-1),j:(j+X-1)] <- Sigma11
  mrca22[k:(k+Y-1),k:(k+Y-1)] <- Sigma22
  mrca12[j:(j+X-1),k:(k+Y-1)] <- Sigma12
  n.present <- c(n.present,Y)
  n.absent <- c(n.absent,X)
  i <- i + S
  j <- j + X
  k <- k + Y
}

lon <- lon[!is.na(lon)]
lat <- lat[!is.na(lat)]

F <- max(as.numeric(as.factor(fam.ID)))

rownames(all.states) <- NULL
colnames(all.states) <- NULL

data.list <- list(N=N,
                  P=P,
                  M=M,
                  T=T,
                  B=B,
                  F=F,
                  D=D,
                  n_present=n.present,
                  n_missing=n.absent,
                  child=all.child,
                  parent=all.parent,
                  brlen=all.b.lens,
                  fam_ID=as.numeric(as.factor(fam.ID)),
                  tiplik=all.states,
                  lon=lon,
                  lat=lat,
                  ancestor_lens=lens.to.root,
                  mrca22=mrca22,
                  roots=roots
)

fit <- stan(file=paste('models/',model.path,sep=''),data=data.list,include=TRUE)

save.image(file=paste('output/fit_',model.path,'_',curr.seed,'.Rdata',sep=''))
