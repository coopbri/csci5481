# Run with:
# Rscript hw3-plot-edges.r edges.txt tip-labels.txt
# or 
# Rscript hw3-plot-edges.r edges.txt tip-labels.txt bootstrap.txt
library('ape')
library('RColorBrewer')



args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2 || length(args) > 3) stop('Need 2 or 3 args.')

make.newick <- function(edge,edge.length,tip.label,root.ix=NULL){
	# assumes tips are the first indices (i.e. 1, 2, 3, ...)
	# and that those indices correspond to tip.label string vector
	# If root is NULL, then this is the initial call to the function.
	# Find the root node as the internal node that is never a child
	# Otherwise, this is a recursive call; print postorder traversal 
	# of subtree starting with root.ix.

	terminal.semicolon <- FALSE
	if(is.null(root.ix)){
		terminal.semicolon <- TRUE
		root.ix <- setdiff(edge[,1], edge[,2])
	}
	if(!any(edge[,1] == root.ix)){
		# Base case: If root is a tip, print its name
		row.ix <- which(edge[,2] == root.ix)
		ret.val <- paste(tip.label[root.ix],':',edge.length[row.ix],sep='')
	} else {
		# Recursion: get trees for children, then add this node
		# like this: (child1tree,child2tree,...):edgelength.to.parent
		ret.val <- '('
		for(row.ix in which(edge[,1] == root.ix)){
			subtree.val <- make.newick(edge, edge.length, tip.label, root.ix=edge[row.ix,2])
			delimiter <- ',' 
			if(ret.val == '(') delimiter <- '' # only add delimiter if not first child
			ret.val <- paste(ret.val,delimiter,subtree.val,sep='')
		}
		ret.val <- paste(ret.val,')',sep='')

		# test whether this is a child node and has a length to parent
		if(any(edge[,2] == root.ix)){
			ret.val <- paste(ret.val,':',edge.length[which(edge[,2] == root.ix)],sep='')
		}
	}
	if(terminal.semicolon) ret.val <- paste(ret.val,';',sep='')
	return(ret.val)
}

make.phylo.from.treelike.list <- function(t){
	# used to convert tree-like list to actual phylo object
	# uses intermediate Newick format because I couldn't 
	# manage to cast a list directly as a phylo object
	t <- read.tree(text=make.newick(t$edge,t$edge.length,t$tip.label))
	return(t)
}



# load edges
edges <- as.matrix(read.table(args[1],sep='\t',head=F,check=F,comment=''))
cat('\nLoaded edge matrix:\n')
print(edges)

# load tip labels and colors
tip.labels <- read.table(args[2], sep='\t',head=F,row=1,check=F,comment='')
cat('\nLoaded tip label table:\n')
print(tip.labels)

cols <- as.character(tip.labels[,2])
tip.labels <- rownames(tip.labels)
cat('\nTip colors are:\n')
print(cols)
cat('\nTip labels are:\n')
print(tip.labels)

# construct new tree
newtree <- list(edge=as.matrix(edges[,1:2]),
                edge.length=edges[,3],
                tip.label=tip.labels,
                Nnode=as.integer(max(edges[,1]) - min(edges[,1]) + 1))

class(newtree) <- 'phylo'

newtree <- make.phylo.from.treelike.list(newtree)

cat('\nLoaded tree:\n')
print(newtree)

cat('\nStructure of tree object:\n')
print(str(newtree))

# bootstrap support
boots <- NULL
if(length(args) == 3) boots <- scan(args[3])
cat('Bootstrap values are:\n')
print(boots)

# plot tree with labeled nodes (boostrap) and tips (phylum)
cat('Plotting tree...\n')
pdf('tree.pdf',width=8,height=8)
par(oma=c(0,0,0,0), mar=c(1,1,1,1))
plot.phylo(newtree,show.tip.label=FALSE,type='fan')
if(!is.null(boots)){
    nodelabels(col='black',frame='none',pie=boots,width=10,height=10,cex=.3)
}
tiplabels(tip.labels,col=cols, frame='none', cex=.5)
dev.off()

