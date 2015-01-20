# social network analysis: manipulation, visualization and model fitting 

install.packages("igraph")
install.packages('latentnet')
install.packages('statnet')
library(statnet,quietly=T)
require(igraph)
require(latentnet)

################# 
#### statnet ####
#################
network <- read.table("http://moreno.ss.uci.edu/ckm.dat",skip=9)
attributes <- read.table('http://moreno.ss.uci.edu/attributes.dat',skip=18)
network <- as.matrix(network)
dim(network)

net.a <- network(network[1:246,],directed=TRUE) # advise
net.d <- network(network[247:492,],directed=TRUE) # discussion
net.f <- network(network[493:738,],directed=TRUE) # friend

# some summary
summary(net.a)
network.dyadcount(net.a) # choose(246,2)*2
network.edgecount(net.a)
network.size(net.a)

# visualization
plot(net.a,displaylabels=F) # nice graph!
gplot(net.a)

plot(net.a,displaylabels=F,mode='circle') # circle arrangement
?plot

# add vertex covariates
net.a %v% 'city' <- attributes[,1]
net.a %v% 'adoption_date' <- attributes[,2]
net.a %v% 'med_sch_yr' <- attributes[,3]
net.a %v% 'meetings' <- attributes[,4]
net.a %v% 'jours' <- attributes[,5]
net.a %v% 'free_time' <- attributes[,6]
net.a %v% 'discuss' <- attributes[,7]
net.a %v% 'clubs' <- attributes[,8]
net.a %v% 'friends' <- attributes[,9]
net.a %v% 'community' <- attributes[,10]
net.a %v% 'patients' <- attributes[,11]
net.a %v% 'proximity' <- attributes[,12]
net.a %v% 'specialty' <- attributes[,13]

net.d %v% 'city' <- attributes[,1]
net.d %v% 'adoption_date' <- attributes[,2]
net.d %v% 'med_sch_yr' <- attributes[,3]
net.d %v% 'meetings' <- attributes[,4]
net.d %v% 'jours' <- attributes[,5]
net.d %v% 'free_time' <- attributes[,6]
net.d %v% 'discuss' <- attributes[,7]
net.d %v% 'clubs' <- attributes[,8]
net.d %v% 'friends' <- attributes[,9]
net.d %v% 'community' <- attributes[,10]
net.d %v% 'patients' <- attributes[,11]
net.d %v% 'proximity' <- attributes[,12]
net.d %v% 'specialty' <- attributes[,13]

net.f %v% 'city' <- attributes[,1]
net.f %v% 'adoption_date' <- attributes[,2]
net.f %v% 'med_sch_yr' <- attributes[,3]
net.f %v% 'meetings' <- attributes[,4]
net.f %v% 'jours' <- attributes[,5]
net.f %v% 'free_time' <- attributes[,6]
net.f %v% 'discuss' <- attributes[,7]
net.f %v% 'clubs' <- attributes[,8]
net.f %v% 'friends' <- attributes[,9]
net.f %v% 'community' <- attributes[,10]
net.f %v% 'patients' <- attributes[,11]
net.f %v% 'proximity' <- attributes[,12]
net.f %v% 'specialty' <- attributes[,13]

# descriptive statistics
summary(net.a)

table(component.dist(net.f)$csize) # component distribution
?component.dist

net.a.1 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==1))
class(net.a.1)
net.a.2 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==2))
net.a.3 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==3))
net.a.4 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==4))

table(net.a%v%'city')
table(net.a%v%'adoption_date')

mixingmatrix(net.a,'city')
mixingmatrix(net.a,'med_sch_yr')

# ergm models
ergm.1 <- ergm(net.a.1 ~ edges); summary(ergm.1)
names(ergm.1) # check model info

ergm.2 <- ergm(net.a.1 ~ edges + nodematch('discuss'));summary(ergm.2)

ergm.3 <- ergm(net.a ~ edges + nodematch("city") + nodematch("med_sch_yr") + nodematch('meetings') +
                 nodematch('jours') + nodematch('free_time') + nodematch('discuss') +
                 nodematch('clubs') + nodematch('friends') + nodematch('community') +
                 nodematch('patients') + nodematch('proximity') + 
                 nodematch('specialty'));summary(ergm.3)


# latent space models
lsm.1 <- ergmm(net.a.1 ~ euclidean(d = 2), tofit = c("mle"));summary(lsm.1)
m <- lsm.1$mle$Z
dim(m)
plot(m[,1],m[,2])
dim(network)


################
#### igraph ####
################

g <- graph.formula(1-2, 1-3, 2-3, 2-4, 3-5, 4-5, 4-6, 4-7, 5-6, 6-7) # create a undirected graph with V=7 and E specified
g$name <- "Toy Graph" # give g a name
V(g); E(g) # vertex and edge sequence
str(g)
plot(g)

dg <- graph.formula(1-+2,1-+3,2++3) # create a directed graph
str(dg)
plot(dg)

get.adjlist(g) # return adjacency list of the network
get.edgelist(g) # get edgelist of the network
get.adjacency(g) # get adjacency matrix of the network


h <- induced.subgraph(g, 1:5) # induced subgraph
h <- g - vertices(c(6,7)) # remove nodes 6 and 7 from the graph g

h <- h + vertices(c(6,7)) # recover nodes 6 and 7 back
g <- h + edges(c(4,6),c(4,7),c(5,6),c(6,7)) # recover edges associated with nodes 6 and 7

h1 <- h
h2 <- graph.formula(4-6, 4-7, 5-6, 6-7)
g1 <- graph.union(h1,h2) # set operation: union
g2 <- graph.intersection(h1,h2) # intersection
g3 <- graph.difference(h1,h2) # difference

V(dg)$name <- c("Sam", "Mary", "Tom") # adding attributes to vertices in graph dg
V(dg)$gender <- c('M','F','M') 

wg <- g
E(wg)$weight <- runif(ecount(wg)) # adding weight as attributes of edges
is.weighted(wg)
E(wg)$weight

install.packages('sand') # read in external data frame
library(sand)
g.lazega <- graph.data.frame(elist.lazega,directed="FALSE",vertices=v.attr.lazega)
g.lazega$name <- "Lazega Lawyers"

vcount(g.lazega) 
ecount(g.lazega)
str(g.lazega)
list.vertex.attributes(g.lazega)
list.edge.attributes(g.lazega)

is.simple(g) # simple graph: no loop, single edge multi-graph: not simple
neighbors(g,5) # neighbors of node 5
degree(dg, mode="in") # in-degree
degree(dg, mode="out") # out-degree
is.connected(g)
clusters(g) # components

diameter(g, weights=NA) # length of longest path in graph g

# network visualization

g.1 <- graph.lattice(c(5,5,5))
data(aidsblog)
summary(aidsblog)

igraph.options(vertex.size=3,vertex.label=NA,edge.arrow.size=.5)
par(mfrow=c(1,2))

plot(g.1,layout=layout.circle) # circle layout
title("5x5x5 Lattice")
plot(aidsblog,layout=layout.circle)
title("Blog Network")

plot(g.1, layout=layout.fruchterman.reingold) # spring-embedded method
title("5x5x5 Lattice")
plot(aidsblog, layout=layout.fruchterman.reingold)
title("Blog Network")

g.tree <- graph.formula(1-+2,1-+3,1-+4,2-+5,2-+6,2-+7,3-+8,3-+9,4-+10)
par(mfrow=c(1, 3))
igraph.options(vertex.size=30, edge.arrow.size=0.5,vertex.label=NULL)
plot(g.tree, layout=layout.circle)
plot(g.tree, layout=layout.reingold.tilford(g.tree,circular=T))
plot(g.tree, layout=layout.reingold.tilford)

############################# statnet

install.packages('statnet')
require (statnet)
require (ergm)

data(package='ergm') # check the dataset in ergm

data('faux.magnolia.high') # vertex attributes: name,race,sex,grade
model.covariates <- ergm(faux.mesa.high ~ edges + nodematch(('Grade'),diff=TRUE) + nodefactor('Sex'))
summary(model.covariates)

fmh <- faux.magnolia.high
fmh
summary(fmh)

plot(fmh,displayisolates=FALSE,vertex.cex=.7)

table(component.dist(fmh)$csize)

fmh.degreedist <- table(degree(fmh,cmode='indegree'))
fmh.degreedist

summary(fmh ~ degree(0:8,'Sex'))
summary(fmh ~ edges + triangle)

mixingmatrix(fmh,'Grade') # edges between and within groups defined by covariates
mixingmatrix(fmh,'Sex')

gr <- fmh %v% 'Grade' # a vector of attribute values for all actors
table(gr)

# ERGMs

m1 <- ergm(fmh ~ edges) # Bernoulli model
summary(m1) # interpretation of coefficient: logit(p)=log(p/1-p)=theta * delta statistics given the rest of the network

m2 <- ergm(fmh ~ edges + nodematch('Grade') + nodematch('Race') + nodematch('Sex'))
summary(m2) # completely heterogeneous -> completely homogeneous

sim2 <- simulate(m2,burnin=1e+6,verbose=TRUE,seed=123)
mixingmatrix(sim2,'Race')

m3 <- ergm(fmh ~ edges + triangle, seed=124) # model degeneracy: no convergency, bad fit

2^choose(9,2)/1e+9 # 68 billion

m4 <- ergm(fmh ~ edges + nodematch('Grade') + nodematch('Race') + nodematch('Sex') +
                  gwesp(0,fixed=TRUE),MCMCsamplesize=1e+5,maxit=15,verbose=TRUE,
                control=control.ergm(steplength=.25),seed=125)


####################### physician's new drug adoption data
library(igraph)

# read adjacency matrix and covariates
graph <- read.table("http://moreno.ss.uci.edu/ckm.dat",skip=9)
attributes <- read.table('http://moreno.ss.uci.edu/attributes.dat',skip=18)
graph <- as.matrix(graph)

head(attributes)
dim(attributes)


# divide by three different relations
graph.advise <- graph[1:246,]
graph.discussion <- graph[247:492,]
graph.friend <- graph[493:738,]

# create graph subjects
g.advise <- graph.adjacency(graph.advise)
g.discussion <- graph.adjacency(graph.discussion)
g.friend <- graph.adjacency(graph.friend)
str(g.advise)

# add in node-level covariates 
V(g.advise)$city <- attributes[,1]
V(g.advise)$adoption_date <- attributes[,2]
V(g.advise)$med_sch_yr <- attributes[,3]
V(g.advise)$meetings <- attributes[,4]
V(g.advise)$jours <- attributes[,5]
V(g.advise)$free_time <- attributes[,6]
V(g.advise)$discuss <- attributes[,7]
V(g.advise)$clubs <- attributes[,8]
V(g.advise)$friends <- attributes[,9]
V(g.advise)$community <- attributes[,10]
V(g.advise)$patients <- attributes[,11]
V(g.advise)$proximity <- attributes[,12]
V(g.advise)$specialty <- attributes[,13]

plot.igraph(g.advise,vertex.size=5,layout=layout.fruchterman.reingold)




