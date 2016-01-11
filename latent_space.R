
#################################################################
#  Identification and Bias-Amplification:                       #  
#  Latnet Space Approach to Contagion on Observational Networks #
#################################################################

# Search for (Ctrl + F) the following key words 
# to jump to the cooresponding section

# 1. read data
# 2. descriptive statistics
# 3. missing data
# 4. contagion variable
# 5. ERGMs and latent space models
# 6. survival analysis
# 7. simulation

install.packages("igraph")
install.packages('latentnet')
install.packages('statnet')
install.packages('Rtools')
install.packages('devtools')
install.packages('ggplot2')
install.packages('survival')
install.packages('stargazer')
install.packages('PCIT')
install.packages('ineq')
install.packages('knitr')

library(MASS)
library(statnet)
#library(network)
library(latentnet)
#require(igraph)
#library(sna)
library(ggplot2)
#library(devtools)
library(survival)
#library(PCIT)
library(stargazer)
#library(ineq)
#library(knitr)

detach("package:igraph", unload=TRUE)

#############################################
#### read data and creat netowrk objects ####
#############################################

network <- read.table("http://moreno.ss.uci.edu/ckm.dat",skip=9)
attributes <- read.table('http://moreno.ss.uci.edu/attributes.dat',skip=18)
# network <- as.matrix(network)
# dim(network)

net.a <- network(network[1:246,],directed=TRUE) # advisor
net.d <- network(network[247:492,],directed=TRUE) # discussion
net.f <- network(network[493:738,],directed=TRUE) # friend

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

net.a %v% 'missing' <- attributes$missing

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

net.d %v% 'missing' <- attributes$missing

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

net.f %v% 'missing' <- attributes$missing

dim(attributes)

# split into cities
net.a.1 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==1))
net.a.2 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==2))
net.a.3 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==3))
net.a.4 <- get.inducedSubgraph(net.a,which(net.a%v%'city'==4))

net.d.1 <- get.inducedSubgraph(net.d,which(net.a%v%'city'==1))
net.d.2 <- get.inducedSubgraph(net.d,which(net.a%v%'city'==2))
net.d.3 <- get.inducedSubgraph(net.d,which(net.a%v%'city'==3))
net.d.4 <- get.inducedSubgraph(net.d,which(net.a%v%'city'==4))

net.f.1 <- get.inducedSubgraph(net.f,which(net.a%v%'city'==1))
net.f.2 <- get.inducedSubgraph(net.f,which(net.a%v%'city'==2))
net.f.3 <- get.inducedSubgraph(net.f,which(net.a%v%'city'==3))
net.f.4 <- get.inducedSubgraph(net.f,which(net.a%v%'city'==4))

table(net.a%v%'city')
table(net.a%v%'adoption_date')

mixingmatrix(net.a,'city') # it mmeans no cross-city connections 
mixingmatrix(net.d,'city')
mixingmatrix(net.f,'city')

mixingmatrix(net.a,'specialty')


# attributes (covariates)
names(attributes) <- c('city','adoption_date','med_sch_yr','meetings','jours','free_time',
                       'discuss','clubs','friends','community','patients','proximity','specialty')
attributes$id <- 1:246

attributes$med_sch_yr[which(attributes$med_sch_yr==9)] <- NA
attributes$meetings[which(attributes$meetings==9)] <- NA
attributes$jours[which(attributes$jours==9)] <- NA
attributes$free_time[which(attributes$free_time==9)] <- NA
attributes$discuss[which(attributes$discuss==9)] <- NA
attributes$clubs[which(attributes$clubs==9)] <- NA
attributes$friends[which(attributes$friends==9)] <- NA
attributes$community[which(attributes$community==9)] <- NA
attributes$patients[which(attributes$patients==9)] <- NA
attributes$proximity[which(attributes$proximity==9)] <- NA
attributes$specialty[which(attributes$specialty==9)] <- NA
# head(attributes,30)

table(attributes$adoption_date)
attributes$adoption_date[which(attributes$adoption_date==98)] <- NA
attributes$adoption_date[which(attributes$adoption_date==99)] <- NA

# attributes$adoption_date <- as.numeric(attributes$adoption_date)
# attributes$event <- as.numeric(attributes$event)

table(attributes$adoption_date,exclude=NULL)

attributes.1 <- attributes[which(attributes$city==1),]
attributes.2 <- attributes[which(attributes$city==2),]
attributes.3 <- attributes[which(attributes$city==3),]
attributes.4 <- attributes[which(attributes$city==4),]

################################
#### descriptive statistics ####
################################

summary(net.a)
network.dyadcount(net.a) # choose(246,2)*2
network.edgecount(net.a)
network.size(net.a)
as.sociomatrix(net.a)
list.vertex.attributes(net.a)
list.edge.attributes(net.a)

sum(degree(net.a,cmode='indegree'))
table(degree(net.a,cmode='outdegree'))
table(degree(net.d,cmode='outdegree'))
table(degree(net.f,cmode='outdegree'))

g.a <- degree(net.a,cmode='indegree')
g.d <- degree(net.d,cmode='indegree')
g.f <- degree(net.f,cmode='indegree')

table(degree(net.d,cmode='indegree'))
table(degree(net.f,cmode='indegree'))

id <- 1:246
table(degree(net.a,c='indegree'))
out.degree.net.a <- degree(net.a,cmode='outdegree')
in.degree.net.a <- degree(net.a,cmode='indegree')
degree.net.a <- cbind(id,in.degree.net.a,out.degree.net.a)

out.degree.net.d <- degree(net.d,cmode='outdegree')
in.degree.net.d <- degree(net.d,cmode='indegree')
degree.net.d <- cbind(id,in.degree.net.d,out.degree.net.d)

out.degree.net.f <- degree(net.f,cmode='outdegree')
in.degree.net.f <- degree(net.f,cmode='indegree')
degree.net.f <- cbind(id,in.degree.net.f,out.degree.net.f)

# centrality measure (degree,betweenness,closeness)
d.a.1 <- degree(net.f.1,cmode='indegree')
d.a.2 <- degree(net.f.2,cmode='indegree')
d.a.3 <- degree(net.f.3,cmode='indegree')
d.a.4 <- degree(net.f.4,cmode='indegree')

b.a.1 <- betweenness(net.f.1)
b.a.2 <- betweenness(net.f.2)
b.a.3 <- betweenness(net.f.3)
b.a.4 <- betweenness(net.f.4)

c.f.1 <- closeness(net.f.1,cmode="suminvdir")
c.f.2 <- closeness(net.f.2,cmode="suminvdir")
c.f.3 <- closeness(net.f.3,cmode="suminvdir")
c.f.4 <- closeness(net.f.4,cmode="suminvdir")

cor(d.a.1,c.a.1)
cor(d.a.2,c.a.2)
cor(d.a.3,c.a.3)
cor(d.a.4,c.a.4)

c.a <- c(c.a.1,c.a.2,c.a.3,c.a.4)
c.d <- c(c.d.1,c.d.2,c.d.3,c.d.4)
c.f <- c(c.f.1,c.f.2,c.f.3,c.f.4)

attributes$a_closeness <- c.a
attributes$d_closeness <- c.d
attributes$f_closeness <- c.f


#components(net.a.1,connected='weak')

component.dist(net.a.1,connected='weak')$csize
component.dist(net.a.2,connected='weak')$csize
component.dist(net.a.3,connected='weak')$csize
component.dist(net.a.4,connected='weak')$csize

component.dist(net.d.1,connected='weak')$csize
component.dist(net.d.2,connected='weak')$csize
component.dist(net.d.3,connected='weak')$csize
component.dist(net.d.4,connected='weak')$csize

component.dist(net.f.1,connected='weak')$csize
component.dist(net.f.2,connected='weak')$csize
component.dist(net.f.3,connected='weak')$csize
component.dist(net.f.4,connected='weak')$csize

#component.largest(net.a.1,connected='weak')

summary(net.f.4)

# random graph
r.a.1 <- rgraph(106,tprob=0.01783083)
r.a.2 <- rgraph(41,tprob=0.03836735)
r.a.3 <- rgraph(35,tprob=0.03699789)
r.a.4 <- rgraph(33,tprob=0.06218487) 

r.d.1 <- rgraph(113,tprob=0.02040967)
r.d.2 <- rgraph(44,tprob=0.04571429)
r.d.3 <- rgraph(41,tprob=0.04756871)
r.d.4 <- rgraph(33,tprob=0.07226891)

r.f.1 <- rgraph(110,tprob=0.01679929)
r.f.2 <- rgraph(46,tprob=0.04489796)
r.f.3 <- rgraph(43,tprob=0.04334038)
r.f.4 <- rgraph(34,tprob=0.07226891)

# average path length

mean(geodist(net.f.1,inf.replace=NA)$gdist,na.rm=T);
mean(geodist(net.f.2,inf.replace=NA)$gdist,na.rm=T);
mean(geodist(net.f.3,inf.replace=NA)$gdist,na.rm=T);
mean(geodist(net.f.4,inf.replace=NA)$gdist,na.rm=T)

mean(geodist(r.f.1,inf.replace=NA)$gdist,na.rm=T);mean(geodist(r.f.2,inf.replace=NA)$gdist,na.rm=T);mean(geodist(r.f.3,inf.replace=NA)$gdist,na.rm=T);mean(geodist(r.f.4,inf.replace=NA)$gdist,na.rm=T)

mean(clusteringCoefficient(as.matrix(r.f.1)),na.rm=T);
mean(clusteringCoefficient(as.matrix(r.f.2)),na.rm=T);
mean(clusteringCoefficient(as.matrix(r.f.3)),na.rm=T);
mean(clusteringCoefficient(as.matrix(r.f.4)),na.rm=T)

mean(geodist(r.a.1,inf.replace=NA)$gdist,na.rm=T)
mean(clusteringCoefficient(r.a.1),na.rm=T)

gtrans(net.f.1,measure='weakcensus')
gtrans(net.f.2,measure='weakcensus')
gtrans(net.f.3,measure='weakcensus')
gtrans(net.f.4,measure='weakcensus')

gtrans(r.f.1,measure='weakcensus')
gtrans(r.a.3,measure='weakcensus')
gtrans(r.f.3,measure='weakcensus')
gtrans(r.f.4,measure='weakcensus')

kcycle.census(net.f.1,maxlen=3)
kcycle.census(net.f.2,maxlen=3)
kcycle.census(net.f.3,maxlen=3)
kcycle.census(net.f.4,maxlen=3)

kcycle.census(r.f.1,maxlen=3)
kcycle.census(r.f.2,maxlen=3)
kcycle.census(r.f.3,maxlen=3)
kcycle.census(r.f.4,maxlen=3)

# average clustering coefficient

id1 <- 1:117 #117
id2 <- 118:167 # 50
id3 <- 168:211 # 44
id4 <- 212:246 # 35

plot(density(net.a%v%'med_sch_yr'))
table(net.a.1%v%'med_sch_yr')

table(component.dist(net.f)$csize) # component distribution

# visualization
plot(net.f.4) 
plot(net.d,displaylabels=F)
plot(net.f,displaylabels=F)

# x <- m$mle$Z[,1]
# y <- m$mle$Z[,2]
# t <- net.a.1%v%'adoption_date'
# length(x)
# d <- data.frame(x,y,t)
# qplot(x,y,data=d,color=t)
# plot(net.a,displaylabels=F,mode='circle') # circle arrangement


#######################
#### missing data ####
#######################

missing_patient <- as.numeric(is.na(attributes$patient))
missing_med_sch_yr <- as.numeric(is.na(attributes$med_sch_yr))
missing_meetings <- as.numeric(is.na(attributes$meetings))

attributes$id[which(attributes$specialty<=3)]

tapply(data$missing,data$specialty,sum)

data$missing_patient <- missing_patient
data$missing_med_sch_yr <- as.numeric(is.na(attributes$med_sch_yr))
data$missing_meetings <- as.numeric(is.na(attributes$meetings))

head(attributes)
table(attributes$adoption_date,exclude=NULL)
table(attributes$missing)

attributes$missing <- 1:nrow(attributes)

for (i in 1:nrow(attributes)){
  if(attributes$adoption_date[i]%in%NA){
    attributes$missing[i] <- 1     
  } else {
    attributes$missing[i] <- 0
  }
}  

c <- net.a.1%v%'missing'
c <- net.a.2%v%'missing'
c <- net.a.3%v%'missing'
c <- net.a.4%v%'missing'

c <- net.d.1%v%'missing'
c <- net.d.2%v%'missing'
c <- net.d.3%v%'missing'
c <- net.d.4%v%'missing'

c <- net.f.1%v%'missing'
c <- net.f.2%v%'missing'
c <- net.f.3%v%'missing'
c <- net.f.4%v%'missing'

for (i in 1:length(c)){
  if(c[i]==1){
    c[i] <- 'white'
  } else {
    c[i] <- 'red'
  }
}

gplot(net.a.1,vertex.col=c)
gplot(net.a.2,vertex.col=c)
gplot(net.a.3,vertex.col=c)
gplot(net.a.4,vertex.col=c)

gplot(net.d.1,vertex.col=c)
gplot(net.d.2,vertex.col=c)
gplot(net.d.3,vertex.col=c)
gplot(net.d.4,vertex.col=c)

gplot(net.f.1,vertex.col=c)
gplot(net.f.2,vertex.col=c)
gplot(net.f.3,vertex.col=c)
gplot(net.f.4,vertex.col=c)


table(net.a.2%v%'missing')
d <- attributes$adoption_date
c <-d

for (i in 1:length(d)){
  if(d[i] %in% NA){
    c[i] <- 'white'
  } else if(d[i]>=1 & d[i]<=8){
    c[i] <- 'red'
  } else if (d[i]>=7 & d[i]<=18) {
    c[i] <- 'orange'
  } 
}

ta
ble(c)
table(attributes$adoption_date,exclude=NULL)
gplot(net.d.3,vertex.col=c)


net.a.4%v%'missing'
cor(attributes,use='complete.obs')
tapply(attributes$missing, attributes$city, sum)

model <-glm(missing ~ factor(city) + med_sch_yr + factor(meetings) + 
              factor(free_time) +
              factor(discuss) + factor(friends) + community + factor(specialty) + 
              log(jours) + factor(clubs) + factor(patients) + 
              a_in_degree + d_in_degree + f_in_degree + 
              a_closeness + d_closeness + f_closeness + 1, family="binomial", data=data)
summary(model)

model <-glm(missing_patient ~ med_sch_yr + factor(discuss) +
              log(jours) + patients + factor(specialty) + 
              factor(meetings) + factor(free_time) + friends +
              factor(clubs) + community + factor(proximity) + 1,
              family='binomial',data=data)
summary(model)

model <-glm(missing ~ med_sch_yr + 
              log(jours) + patients + factor(specialty) + 
              community + factor(proximity) + 1,
            family='binomial',data=data)

###########################
#### survival analysis ####
###########################

# numeric:med_sch_yr,jours,friends,community,patients,proximity
# categorical:meetings,free_time,discuss,clubs,specialty (~ group)

data <- attributes.1
data <- attributes.2
data <- attributes.3
data <- attributes.4

dim(data);tail(data)

dim(data.a)
dim(data.d)
dim(data.f)

head(data.a)
test <- data.a
test <- data.d
test <- data.f

m.3 <- lm(contagion ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
          factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
          log(jours+1) + factor(clubs) + factor(patients) +z1+z2+1,data=test)
summary(m.3)
stargazer(m,report='vc*')

summary(model.a.3)

model.a <- glm(event ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
                 factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
                 log(jours+1) + factor(clubs) + factor(patients) +
                 contagion + data.d$contagion + data.f$contagion +1, family="binomial", data=test)




data <- attributes
head(data)

dim(data)

data$event <- rep(1,nrow(data))
data$event[which(data$adoption_date==18)] <- 0

head(attributes)
dim(attributes)
dim(data)

# m <- network[1:117,1:117]
# 
# dim(data)
# x <- data[,c(3:13)]
# x <- apply(x,2,as.numeric)
# dim(x)

#nonparametric exploration
prop.table(table(data$adoption_date,exclude=NULL))
table(data$adoption_date)

table(data$city)

head(data)

kmsurvival <- survfit(Surv(data$adoption_date,data$event) ~ 1)
kmsurvival <- survfit(Surv(data$adoption_date,data$event) ~ data$city)

# head(data.a)
# coxsurvival <- coxph(Surv(data.a$period,data.a$event) ~ 
#                        as.matrix(data.a[,c(11,14:24)]))
# summary(coxsurvival)

stargazer(kmsurvival)
summary(kmsurvival)
names(kmsurvival)
class(kmsurvival)

d <- cbind(kmsurvival$time,kmsurvival$n.risk,kmsurvival$n.event,kmsurvival$surv,
      kmsurvival$std.err,kmsurvival$lower,kmsurvival$upper)
kmsurvival$std.err
head(d)

d <- as.data.frame(d)

apply(d,2,as.character)

stargazer(d)

stargazer(attributes)

plot(kmsurvival,xlab='time',ylab='survial rate')
par(mfrow=c(1,1))
plot(kmsurvival,conf.int='both',col=c('red','pink','orange','black'),
     xlab='time',ylab='survival rate')


x <- data[,c(3:13,16:21)]
head(x);dim(x)
x <- as.matrix(x)

exponential <- survreg(Surv(data$adoption_date,data$event)~ x,dist='exponential')

summary(exponential)
View(data)

# prop.table(table(attributes$adoption_date,exclude=NULL))

# head(m)
# table(rowSums(m))
# table(colSums(m))
# 
# matrix.a[1,110]==1
# 
# network.size(net.a.1)
# head(data)
# 
# get.neighborhood(net.a.1,10,'out')
# attributes$adoption_date[i]

install.packages('car')
library(car)
vif(model)
stargazer(model)
?stargazer
names(model)

a_in_degree <- degree(net.a,cmode='indegree')
d_in_degree <- degree(net.d,cmode='indegree')
f_in_degree <- degree(net.f,cmode='indegree')

cor(d_in_degree,f_in_degree)

attributes$a_in_degree <- a_in_degree
attributes$d_in_degree <- f_in_degree
attributes$f_in_degree <- d_in_degree

dim(data)
head(data)
head(data.a)

normalize <- function (x){
  x <- as.matrix(x)
  (x - mean(x,na.rm=T))/sd(x,na.rm=T)
}

data.a.nor <- data.a
colnames(data.a)

head(data.a.nor)

data.a.nor[,c(11,29,30)] <- normalize(data.a.nor[,c(11,29,30)]
)

model <- lm(log(adoption_date) ~ factor(med_sch_yr) + meetings + free_time +
              log(jours) + discuss + factor(clubs) + friends +
              community + factor(patients) + proximity + specialty +
              log(a_in_degree+1) + log(d_in_degree+1) + log(f_in_degree+1) +
              a_closeness + d_closeness + f_closeness + 1,data=data)
summary(model)

model <- lm(log(adoption_date) ~ med_sch_yr + 
              factor(clubs)  +
              community + patients  + 
              log(a_in_degree+1) + log(d_in_degree+1) + log(f_in_degree+1) +
              a_closeness + d_closeness + f_closeness + 1
            ,data=data)
summary(model)

stargazer(model,single.row=T)
?stargazer


dim(attributes)
l <- lm(adoption_date ~ 1+ factor(med_sch_yr) + factor(meetings) + log(jours) + factor(free_time) + 
          factor(discuss) + factor(clubs) + factor(friends) + factor(community) + 
          factor(patients) + factor(proximity) + factor(specialty), data=attributes)
l <- lm(adoption_date ~ 1+ med_sch_yr + meetings + log(jours) + free_time + 
          discuss + clubs + friends + community + 
          patients + proximity + specialty, data=attributes)
summary(l)

d <- as.data.frame(cbind(attributes.2$adoption_date,m.a.2$mle$Z))
l <- lm(V1~V2+V3,data=d)
summary(l)

# fit logistic model

data.a.1 <- read.csv('network.data.a.1.csv')
data.a.2 <- read.csv('network.data.a.2.csv')
data.a.3 <- read.csv('network.data.a.3.csv')
data.a.4 <- read.csv('network.data.a.4.csv')

data.d.1 <- read.csv('network.data.d.1.csv')
data.d.2 <- read.csv('network.data.d.2.csv')
data.d.3 <- read.csv('network.data.d.3.csv')
data.d.4 <- read.csv('network.data.d.4.csv')

data.f.1 <- read.csv('network.data.f.1.csv')
data.f.2 <- read.csv('network.data.f.2.csv')
data.f.3 <- read.csv('network.data.f.3.csv')
data.f.4 <- read.csv('network.data.f.4.csv')

dim(data.a.1)
head(data.d.2)

data.a <- as.data.frame(rbind(data.a.1,data.a.2,data.a.3,data.a.4))
data.d <- as.data.frame(rbind(data.d.1,data.d.2,data.d.3,data.d.4))
data.f <- as.data.frame(rbind(data.f.1,data.f.2,data.f.3,data.f.4))

data.a <- merge(data.a,degree.net.a,by='id')
data.d <- merge(data.d,degree.net.d,by='id')
data.f <- merge(data.f,degree.net.f,by='id')

dim(data.a.1)
dim(data.d.2)
dim(data.f.1)

head(data.a)

data.a[,c(11,29,30)] <- normalize(data.a[,c(11,29,30)])
data.d[,c(11,29,30)] <- normalize(data.d[,c(11,29,30)])
data.f[,c(11,29,30)] <- normalize(data.f[,c(11,29,30)])

# full model 
test <- data.a.1
test <- data.d
test <- data.f

head(data.all)
dim(data.all)
data.all = data.all[,-c(31,32)]

cor(data.all$z1.f,data.all$z1.d)


data.all$contagion.d <- data.d$contagion
data.all$contagion.f <- data.f$contagion

data.all$z1.d <- data.d$z1
data.all$z2.d <- data.d$z2

data.all$z1.f <- data.f$z1
data.all$z2.f <- data.f$z2


data.all <- as.data.frame(cbind(data.a,data.d$contagion,data.f$contagion))
head(data.all)
test <- data.all

head(data.a.1)
test <- as.data.frame(rbind(data.a.1,data.a.2,data.a.4))
test <- as.data.frame(rbind(data.d.1,data.d.2,data.d.4))
test <- as.data.frame(rbind(data.f.1,data.f.2,data.f.4))

dim(test)

m <- lm(contagion ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
          factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
          log(jours+1) + factor(clubs) + factor(patients) +
          z1 + z2 + 1,data=test)
summary(m)

m <- glm(contagion ~ factor(city)+z1+z2+1,family='binomial',data=test)
summary(m)

summary(model.a.3)

model.a <- glm(event ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
                 factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
                 log(jours+1) + factor(clubs) + factor(patients) +
                 contagion + data.d$contagion + data.f$contagion +1, family="binomial", data=test)

model <- glm(event ~ factor(med_sch_yr) + factor(meetings) + factor(free_time)+
               factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
               log(jours+1) + factor(clubs) + factor(patients) +1,
             family='binomial',data=test)

model <- glm(event ~ med_sch_yr + meetings + free_time +
               discuss + friends + community + specialty + 
               log(jours+1) + clubs + patients +1,
             family='binomial',data=test)
summary(model)


model.a <- glm(event ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
                 factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
                 log(jours+1) + factor(clubs) + factor(patients) +
                 contagion + data.d$contagion + data.f$contagion +
                 z1.d + z2.d + z1.f + z2.f + z1 + z2 + 1, family="binomial", data=test)

stargazer(model.a.2,model.a.3,single.row=T,report='vc*')


model.d <- glm(event ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
                 factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
                 log(jours+1) + factor(clubs) + factor(patients) +
                 contagion + z1 + z2 + 1, family="binomial", data=test)

model.f <- glm(event ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
                 factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
                 log(jours+1) + factor(clubs) + factor(patients) +
                 contagion + z1 + z2 + 1, family="binomial", data=test)

stargazer(model.a,model.d,model.f,single.row=T,report='vc*')
stargazer(model.a.2,model.a,single.row=T,report='vc*')

model <- glm(event ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
               factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
               log(jours+1) + factor(clubs) + factor(patients) +
               data.d$contagion + data.f$contagion +
               contagion + 1, family="binomial", data=test)

model <-glm(event ~ factor(city)+contagion+1,family='binomial',data=test)
summary(model)

model.z <-glm(event ~ factor(city) + factor(med_sch_yr) + factor(meetings) + factor(free_time)+
                factor(discuss) + factor(friends) + factor(community) + factor(specialty) + 
                log(jours+1) + factor(clubs) + factor(patients) +
                contagion + data.d$contagion + data.f$contagion +
                z1 + z2 + 1, family="binomial", data=test)

summary(model.z)

model <- glm(event ~ z1 + z2 + 1, family='binomial',data=test)
summary(model)

# selected model
model <- glm(event ~ factor(med_sch_yr) + log(jours) + factor(patients) + 
               contagion +1 , family="binomial", data=test)
summary(model)

model <- glm(event ~ factor(city) + factor(med_sch_yr) + log(jours) + factor(patients) + 
               contagion +z1 + z2 + 1 , family="binomial", data=test)
summary(model)

model <- glm(event ~ factor(city) + factor(med_sch_yr) + log(jours) + factor(patients) + 
               contagion + z1 +z2 +z3 +1 , family="binomial", data=test)
summary(model)

# table(attributes.1$adoption_date)
# table(attributes.2$adoption_date)
# table(attributes.3$adoption_date)
# table(attributes.4$adoption_date)

# attributes$event <- rep(1,246)
# attributes$event[which(attributes$adoption_date==19)] <- 0
# attributes$event[which(attributes$adoption_date==98)] <- NA


# for (i in 1:length(attributes.1)){
#   if(attributes.1$adoption_date[i]=='19'){
#     attributes.1$event[i] <- 0
#   } else if (attributes.1$adoption_date[i]==NA){
#     attributes.1$event[i] <- NA
#   } else {
#     attributes.1$event[i] <- 1
#   }  
# }

#########################################
#### creating the contagion variable ####
#########################################

# get neighbor ids
# specify id and size
id <- 1:117 #117
id <- 118:167 # 50
id <- 168:211 # 44
id <- 212:246 # 35

size <- 117
size <- 50
size <- 44
size <- 35

neighbor.1 <- rep(-1,size)
neighbor.2 <- rep(-1,size)
neighbor.3 <- rep(-1,size)

# get neighbors for an id (up to three neighbors) 

id
length(id)
neighbor.1 <- rep(NA,107)
neighbor.2 <- rep(NA,107)
neighbor.3 <- rep(NA,107)

for (i in 1:107){
  neighbor.1[i] <- get.neighborhood(net.a.1,id[i],'out')[1]
  neighbor.2[i] <- get.neighborhood(net.a.1,id[i],'out')[2]
  neighbor.3[i] <- get.neighborhood(net.a.1,id[i],'out')[3]
}

get.neighborhood(s)
length(neighbor.1)
neighbor.1[1:10]
neighbor.2[1:10]
neighbor.3[1:10]

rm(d)
d <- as.data.frame(cbind(id,neighbor.1,neighbor.2,neighbor.3))
d
d$id%in%id
dim(d);head(d)

# neighbor.1.adoption <- rep(-1,size)
# length(neighbor.1.adoption)
# 
# table(data$adoption_date,exclude=NULL)

# use ids to get their adoption date by checking the attributes table

neighbor.1.adoption <- attributes$adoption_date[neighbor.1]
neighbor.2.adoption <- attributes$adoption_date[neighbor.2]
neighbor.3.adoption <- attributes$adoption_date[neighbor.3]

# sanity check
neighbor.1.adoption[1:10]
neighbor.2.adoption[1:10]
neighbor.3.adoption[1:10]

attributes[which(attributes$id==5),'adoption_date']

d <- cbind(d,neighbor.1.adoption,neighbor.2.adoption,neighbor.3.adoption)
head(d,30)

# creat adoption dummy (before time =1, after/on time =0) 

size <- 107
time <- 1

neighbor.1.dummy <- rep(-1,size)
neighbor.1.dummy <- (neighbor.1.adoption - time)
neighbor.2.dummy <- rep(-1,size)
neighbor.2.dummy <- (neighbor.2.adoption - time)
neighbor.3.dummy <- rep(-1,size)
neighbor.3.dummy <- (neighbor.3.adoption - time)

for (i in 1:size){
  if (neighbor.1.dummy[i] %in% NA){
    neighbor.1.dummy[i] <- NA
  }else if (neighbor.1.dummy[i] < 0){
    neighbor.1.dummy[i] <- 1
  }else {neighbor.1.dummy[i] <- 0}
}

for (i in 1:size){
  if (neighbor.2.dummy[i] %in% NA){
    neighbor.2.dummy[i] <- NA
  }else if (neighbor.2.dummy[i] < 0){
    neighbor.2.dummy[i] <- 1
  }else {neighbor.2.dummy[i] <- 0}
}

for (i in 1:size){
  if (neighbor.3.dummy[i] %in% NA){
    neighbor.3.dummy[i] <- NA
  }else if (neighbor.3.dummy[i] < 0){
    neighbor.3.dummy[i] <- 1
  }else {neighbor.3.dummy[i] <- 0}
}

t <- as.data.frame(cbind(d,neighbor.1.dummy,neighbor.2.dummy,neighbor.3.dummy))
t

for (i in 1:size){
  t$contagion[i] <- mean(c(t$neighbor.1.dummy[i],t$neighbor.2.dummy[i],t$neighbor.3.dummy[i]),na.rm=T)
}

period <- rep(time,size)
t <- cbind(t,period)
all <- t
head(all)

for (time in 2:5){
  
  neighbor.1.dummy <- rep(-1,size)
  neighbor.1.dummy <- (neighbor.1.adoption - time)
  neighbor.2.dummy <- rep(-1,size)
  neighbor.2.dummy <- (neighbor.2.adoption - time)
  neighbor.3.dummy <- rep(-1,size)
  neighbor.3.dummy <- (neighbor.3.adoption - time)
  
  for (i in 1:size){
    if (neighbor.1.dummy[i] %in% NA){
      neighbor.1.dummy[i] <- NA
    }else if (neighbor.1.dummy[i] < 0){
      neighbor.1.dummy[i] <- 1
    }else {neighbor.1.dummy[i] <- 0}
  }
  
  for (i in 1:size){
    if (neighbor.2.dummy[i] %in% NA){
      neighbor.2.dummy[i] <- NA
    }else if (neighbor.2.dummy[i] < 0){
      neighbor.2.dummy[i] <- 1
    }else {neighbor.2.dummy[i] <- 0}
  }
  
  for (i in 1:size){
    if (neighbor.3.dummy[i] %in% NA){
      neighbor.3.dummy[i] <- NA
    }else if (neighbor.3.dummy[i] < 0){
      neighbor.3.dummy[i] <- 1
    }else {neighbor.3.dummy[i] <- 0}
  }
  
  t <- as.data.frame(cbind(d,neighbor.1.dummy,neighbor.2.dummy,neighbor.3.dummy))
  
  for (i in 1:size){
    t$contagion[i] <- mean(c(t$neighbor.1.dummy[i],t$neighbor.2.dummy[i],t$neighbor.3.dummy[i]),na.rm=T)
  }
  
  period <- rep(time,size)
  t <- cbind(t,period)
  all <- rbind(all,t)
  
}

dim(all);head(all);head(X)
X <- X[,-1]
X$id <- id
head(X)
data <- merge(all,X,by.x='id')
dim(data)

for (i in 1:nrow(data)){
  if(data$adoption[i] %in% NA){
    data$event[i] <- NA
  } else if(data$period[i] < data$adoption[i]){
    data$event[i] <- 0
  } else {
    data$event[i] <- 1
  }
}

# table(all$neighbor.3.dummy)
# all$neighbor.1.dummy[which(all$neighbor.1.dummy<0)] <- 0
# all$neighbor.2.dummy[which(all$neighbor.2.dummy<0)] <- 0
# all$neighbor.3.dummy[which(all$neighbor.3.dummy<0)] <- 0
# 
# all$contagion <- rep(-1,nrow(all))
# for (i in 1:nrow(all)){
# all$contagion[i] <- sum(all$neighbor.1.dummy[i],all$neibor.2.dummy[i],all$neighbor.3.dummy[i],na.rm=1)
# }
# head(all)
# table(all$contagion)

# creating event variable (0,1)s

all$event <- rep(-1,nrow(all))
head(all)

for (i in 1:nrow(all)){
  if(all$adoption_date[i] %in% NA){
    all$event[i] <- NA
  } else if(all$period[i] < all$adoption_date[i]){
    all$event[i] <- 0
  } else {
    all$event[i] <- 1
  }
}


head(data)
table(data$event)


all$event[which(all$adoption_date==18)] <- 0

# d$contagion <- (d$contagion - mean(d$contagion,na.rm=T))/sd(d$contagion,na.rm=T)
# head(all[which(all$id==125),])

# only keep non-missing adoption dates and get rid of obs when period > adoption date

test <- test[which(complete.cases(test$event)),]
test <- test[which(test$period<=test$adoption_date),]

test$contagion[which(test$neighbor.1.dummy%in%NA & test$neighbor.2.dummy%in%NA & test$neighbor.3.dummy%in%NA)] <- NA

table(test$contagion,exclude=NULL)
table(test$event,exclude=NULL)
table(test$adoption_date,exclude=NULL)

prop.table(table(test$adoption_date))

# sanity check

tapply(test$event,test$id,sum)
test[test$id==108,'adoption_date']

# head(test)
# test[10:20,c(1,8,9,10,11,27)]

# creating two new types of treatment: contagion dummy and contagion percent

for (i in 1:nrow(test)){
  test$contagion.dummy[i] <- sum(test$neighbor.1.dummy[i],test$neighbor.2.dummy[i],test$neighbor.3.dummy[i],na.rm=T)
}

test$contagion.percent <- test$contagion * 100
head(test)

write.csv(data.a.1,'network.data.a.1.csv',row.names=F)
write.csv(data.a.2,'network.data.a.2.csv',row.names=F)
write.csv(data.a.3,'network.data.a.3.csv',row.names=F)
write.csv(data.a.4,'network.data.a.4.csv',row.names=F)

write.csv(data.d.1,'network.data.d.1.csv',row.names=F)
write.csv(data.d.2,'network.data.d.2.csv',row.names=F)
write.csv(data.d.3,'network.data.d.3.csv',row.names=F)
write.csv(data.d.4,'network.data.d.4.csv',row.names=F)

write.csv(data.f.1,'network.data.f.1.csv',row.names=F)
write.csv(data.f.2,'network.data.f.2.csv',row.names=F)
write.csv(data.f.3,'network.data.f.3.csv',row.names=F)
write.csv(data.f.4,'network.data.f.4.csv',row.names=F)


####################
#### simulation ####
####################

x1 <- attributes.1$med_sch_yr
x2 <- attributes.1$jours
z1 <- m.a.1$mle$Z[,1]
z2 <- m.a.1$mle$Z[,2]
u <- rnorm(117)

id <- c(1:91,93:108)

alpha <- -10
beta.x1 <- .5
beta.x2 <- 1.5
beta.z1 <- 2
beta.u <- 2
  
# bias reduction!
alpha <- -10
beta.x1 <- 2
beta.x2 <- 2
beta.z1 <- .5
beta.u <- 2

alpha <- -10
beta.x1 <- 2
beta.x2 <- 2
beta.z1 <- 3
beta.u <- 5

rm(data)
rm(all)
rm(X)

adoption <- rep(0,length(x1))

beta <- cbind(alpha,beta.x1,beta.x2,beta.z1,beta.u)
X <- as.data.frame(cbind(rep(1,length(x1)),x1,x2,z1,u,adoption))
dim(X)
X <- X[complete.cases(X),]
head(X)

logodds <- beta%*%t(X[,1:5])
logodds <- logodds[!is.na(logodds)]

prob <- function(x){
  exp(x)/(1+exp(x))
}

p <- prob(logodds)
t1 <- rbinom(107,1,p)
t2 <- rbinom(107,1,p)
t3 <- rbinom(107,1,p)
t4 <- rbinom(107,1,p)

X$adoption[t1==1] <- 1
X$adoption[t1==0&t2==1] <- 2
X$adoption[t1==0&t2==0&t3==1] <- 3
X$adoption[t1==0&t2==0&t3==0&t4==1] <- 4

X$adoption[X$adoption==0] <- 5
table(X$adoption)

r <- rgraph(106,tprob=0.01783083)
get.neighborhood(net.a.1,id,'out')
s.adoption <- rep(0,106)

# get neighbors
neighbor.1.adoption <- X$adoption[neighbor.1]
neighbor.2.adoption <- X$adoption[neighbor.2]
neighbor.3.adoption <- X$adoption[neighbor.3]

# sanity check
neighbor.1.adoption[1:10]
neighbor.2.adoption[1:10]
neighbor.3.adoption[1:10]

d <- cbind(d,neighbor.1.adoption,neighbor.2.adoption,neighbor.3.adoption)
head(d,30)

# creat adoption dummy (before time =1, after/on time =0)
size <- 107
time <- 1

neighbor.1.dummy <- rep(-1,size)
neighbor.1.dummy <- (neighbor.1.adoption - time)
neighbor.2.dummy <- rep(-1,size)
neighbor.2.dummy <- (neighbor.2.adoption - time)
neighbor.3.dummy <- rep(-1,size)
neighbor.3.dummy <- (neighbor.3.adoption - time)

for (i in 1:size){
  if (neighbor.1.dummy[i] %in% NA){
    neighbor.1.dummy[i] <- NA
  }else if (neighbor.1.dummy[i] < 0){
    neighbor.1.dummy[i] <- 1
  }else {neighbor.1.dummy[i] <- 0}
}

for (i in 1:size){
  if (neighbor.2.dummy[i] %in% NA){
    neighbor.2.dummy[i] <- NA
  }else if (neighbor.2.dummy[i] < 0){
    neighbor.2.dummy[i] <- 1
  }else {neighbor.2.dummy[i] <- 0}
}

for (i in 1:size){
  if (neighbor.3.dummy[i] %in% NA){
    neighbor.3.dummy[i] <- NA
  }else if (neighbor.3.dummy[i] < 0){
    neighbor.3.dummy[i] <- 1
  }else {neighbor.3.dummy[i] <- 0}
}

t <- as.data.frame(cbind(d,neighbor.1.dummy,neighbor.2.dummy,neighbor.3.dummy))
t

for (i in 1:size){
  t$contagion[i] <- mean(c(t$neighbor.1.dummy[i],t$neighbor.2.dummy[i],t$neighbor.3.dummy[i]),na.rm=T)
}

period <- rep(time,size)
t <- cbind(t,period)
all <- t
head(all)

for (time in 2:5){
  
  neighbor.1.dummy <- rep(-1,size)
  neighbor.1.dummy <- (neighbor.1.adoption - time)
  neighbor.2.dummy <- rep(-1,size)
  neighbor.2.dummy <- (neighbor.2.adoption - time)
  neighbor.3.dummy <- rep(-1,size)
  neighbor.3.dummy <- (neighbor.3.adoption - time)
  
  for (i in 1:size){
    if (neighbor.1.dummy[i] %in% NA){
      neighbor.1.dummy[i] <- NA
    }else if (neighbor.1.dummy[i] < 0){
      neighbor.1.dummy[i] <- 1
    }else {neighbor.1.dummy[i] <- 0}
  }
  
  for (i in 1:size){
    if (neighbor.2.dummy[i] %in% NA){
      neighbor.2.dummy[i] <- NA
    }else if (neighbor.2.dummy[i] < 0){
      neighbor.2.dummy[i] <- 1
    }else {neighbor.2.dummy[i] <- 0}
  }
  
  for (i in 1:size){
    if (neighbor.3.dummy[i] %in% NA){
      neighbor.3.dummy[i] <- NA
    }else if (neighbor.3.dummy[i] < 0){
      neighbor.3.dummy[i] <- 1
    }else {neighbor.3.dummy[i] <- 0}
  }
  
  t <- as.data.frame(cbind(d,neighbor.1.dummy,neighbor.2.dummy,neighbor.3.dummy))
  
  for (i in 1:size){
    t$contagion[i] <- mean(c(t$neighbor.1.dummy[i],t$neighbor.2.dummy[i],t$neighbor.3.dummy[i]),na.rm=T)
  }
  
  period <- rep(time,size)
  t <- cbind(t,period)
  all <- rbind(all,t)
  
}

dim(all);head(all);head(X)
X <- X[,-1]
X$id <- id
head(X)
data <- merge(all,X,by.x='id')
dim(data)

for (i in 1:nrow(data)){
  if(data$adoption[i] %in% NA){
    data$event[i] <- NA
  } else if(data$period[i] < data$adoption[i]){
    data$event[i] <- 0
  } else {
    data$event[i] <- 1
  }
}

data$event[data$adoption==5] <- 0

model.1 <- glm(event ~ x1 + x2 + z1 + u + 1, family='binomial',data=data)
model.2 <- glm(event ~ x1 + x2 + contagion + 1, family='binomial',data=data)
model.3 <- glm(event ~ x1 + x2 + contagion + z1 +1, family='binomial',data=data)

summary(model.1)
summary(model.2)
summary(model.3)

stargazer(model.1,model.2,model.3,report='vc*')


##############################################
#### fitting ERGMs and latent space model #### 
##############################################

# adding in latent coordinates
# latent space models
# edges + mutual + ctriple +ttriple + triangle
net <- net.a.1
net <- net.a.2
net <- net.a.3
net <- net.a.4

net <- net.d.1
net <- net.d.2
net <- net.d.3
net <- net.d.4

net <- net.f.1
net <- net.f.2
net <- net.f.3
net <- net.f.4

table(degree(net.f.1,cmode='outdegree'))

# selected model 
e.a.1 <- ergm(net ~ edges + mutual + ostar(3) +
            nodematch("med_sch_yr") + nodefactor('med_sch_yr') +
            nodematch('proximity') + nodefactor('proximity') +
            nodematch('specialty') + nodefactor('specialty') +
            nodefactor('friends') + 
            nodefactor('community') + 
            nodefactor('jours') +
            nodefactor('free_time'))

e.s <- simulate(e.f.1)
summary(e.s)
plot(e.s)
gof(m.a.1)
plot(gof(e.a.1))

summary(net.a.4~triangle)
mixingmatrix(net.a,'specialty')

stargazer(e.f.1,e.f.2,e.f.3,e.f.4,single.row=T,report=('vc*'))
stargazer(e.a.1,single.row=T,report=c('vc*'))

mcmc.diagnostics(e.f.2)
stargazer(gof(e.a.1~idegree))

stargazer

# full model
e <- ergm(net ~ edges + mutual +
             nodematch("med_sch_yr") + nodefactor('med_sch_yr') + 
             nodematch('meetings') + nodefactor('meetings')  +
             nodematch('jours') + nodefactor('jours') +
             nodematch('free_time') + nodefactor('free_time') +
             nodematch('discuss') + nodefactor('discuss') +
             nodematch('clubs') + nodefactor('clubs') +
             nodematch('friends') + nodefactor('friends') + 
             nodematch('community') + nodefactor('community') +
             nodematch('patients') + nodefactor('patients') +
             nodematch('proximity') + nodefactor('proximity') + 
             nodematch('specialty') + nodefactor('specialty'))

e <- ergm(net ~ edges + mutual + 
            nodematch("med_sch_yr") + 
            nodematch('meetings') + 
            nodematch('jours') + 
            nodematch('free_time') + 
            nodematch('discuss') + 
            nodematch('clubs') + 
            nodematch('friends') +  
            nodematch('community') + 
            nodematch('patients') + 
            nodematch('proximity') +  
            nodematch('specialty'))
summary(e.a.1)

gof(e.a.1,GOF=~idegree + odegree + distance + triadcensus)

par(mfrow=c(1,1))
plot(gof(m.a.1,GOF=~idegree + odegree + distance + triadcensus),main=NULL)

?plot
e <- ergm(net ~ edges + 
            nodematch("med_sch_yr") + nodefactor('med_sch_yr') + 

            nodematch('jours') + nodefactor('jours') +
            nodematch('free_time') + nodefactor('free_time') +

            nodematch('friends') + nodefactor('friends') + 
            nodematch('community') + nodefactor('community') +

            nodematch('proximity') + nodefactor('proximity') + 
            nodematch('specialty') + nodefactor('specialty'))
summary(e.a.1)

# m <- ergmm(net.a.1~euclidean(d=2)+latentcov('med_sch_yr')+latentcov('meetings')+
#              latentcov('jours')+latentcov('free_time')+latentcov('discuss')+
#              latentcov('clubs')+latentcov('friends')+latentcov('community')+
#              latentcov('patients')+latentcov('proximity')+latentcov('specialty'),verbose=T)

m.d.1 <- ergmm(net ~ euclidean(d=3) + edges + 
                 nodematch("med_sch_yr",diff=T) + nodefactor('med_sch_yr') + 
                 nodematch('meetings',diff=T) + nodefactor('meetings')  +
                 nodematch('jours',diff=T) + nodefactor('jours') +
                 nodematch('free_time',diff=T) + nodefactor('free_time') +
                 nodematch('discuss',diff=T) + nodefactor('discuss') +
                 nodematch('clubs',diff=T) + nodefactor('clubs') +
                 nodematch('friends',diff=T) + nodefactor('friends') + 
                 nodematch('community',diff=T) + nodefactor('community') +
                 nodematch('patients',diff=T) + nodefactor('patients') +
                 nodematch('proximity',diff=T) + nodefactor('proximity') + 
                 nodematch('specialty',diff=T) + nodefactor('specialty'),
                 tofit=c('mle'),verbose=T)

######################################
#### computing latent coordinates ####
######################################
net <- net.d.1
m.d.1 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') + 
                 nodefactor('specialty'),
                 tofit=c('mle'),verbose=T)
net <- net.d.2
m.d.2 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') + 
                 nodefactor('specialty'),
                 tofit=c('mle'),verbose=T)
net <- net.d.3
m.d.3 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') + 
                 nodefactor('specialty'),
                 tofit=c('mle'),verbose=T)
net <- net.d.4
m.d.4 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') + 
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
print('net.d is done!')

net <- net.f.1
m.f.1 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
net <- net.f.2
m.f.2 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
net <- net.f.3
m.f.3 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
net <- net.f.4
m.f.4 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
print ('net.f is done!')

net <- net.a.1
m.a.1 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + ostar(1:3) +
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)

net <- net.a.2
m.a.2 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
net <- net.a.3
m.a.3 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
net <- net.a.4
m.a.4 <- ergmm(net ~ euclidean(d=2) + edges + rreceiver + 
                 mutual + ostar(3) +
                 nodefactor('med_sch_yr') + 
                 nodefactor('meetings')  +
                 nodefactor('jours') +
                 nodefactor('free_time') +
                 nodefactor('discuss') +
                 nodefactor('clubs') +
                 nodefactor('friends') + 
                 nodefactor('community') +
                 nodefactor('patients') +
                 nodefactor('proximity') +
                 nodefactor('specialty'),
               tofit=c('mle'),verbose=T)
summary(m.a.4)


sm <- simulate(e.a.4)
plot(m.f.4,what='mle',main='',print.formula=F)
m.a.3$mle

summary(e.a.1)
summary(m.a.1)

mixingmatrix(net.a,'med_sch_yr')
mixingmatrix(net.a,'proximity')
mixingmatrix(net.a,'specialty')

# merge coordinates with data set

data.a.1 <- data.a.1[,-c(29,30)]
data.a.2 <- data.a.2[,-c(29,30)]
data.a.3 <- data.a.3[,-c(29,30)]
data.a.4 <- data.a.4[,-c(29,30)]

data.d.1 <- data.d.1[,-c(29,30)]
data.d.2 <- data.d.2[,-c(29,30)]
data.d.3 <- data.d.3[,-c(29,30)]
data.d.4 <- data.d.4[,-c(29,30)]

data.f.1 <- data.f.1[,-c(29,30)]
data.f.2 <- data.f.2[,-c(29,30)]
data.f.3 <- data.f.3[,-c(29,30)]
data.f.4 <- data.f.4[,-c(29,30)]

id <- 1:117 #117
id <- 118:167 # 50
id <- 168:211 # 44
id <- 212:246 # 35

head(data.d.1)
head(data.a.1)
m.a.1$mle$Z

m <- m.f.4
data <- data.f.4

z <- as.data.frame(cbind(id,m$mle$Z))
colnames(z) <- c('id','z1','z2')
data <- merge(data,z,by='id')

data.f.4 <- data
head(data.a.1)

m.d.2 <- ergmm(net ~ euclidean(d=2),tofit=c('mle'))
plot(m.d.2)

m.d.1$mle$Z
m.d.2$mle$Z

summary(m.d.2)
summary(m.d.1)

# adding in latent coordinates

id <- 1:117 #117
id <- 118:167 # 50
id <- 168:211 # 44
id <- 212:246 # 35

coordinate <- m$mle$Z

coordinate <- cbind(id,coordinate)
coordinate <- as.data.frame(coordinate)

names(coordinate) <- c('id','d1','d2','d3')
dim(coordinate)
head(coordinate)

dim(test)
test[1:10,29]
test <- test[,-c(29:34)]
test <- merge(test,coordinate,by='id')
head(test)
dim(test)

apply(m$mle$Z,1,mean)

# treatmetn process 
model <- glm(event ~ factor(med_sch_yr) + log(jours) + factor(patients) +
               contagion + d1 + d2 + d3 + 1, family="binomial", data=test)
summary(model)
l <- polr(factor(contagion.dummy) ~ d1 + d2 + d3, data=test,Hess=T)
summary(l)

#  + factor(med_sch_yr) +
#    factor(jours) + factor(specialty)
#  + factor(meetings) +
#    factor(discuss) + factor(clubs) +
#   factor(friends) + factor(patients) +
#   factor(proximity)
