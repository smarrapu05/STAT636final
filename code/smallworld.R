library(igraph)
whale<-read_graph("636/6BMG_edge.txt",format="edgelist")
c.w<-transitivity(whale)
l.w<-mean_distance(whale)

l.rand.w.store<-rep(0,1000)
c.rand.w.store<-rep(0,1000)
for (i in 1:1000){
  GG<-sample_gnm(n=154,m=length(E(whale)))
  l.rand.w.store[i]<-mean_distance(GG)
  c.rand.w.store[i]<-transitivity(GG)
}
c.rand.w<-mean(c.rand.w.store)
l.rand.w<-mean(l.rand.w.store)
S.w <- c.w * l.rand.w / l.w / c.rand.w

calm<-read_graph("636/1CLM_edge.txt",format="edgelist")
calm<-delete.vertices(calm,1:4)

c.c<-transitivity(calm)
l.c<-mean_distance(calm)

l.rand.c.store<-rep(0,1000)
c.rand.c.store<-rep(0,1000)
for (i in 1:1000){
  GG<-sample_gnm(n=length(V(calm)),m=length(E(calm)))
  l.rand.c.store[i]<-mean_distance(GG)
  c.rand.c.store[i]<-transitivity(GG)
}
c.rand.c<-mean(c.rand.c.store)
l.rand.c<-mean(l.rand.c.store)
S.c <- c.c * l.rand.c / l.c / c.rand.c

pgk<-read_graph("636/PGK_edge.txt",format="edgelist")
pgk<-delete.vertices(pgk,1)

c.p<-transitivity(pgk)
l.p<-mean_distance(pgk)

l.rand.p.store<-rep(0,1000)
c.rand.p.store<-rep(0,1000)
for (i in 1:1000){
  GG<-sample_gnm(n=length(V(pgk)),m=length(E(pgk)))
  l.rand.p.store[i]<-mean_distance(GG)
  c.rand.p.store[i]<-transitivity(GG)
}
c.rand.p<-mean(c.rand.p.store)
l.rand.p<-mean(l.rand.p.store)
S.p <- c.p * l.rand.p / l.p / c.rand.p

c.w
l.w
c.rand.w
l.rand.w
S.w

c.c
l.c
c.rand.c
l.rand.c
S.c

c.p
l.p
c.rand.p
l.rand.p
S.p