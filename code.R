# Carga de datos correspondientes al arból obtenido por ML

phy<-read.newick(file="arbol.nwk")
phy2<-compute.brlen(phy, power = 0.7)
phy3<-compute.brlen(phy, power = 1)
plot(phy2)



# Visualización inicial

plotTree(phy2, fsize = 0.6, ftype = "i", lwd=0.7,colors=cols)

dotTree(phy, hongosclass, fsize=0.75, pt.cex=2, colors = color1, legend = F)
legend("bottomleft",c("Bryophytic", "Carton", "Domatium", "DSE", "Epilithic", "Epiphytic",
                       "Opportunistic", "Other"), pch=21,pt.bg=color1,pt.cex=1.5, cex =0.75)
nodelabels(phy$node.label,
           adj=c(1,-0.2),frame="none", cex =0.5)

## Branch Lengths Computation

plotTree(compute.brlen(phy), fsize = 0.6, ftype = "i")

dotTree(compute.brlen(phy, power = 1), hongosclass, fsize=0.75, pt.cex=0.75, colors = color1, legend = F)
legend("bottomleft",c("Bryophytic", "Carton", "Domatium", "DSE", "Epilithic", "Epiphytic",
                       "Opportunistic", "Other"), pch=21,pt.bg=color1,pt.cex=1.5, cex =0.75)


#######################################################################

# Corregir la longitud de las ramas para que sean detectadas por R

phy$edge.length[phy$edge.length==0]<-max(nodeHeights(phy))*1e-6

#######################################################################

# Acondicionar los datos

cepas$Ecology <- as.factor(cepas$Ecology)


rownames(cepas) <- cepas$ITS

## Extract character of interest

hongosclass<-setNames(cepas$Ecology,rownames(cepas))

## plot traitgram

phenogram(phy,hongosclass, colors = "green", ylim = c(8,1), spread.range = c(5,1), ftype="i",
          fsize=0.6, xlab="Relative Times", ylab="Ecology",  par(mar=c(4.2,4.1,1.1,0)))
nodelabels(phy$node.label,
           adj=c(1,-0.2),frame="none", cex =0.5)



##
fit.BM<-fastAnc(phy,hongosclass,CI=TRUE)
print(fit.BM)

## estimate ancestral state under BM model
fit.BM<-anc.ML(phy,hongosclass)
print(fit.BM)

#######################################################################

## LTT plots
dark.tree<-multi2di(phy2)



tiff("/home/federico/Escritorio/FiguraS2.tiff",
     width=10,height=9,units="cm",
     bg="white",res=320)
par(mfrow=c(1,1),mar=c(4.2,4.1,1.1,0))
obj<-ltt(phy3,gamma=T,show.tree=TRUE,log.lineages=FALSE,log="y",ylim=c(2,Ntip(phy2)),
         bty="l", col = "black", cex =10)
dev.off()




darter.tree<-multi2di(phy)
darter.tree

obj<-ltt(dark.tree,plot=FALSE)
obj

plot(obj,log.lineages=FALSE,log="y",main="LTT plot for darters",
     ylim=c(2,Ntip(darter.tree)))

h<-max(nodeHeights(darter.tree))
x<-seq(0,h,by=h/100)
b<-(log(Ntip(darter.tree))-log(2))/h
lines(x,2*exp(b*x),col="red",lty="dashed",lwd=2)


trees<-pbtree(b=b,n=Ntip(darter.tree),t=h,nsim=100,method="direct",
              quiet=TRUE)
obj<-ltt(trees,plot=FALSE)
plot(obj,col="grey",main="LTT of darters compared to simulated LTTs")
lines(c(0,h),log(c(2,Ntip(darter.tree))),lty="dashed",lwd=2,col="red")
## now let's overlay our original tree
ltt(darter.tree,add=TRUE,lwd=2)

-gammatest(obj)


## GTT plot

cha.gtt<-gtt(phy)
plot(cha.gtt)

######################################################################

# Evaluación de modelos ER, SYM y ARD por AIC
fitER<-fitMk(phy,hongosclass,model="ER")
fitARD<- fitMk(phy,hongosclass,model="ARD")
fitSYM<- fitMk(phy,hongosclass,model="SYM")


AIC(fitER, fitARD, fitSYM, fitSYM2)
summary(fitARD)

tiff("/home/federico/Escritorio/Figura2.tiff",
     width=20,height=16,units="cm",
     bg="white",res=320)

plot.fitMk(fitARD, lwd=0.4, par(mar=c(0,0,0,0),offset=1 ))
dev.off()




## do stochastic mapping
smap.trees <-make.simmap(phy2,hongosclass,model="SYM",
                                 nsim=1000, Q = "empirical")


## print a summary of the stochastic mapping
summary(smap.trees)
trees <- smap.trees


plot(trees2)
obj<-describe.simmap(trees,plot=T)
obj

# plot fan
cols <- setNames(c("red", "blue", "yellow", "pink", "purple", "green", "brown", "orange"),levels(cepas$Ecology))

plotSimmap(trees[[1]],type = "fan", fsize = 0.8, colors = cols,lwd=3)
nodelabels(pie=obj$ace,piecol=cols,cex=0.2)
add.simmap.legend(colors=cols,x=0.1*par()$usr[1],
                  y=0.1*par()$usr[4],prompt=FALSE)


## Plots stochastic mapping
color1 <- setNames(c("red", "blue", "yellow", "pink", "purple", "green", "brown", "orange"),levels(cepas$Ecology))

plot(smap.trees,color1,fsize=0.5,ftype="i",outline=F,
     lwd=3)

par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default

## plot a posterior probabilities of ancestral states
par(mar=c(5.1,4.1,4.1,0))
par(mar=c(4.2,4.1,1.1,0))

tiff("/home/federico/Escritorio/Figura1.tiff",
     width=11,height=22,units="cm",
     bg="white",res=320)
plot(summary(smap.trees),colors=color1,ftype="i", fsize =0.4, cex=0.4,no.margin=T,
     type="phylogram", lwd=0.4, offset=0.01)
legend(0.1,18,c("Bryophytic", "Carton", "Domatium", "DSE", "Epilithic", "Epiphytic", "Opportunistic", "Other"),
       pch=21,pt.bg=color1,pt.cex=1.2, cex =0.6,)
#par(family="CM Roman")
dev.off()


head(fonts())
help("extrafont")
#######################################

## reset par to default values
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

## create "contMap" object
ChaCont<-contMap(phy, hongosclass,plot=FALSE,res=200)

## change color scheme
ChaCont<-setMap(ChaCont, c("red", "blue", "yellow", "pink", "purple", "green", "brown", "orange"))

## Plot
plot(ChaCont, legend = FALSE, lwd =2.5, fsize =0.7)
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

########################################################

