# load libraries
library(stringr)
library(gridExtra)
library(robCompositions)
library(ggplot2)
library(GGally)
source("R/utils.R")


######################################################### Prepare nugent data set #################################################################
# nugent data set example
Xvar <- read.csv("data/otutable.txt", header=T,row.names=1, sep="\t", comment.char="")
Xvar <- t(Xvar)
yvar <- read.csv("data/task-nugent-score.txt",sep="\t")
Xvar <- Xvar[match(yvar$X.SampleID,rownames(Xvar)),]
y <- yvar$Var
X <- Xvar[,colSums(Xvar != 0) >= nrow(Xvar)*0.05] # select OTU taxa present in at least in the 5% of the samples

# exclude all variable with less than 2 values > DL
w <- apply(X, 2, function(x) sum(x != 0) > 1 )
X <- X[, w]
x = X
x = as.data.frame(x)
x[x==0] = 0.65 * min(x[x != 0])

# center
p = dim(x)[2]
N = dim(x)[1]
x = exp(scale(log(x), scale = F))

# set names appropriately
namesstr = lapply(colnames(x),function(str){str_split(str, "_strain")[[1]][[1]]})
namesstr = do.call("c",namesstr)
namesstr = lapply(namesstr,function(str){str_split(str, ".1_")[[1]][2]})
namesstr = do.call("c",namesstr)
namesstr[40]  = str_split(str_split(colnames(x)[40],".2_")[[1]][2],"_strain")[[1]][1]
namesstr[10]  = str_split(namesstr[10],"_3")[[1]][1]
namesstr      = sapply(namesstr, function(u){ loc.names = str_split(u,"_" )[[1]]; paste(loc.names[1],"\n",loc.names[2])  })
colnames(x) = namesstr




######################################################### Fit models and show basic results #########################################################

# compute path of solutions 
nlambdas = 40
mod1  = weights(x, positiv = TRUE, symmetric = TRUE, 
                    Lsq       = TRUE, normp     = "2", nlambdas = nlambdas,
                    lambda_rat = 0.01, opt_params = list(maxoutit = 1e5, epsout = 1e-8, 
                                 lassoMaxinit = 1e6, lassoEpsin = 1e-8,
                                 admmMaxoutit = 1e6, admmEpsout = 1e-8))


mod2  = weights(x, positiv = TRUE, symmetric = F, 
                    Lsq       = TRUE, normp     = "2", nlambdas = nlambdas,
                    lambda_rat = 0.01, opt_params = list(maxoutit = 1e5, epsout = 1e-8, 
                                 lassoMaxinit = 1e6, lassoEpsin = 1e-8,
                                 admmMaxoutit = 1e6, admmEpsout = 1e-8))


# get best model indices
ind1aic = which.min(mod1$ics[1,])
ind1bic = which.min(mod1$ics[2,])
ind2aic = which.min(mod2$ics[1,])
ind2bic = which.min(mod2$ics[2,])


# plot both df vs r2s for bic and with horizontals
r2sdf.data = data.frame(x = c(mod1$dfs, mod2$dfs), y = c( mod1$R2s, mod2$R2s), z  = c(rep(1,nlambdas),rep(2,nlambdas)))
ggplot(r2sdf.data, aes(x = x,y = y, col = as.factor(z))) + 
  scale_color_manual(values = c("1" = "black","2" = "blue"))+ 
  theme(legend.position="none",axis.text = element_text(size=13,face="bold"),axis.title = element_text(size=18,face="bold"))+
  geom_vline(xintercept = mod1$dfs[ind1bic],col = "black") + 
  geom_vline(xintercept = mod2$dfs[ind2bic],col = "blue") + 
  geom_point() + 
  xlab(expression("df("*lambda*")")) + 
  ylab(expression(R^2*"("*lambda*")"))+
  ylim(c(0,1))


# R2s for chosen models
mod1$R2s[ind1bic]
mod2$R2s[ind2bic]

# plot ics statistics
ggplot(data.frame(x = 1:nlambdas, y = mod1$ics[1,]), aes(x,y)) + geom_point()
ggplot(data.frame(x = 1:nlambdas, y = mod1$ics[2,]), aes(x,y)) + geom_point()
ggplot(data.frame(x = 1:nlambdas, y = mod2$ics[1,]), aes(x,y)) + geom_point()
ggplot(data.frame(x = 1:nlambdas, y = mod2$ics[2,]), aes(x,y)) + geom_point()


# plot graphs
plot.graph(mod1,ind1aic, layout = layout.circle)
plot.graph(mod1,ind1bic, layout = layout.circle)
plot.graph(mod2,ind2aic, layout = layout.circle)
plot.graph(mod2,ind2bic, layout = layout.circle)



## fit linear model 
xtrans = pivotCoord(x)
y      = scale(y)
N      = dim(x)[1]
K      = 100

# init errors
errsc = rep(0,K)
errsg1 = rep(0,K)
errsg2 = rep(0,K)


# run linear models
for(i in 1:K){
  
  # sample
  print(i)
  indtest  =  sample.int(N,size = ceiling(N*0.2)) 
  indtrain =  setdiff(1:N,indtest)  
  
   # regular CoDa model
   x.train   = xtrans[indtrain,]
   x.test    = xtrans[indtest,]
   y.train   = y[indtrain]
   y.test    = y[indtest]
  
   mod.lm   = lm(z ~ .,data = data.frame(x.train,z = y.train))
   coeffs   = as.vector(mod.lm$coefficients)
   fit      = as.matrix(x.test) %*% coeffs[-1] + coeffs[1]
   errsc[i] = mean((y.test-fit)^2)

   
   # first model
   ind.model = ind1bic
   mod       = mod1
   gilr2.mat = gilr2(mod$Ls[[ind.model]], Lsq = TRUE)
   x.train   = t(gilr2.mat %*% t(log(x[indtrain,])))
   x.test    = t(gilr2.mat %*%t(log(x[indtest,])))
   y.train   = y[indtrain]
   y.test    = y[indtest]

   mod.lm    = lm(z ~ .,data = data.frame(x.train, z = y.train))
   coeffs    = as.vector(mod.lm$coefficients)
   fit       = x.test %*% coeffs[-1] + coeffs[1]
   errsg1[i] = mean((y.test-fit)^2)
  
   
   # second model
   ind.model = ind2bic
   mod       = mod2
   gilr2.mat = gilr2(mod$Ls[[ind.model]], Lsq = TRUE)
   x.train   = t(gilr2.mat %*% t(log(x[indtrain,])))
   x.test    = t(gilr2.mat %*%t(log(x[indtest,])))
   y.train   = y[indtrain]
   y.test    = y[indtest]

   mod.lm    = lm(z ~ .,data = data.frame(x.train, z = y.train))
   coeffs    = as.vector(mod.lm$coefficients)
   fit       = x.test %*% coeffs[-1] + coeffs[1]
   errsg2[i] = mean((y.test-fit)^2)
   
   
} 


# show solutions of models
boxplot.data   = data.frame(mse = c(errsc,errsg1,errsg2),
                                model = c(rep("CoDa",K),
                                          rep("positive/symmetric",K),
                                          rep("positive/not symmetric",K)))

ggplot(data = boxplot.data, aes(x = reorder(model,-mse), y = mse)) + 
  geom_boxplot() + 
  ggtitle("") + 
  ylab("Mean squared error") + 
  xlab("") + 
  theme(axis.text = element_text(size=13,face="bold"),axis.title = element_text(size=18,face="bold"))



# plot results for best performing model 
ind.model = ind2bic
mod       = mod2

### ---- plot the graph 
graph = mod$graphs[[ind.model]]
V(graph)$label.cex = 0.8
plot(graph, vertex.size = 24, edge.curved= 0.1, layout = layout.circle)

### ---- highest degress
vertices.order = order(mod$degrees.weights[,ind.model],decreasing = T)[1:sum(mod$degrees.weights[,ind.model]!=0)]
dfs = mod$degrees.weights[vertices.order,ind.model]
dfs.names =  namesstr[vertices.order]

ggplot(data.frame(x = dfs.names,y = dfs),aes(x = reorder(x, y),y)) + 
  geom_point() + 
  xlab("") +
  ylab("Sum of weights")+
  coord_flip() + 
  theme(axis.text = element_text(size=13,face="bold"),axis.title = element_text(size=18,face="bold"))


### ---- show the coefficients of the model 
gilr2.mat     = gilr2(mod$Ls[[ind.model]], Lsq = TRUE)
gilr2.mat.inv = gilr2.inv(mod$Ls[[ind.model]], Lsq = TRUE)
nze           = apply(mod$Ls[[ind.model]], 1, function(u){val = max(abs(u)); if(val < 1e-6) val = 0; return(val)})
nze.idx       = which(nze != 0)
x.train       = t(gilr2.mat %*% t(log(x)))
y.train       = y

mod.lm = lm(z ~ .,data = data.frame(x.train, z = y.train))
coeffs = mod.lm$coefficients
coeffs.tr =  gilr2.mat.inv %*% coeffs[-1]
coeffs.gilr1 =  sapply(seq.int(nze.idx), function(idx){ 
                   pivotvar = nze.idx[idx]
                   gilr1.trans = gilr1(mod$Ls[[ind.model]], pivotvar = pivotvar, Lsq = TRUE)
                   vals        = gilr1.trans %*% coeffs.tr
                   sc          = sd(t(gilr1.trans[1,] %*% t(log(x)))) # scale of this coordinate
                   vals[1] * sc
                   })

# plot coefficients
names(coeffs.gilr1) = namesstr[nze.idx]
coeffs.df = data.frame(x = names(coeffs.gilr1), y =  (coeffs.gilr1))
coeffs.graph.pl = ggplot(data = coeffs.df, aes(x = reorder(x, y), y = y)) + 
                          geom_col() +
                          coord_flip() + 
                          xlab("") + 
                          ylab("")+
                          labs(title = "Graphical linear model")+
                          theme(plot.title = element_text(size=15), 
                                axis.text = element_text(size=11,face="bold"),
                                axis.title = element_text(size=11,face="bold"))


plot(coeffs.graph.pl)