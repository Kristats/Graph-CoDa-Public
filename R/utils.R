# load libraries for calculations
dyn.load("fortran/graph_lasso.dll")
library(igraph)
library(plyr)
library(QR)

# for plotting
library(pheatmap)
library(wesanderson)
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(patchwork)


# function to calculate a whole path of weights 
weights <- function(x,
                    positiv = TRUE, 
                    symmetric = TRUE, 
                    Lsq       = TRUE,
                    normp     = "2",
                    nlambdas = 10,
                    lambda_rat = 0.01,
                    opt_params = list(maxoutit = 1e5, epsout = 1e-8, 
                                 lassoMaxinit = 1e5, lassoEpsin = 1e-8,
                                 admmMaxoutit = 1e5, admmEpsout = 1e-8)){
  
  # prepare
  x = exp(scale(log(x), scale = F))
  xlog   = log(x)
  cxlog  = scale(xlog, scale = F)
  p      = dim(x)[2]
  N      = dim(x)[1]
  
  # initialze as input
  Sigma = cov(xlog,use = "pairwise.complete.obs")
  #Sigma[is.na(Sigma)] = 0
  
  # change  normp to numeric
  if(normp == "inf"){ normp = -1 }else{ normp = as.numeric(normp) }
  
  # maximal lambda 
  if(Lsq == TRUE){
          if(symmetric == FALSE){
                  
            # calculate Z %*% ystart
            ystart = .Fortran("leftvec", as.vector(rep(0,p^2)), as.matrix(Sigma), as.double(rep(1,p*(p-1))), as.integer(p))[[1]]
            Zy     = sapply(1:(p*(p-1)), function(idx){
             
             locidx = Reduce("c",.Fortran("rowidx",as.integer(idx),as.integer(0),as.integer(0),as.integer(p)))[2:3]
             (ystart[locidx[1]] - ystart[locidx[2]]) / p^2
             
            })
            
            # take biggest absolute value
            lambda_max = max(abs(Zy))
            
          }else{
            
            # calculate Z %*% ystart
            ystart = .Fortran("leftvecSym", as.vector(rep(0,p^2)), as.matrix(Sigma), as.double(rep(1,p*(p-1)/2)), as.integer(p))[[1]]
            Zy     = sapply(1:(p*(p-1)/2), function(idx){
              
              gidx = Reduce("c",.Fortran("colidxSymToNonsym",as.integer(idx),as.integer(0),as.integer(0),as.integer(p)))[2:3]
              locidx1 = Reduce("c",.Fortran("rowidx",as.integer(gidx[1]),as.integer(0),as.integer(0),as.integer(p)))[2:3]
              locidx2 = Reduce("c",.Fortran("rowidx",as.integer(gidx[2]),as.integer(0),as.integer(0),as.integer(p)))[2:3]
              (ystart[locidx1[1]] - ystart[locidx1[2]] + ystart[locidx2[1]] - ystart[locidx2[2]]) / p^2
              
            })
            
            # take biggest absolute value
            lambda_max = 2 * max(abs(Zy))
          }
  }else{
          
          # calculate Z %*% ystart
            ystart = c(Sigma)
            Zy     = sapply(1:(p*(p-1)/2), function(idx){
              
              gidx = Reduce("c",.Fortran("colidxSymToNonsym",as.integer(idx),as.integer(0),as.integer(0),as.integer(p)))[2:3]
              locidx1 = Reduce("c",.Fortran("rowidx",as.integer(gidx[1]),as.integer(0),as.integer(0),as.integer(p)))[2:3]
              locidx2 = Reduce("c",.Fortran("rowidx",as.integer(gidx[2]),as.integer(0),as.integer(0),as.integer(p)))[2:3]
              (ystart[locidx1[1]] - ystart[locidx1[2]] + ystart[locidx2[1]] - ystart[locidx2[2]]) / p^2
              
            })
            
            # take biggest absolute value
            lambda_max = 2 * max(abs(Zy))
    
  }
  
  
  
  # compute lambdas
  lambdas = exp(seq(log(lambda_max * lambda_rat), log(lambda_max), length.out = nlambdas))
  lambdas = sort(lambdas, decreasing = TRUE)
  
  # set additional opt params
  nlambdas = length(lambdas)
  Ws       = array(0, dim = c(p, p, length(lambdas)))
  
  # solve problem
  if(Lsq == TRUE){
    
          if(symmetric == FALSE){
            
             
            sol = .Fortran("varGraph",
                        Sigma = as.matrix(Sigma), 
                        Ws    = as.array(Ws), 
                        lambdas = as.vector(lambdas), 
                        nlambdas = as.integer(nlambdas),
                        p = as.integer(p), 
                        positiv = as.integer(positiv),
                        as.double(normp),
                        as.integer(opt_params$maxoutit),
                        as.double(opt_params$epsout), 
                        as.double(0),
                        as.integer(opt_params$lassoMaxinit), 
                        as.double(opt_params$lassoEpsin), 
                        as.double(0),
                        as.integer(opt_params$admmMaxoutit), 
                        as.double(opt_params$admmEpsout), 
                        as.double(0),
                        as.integer(0))
            
            
          }
          if(symmetric == TRUE){
            
            sol = .Fortran("varGraphSym",
                        Sigma = as.matrix(Sigma), 
                        Ws    = as.array(Ws), 
                        lambdas = as.vector(lambdas), 
                        nlambdas = as.integer(nlambdas),
                        p = as.integer(p), 
                        positiv = as.integer(positiv),
                        as.double(normp),
                        as.integer(opt_params$maxoutit),
                        as.double(opt_params$epsout), 
                        as.double(0),
                        as.integer(opt_params$lassoMaxinit), 
                        as.double(opt_params$lassoEpsin), 
                        as.double(0),
                        as.integer(opt_params$admmMaxoutit), 
                        as.double(opt_params$admmEpsout), 
                        as.double(0),
                        as.integer(0))
            
          }
            
  }else{
         
            sol = .Fortran("singleVar",
                        Sigma = as.matrix(Sigma), 
                        Ws    = as.array(Ws), 
                        lambdas = as.vector(lambdas), 
                        nlambdas = as.integer(nlambdas),
                        p = as.integer(p), 
                        positiv = as.integer(positiv),
                        as.double(normp),
                        as.integer(opt_params$lassoMaxinit), 
                        as.double(opt_params$lassoEpsin), 
                        as.double(0),
                        as.integer(opt_params$admmMaxoutit), 
                        as.double(opt_params$admmEpsout), 
                        as.double(0),
                        as.integer(0))

  }
  
  ## --- calculate finally some information on resulting weights
  
  # rescale the Ws so that diagonal sums to p 
  sol$Ws = aperm(aaply(sol$Ws,3,function(W){ diag.sum = sum(abs(rowSums(W))); if(diag.sum != 0){ W / diag.sum * p }else{ W } }), perm = c(3,2,1))
    
    
  # calculate laplcians
  Ls = lapply(1:nlambdas, function(idx){ W = sol$Ws[,,idx]; diag(rowSums(W)) - W })
  
  
  # dofs or (sparsity)
  dfs = lapply(1:nlambdas, function(idx) { df = sum(sol$Ws[,,idx] != 0) })
  dfs = Reduce("c",dfs)
  
  if(symmetric == TRUE || Lsq){ dfs = dfs /2 }
  
  # calculate R2s
  R2s = lapply(1:nlambdas, function(idx){
    
          locmod = lm.fit(Ls[[idx]] %*% t(xlog), t(xlog))
          1 - sum((locmod$residuals)^2) / sum( cxlog^2 )
  })
  R2s = Reduce("c",R2s)

   
  #  variable(row)  dofs
  var.dfs = lapply(1:nlambdas, function(idx) { apply(sol$Ws[,,idx], 1, function(u){ sum( u != 0 ) }) })
  var.dfs = Reduce("cbind",var.dfs)
  rownames(var.dfs) = colnames(x)
  colnames(var.dfs) = paste0("lambda", 1:nlambdas)
  
  # degree weights
  degrees.weights = lapply(Ls, function(L){ diag(L) })
  degrees.weights = Reduce("cbind",degrees.weights)
  rownames(degrees.weights) = colnames(x)
  colnames(degrees.weights) = paste0("lambda", 1:nlambdas)
  
  # ICs
  ics = rbind(Reduce("c",lapply(1:nlambdas, function(idx){ -R2s[[idx]] + dfs[[idx]] * 1/ N  })),
              Reduce("c",lapply(1:nlambdas, function(idx){ -R2s[[idx]] + dfs[[idx]] * log(N)/ N  })))
  rownames(ics) = c("aic","bic")
  
  # graphs
  graphs = lapply(1:nlambdas, function(idx) { 
    
  
    # deleteing zero nodes
    nde = unname(which(var.dfs[,idx]!=0))
    
    if(length(nde) == 0){ 
      return(NA) 
      }else if(symmetric == TRUE){
      
        # calculate graph and edge width 
        graph = graph_from_adjacency_matrix((!!sol$Ws[nde,nde,idx]) * 1, 
                                            mode = "undirected", 
                                            diag = F)
        V(graph)$label = colnames(x)[nde]
        edge.list = t(apply(get.edgelist(graph),1,function(u){as.numeric(u)}))
        thickness = apply(edge.list, 1, function(inds){ sol$Ws[inds[1],inds[2],idx] })
        E(graph)$weight = (thickness - min(thickness)) / (max(thickness)-min(thickness) + 1e-5) * 5 + 1
        E(graph)$width = E(graph)$weight

        return(graph)
          
    }else{
      
         # calculate graph and edge width 
        graph = graph_from_adjacency_matrix((!!sol$Ws[nde,nde,idx]) * 1, 
                                            mode = "directed", 
                                            diag = F)
        V(graph)$label = colnames(x)[nde]
        edge.list = t(apply(get.edgelist(graph),1,function(u){as.numeric(u)}))
        thickness = apply(edge.list, 1, function(inds){ sol$Ws[inds[1],inds[2],idx] })
        E(graph)$weight =  (thickness - min(thickness)) / (max(thickness)-min(thickness) + 1e-5) * 5 + 1
        E(graph)$width = E(graph)$weight
        
        return(graph)
      
    }
    
    
    })

   
  
  
   
  return(list(Ws = sol$Ws, 
              Ls = Ls, 
              R2s = R2s, 
              ics = ics,
              dfs = dfs, 
              var.dfs = var.dfs, 
              degrees.weights = degrees.weights,
              graphs = graphs,
              positiv = positiv,
              symmetric = symmetric,
              Lsq = Lsq,
              normp = normp,
              nlambdas = nlambdas,
              lambda_rat = lambda_rat,
              opt_params = opt_params))
  
}



# function to plot the graph of the solution with given index
plot.graph <- function(model, idx, layout){
  
  
  graph = model$graphs[[idx]]
  if(model$symmetric == T){ pl =  plot(graph, layout = layout) }else{ pl = plot(graph, edge.curved= 0.2, layout = layout) }
  
  return(pl)
}



# function to calculate gilr (1) given a L matrix
# Lsq = TRUE -> laplacian is squared
gilr1 <- function(L, pivotvar = 1, Lsq = TRUE, tol = 1e-6){
  
   
   if(Lsq == FALSE & (dim(L)[1]!=dim(L)[2])){ stop("L must be symmetric if not squared!")  }
   if(Lsq == TRUE){  L.prod = t(L)%*%L }else{ L.prod = L }

  
   # compute permuation
   perm = 1:dim(L)[2]
   perm[pivotvar] = 1
   perm[1] = pivotvar
   
   # compute
   L.prod = L.prod[perm,perm]
   L.eig  = eigen(L.prod)
   L.root = L.eig$vectors %*% diag(sqrt(abs(L.eig$values))) %*% t(L.eig$vectors)
   qr.decomp  = QR(L.root)
   C.mat = qr.decomp$R
   C.mat[abs(C.mat) / max(abs(C.mat)) < tol] = 0
   
   # zero rows - so to get a surjective map 
   nzrws     = apply(C.mat, 1, function(u){val = max(abs(u)); if(val < tol){val = 0}; return(val) })
   nzrws.idx = unique(c(1,sort(order(nzrws, decreasing = TRUE)[1:min(qr(C.mat)$rank,length(which(nzrws!=0)))])))#which(nzrws!=0)
   C.mat     = C.mat[nzrws.idx,]
    
   # rename
   colnames(C.mat) = perm
   rownames(C.mat) = perm[nzrws.idx]
      
   return(C.mat)
}


# function to plot matrix of gilr (1) given L matrix 
plot.gilr1 <- function(L, pivotvar = 1, Lsq = TRUE, scaled = TRUE, tol = 1e-6){
 
  
   # compute 
   C.mat = gilr1(L, pivotvar = pivotvar, Lsq = Lsq, tol = tol)
   
   # plot heatmap
   #sc = sapply(diag(C.mat),function(u){ if(u != 0){abs(1/u)}else{0} }) 
   pal = wes_palette("Zissou1", 100, type = "continuous")
   if(scaled == TRUE){ sc = apply(C.mat,1,function(u){sc = max(abs(u)); if(sc != 0){ return(1/sc) }else{ return(0) } }) }else{ sc = 1 }
   pheatmap(sc * C.mat, cluster_rows = FALSE, cluster_cols = FALSE, color = pal)
   
}


# function to plot first gilr for each pivoted 
plot.gilr1.all <- function(L, Lsq = TRUE, scaled = TRUE, cluster.cols = FALSE, tol = 1e-6){
  
  
   if(Lsq == FALSE & (dim(L)[1]!=dim(L)[2])){ stop("L must be symmetric if not squared!")  }
   if(Lsq == TRUE){  L.prod = t(L)%*%L }else{ L.prod = L }

  
   # delete zero rows/cols
   nze = apply(L.prod, 1, function(u){val = max(abs(u)); if(val < tol) val = 0; return(val)})
   nze.idx = which(nze != 0)
  
   #create plot matrix
   plot.mat = t(sapply(nze.idx, function(pivotvar){ 
     
     perm = 1:dim(L)[2]; perm[pivotvar] = 1; perm[1] = pivotvar
     gilr1(L, pivotvar = pivotvar, Lsq = Lsq, tol = tol)[1, perm] 
   
   } ))
   plot.mat = plot.mat[,nze.idx]
   
   # name
   colnames(plot.mat) = colnames(L)[nze.idx] # as.character(1:dim(L)[2])
   rownames(plot.mat) = rownames(L)[nze.idx]
   
   # plot
   #sc = sapply(diag(C.mat),function(u){ if(u != 0){abs(1/u)}else{0} }) 
    pal = wes_palette("Zissou1", 100, type = "continuous")
   if(scaled == TRUE){ sc = apply(plot.mat,1,function(u){sc = max(abs(u)); if(sc != 0){ return(1/sc) }else{ return(0) } }) }else{ sc = 1 }
   pheatmap(sc * plot.mat, cluster_rows = FALSE, cluster_cols = cluster.cols, show_colnames = T, color = pal)
     
  
}


# function to calculate gilr (2) given a L matrix
gilr2 <- function(L,Lsq = TRUE){
  
   if(Lsq == FALSE & (dim(L)[1]!=dim(L)[2])){ stop("L must be symmetric if not squared!")  }
   if(Lsq == TRUE){  L.prod = t(L)%*%L }else{ L.prod = L }

   # compute decomposition
   L.eig  = eigen(L.prod)
   L.rk   = qr(L.prod )$rank
   C.mat = diag(sqrt(abs(L.eig$values[1:L.rk]))) %*% t(L.eig$vectors[,1:L.rk])
 
   # names
   colnames(C.mat) = colnames(L)
   
   # ret
   return(C.mat)
  
}

# function to calculate inverse if gilr (2) given a L matrix
gilr2.inv <- function(L, Lsq = TRUE){
  
  
   if(Lsq == FALSE & (dim(L)[1]!=dim(L)[2])){ stop("L must be symmetric if not squared!")  }
   if(Lsq == TRUE){  L.prod = t(L)%*%L }else{ L.prod = L }

   # compute decomposition
   L.eig  = eigen(L.prod)
   L.rk   = qr(L.prod )$rank
   L.scs  = 1 / sqrt(abs(L.eig$values[1:L.rk]))
   C.mat.inv = L.eig$vectors[,1:L.rk] %*% diag(L.scs) 
 
   # ret
   return(C.mat.inv )

}



# function to plot matrix of gilr (2) given L matrix 
plot.gilr2 <- function(L, Lsq = TRUE, scaled = FALSE, cluster.cols = FALSE, tol = 1e-6){
 
   if(Lsq == FALSE & (dim(L)[1]!=dim(L)[2])){ stop("L must be symmetric if not squared!")  }
   if(Lsq == TRUE){  L.prod = t(L)%*%L }else{ L.prod = L }
     
   # delete zero rows/cols
   nze = apply(L.prod, 1, function(u){val = max(abs(u)); if(val < tol) val = 0; return(val)})
   nze.idx = which(nze != 0)
  

   # compute decomposition
   L.eig  = eigen(L.prod)
   L.rk   = qr(L.prod )$rank
   C.mat = diag(sqrt(abs(L.eig$values[1:L.rk]))) %*% t(L.eig$vectors[,1:L.rk])
   C.mat = C.mat[,nze.idx]
   
   # names
   colnames(C.mat) = colnames(L)[nze.idx]
   
   # plot heatmap
   pal = wes_palette("Zissou1", 100, type = "continuous")
   if(scaled == TRUE){ sc = 1 / sqrt(abs(L.eig$values[1:L.rk])) }else{ sc = 1 }
   p.heatmap = pheatmap(sc * C.mat, cluster_rows = FALSE, cluster_cols = cluster.cols, show_colnames = T, color = pal)
   p.eig     = ggplot(data = data.frame(x = 1:L.rk, y = (abs(L.eig$values[1:L.rk]))) / sum(abs(L.eig$values[1:L.rk])) ,aes(x,y)) + 
               geom_point() + 
               ylab("Sorted Eigenvalues")+
               xlab("Number of Eigenbasis")
   grid.arrange(as.ggplot(p.heatmap), p.eig, nrow = 1)
   
}




























