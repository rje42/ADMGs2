## function that searches through ADMGs with simple edge manipulations
.searchADMG <-
  function(graph, dat, stat, criterion = "AIC", quietly = FALSE, r = TRUE, tol = 1e-6) {
    n = length(graph$v)
    current.graph = graph
    current.stat = stat
    moved = FALSE
    
    for (i in seq_len(n-1)) for (j in seq_len(n-i)+i) {
      if (!(j %in% anc(graph, i)) && !(i %in% pa(graph, j)) && !(i %in% sib(graph, j))) {
        if (!quietly) cat(i,"-->",j,criterion,"= ")
        # try.graph = graph
        try.graph = addEdges(graph, list(directed=list(c(i,j))))
        # try.graph$edges$directed = c(graph$edges$directed, list(c(i,j)))
        f1 = fitADMG(dat, try.graph, tol=tol, sparse=TRUE, quietly=TRUE, r = r)
        try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
        
        if (!quietly) cat(try.stat)
        
        if(try.stat < current.stat - 1e-3) {
          current.graph=try.graph; current.stat=try.stat
          if (!quietly) cat(" *")
          moved = TRUE
        }
        if (!quietly) cat("\n")
      }
      else if (i %in% pa(graph, j)) {
        if (!quietly) cat(i,"-|>",j,criterion,"= ")
        # try.graph = graph
        # try.graph$edges$directed = setdiff(graph$edges$directed, list(c(i,j)))
        try.graph = removeEdges(graph, list(directed=list(c(i,j))))
        f1 = fitADMG(dat, try.graph, tol=tol, sparse=TRUE, quietly=TRUE, r = r)
        try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
        
        if (!quietly) cat(try.stat)
        
        if(try.stat < current.stat - 1e-3) {
          current.graph=try.graph; current.stat=try.stat
          if (!quietly) cat(" *")
          moved = TRUE
        }
        if (!quietly) cat("\n")
      }
      
      if (!(i %in% anc(graph, j)) && !(j %in% pa(graph, i)) && !(j %in% sib(graph, i))) {
        if (!quietly) cat(i,"<--",j,criterion,"= ")
        # try.graph = graph
        # try.graph$edges$directed = c(graph$edges$directed, list(c(j,i)))
      try.graph = addEdges(graph, list(directed=list(c(j,i))))
        f1 = fitADMG(dat, try.graph, tol=tol, sparse=TRUE, quietly=TRUE, r = r)
        try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
        
        if (!quietly) cat(try.stat)
        
        if(try.stat < current.stat - 1e-3) {
          current.graph=try.graph; current.stat=try.stat
          if (!quietly) cat(" *")
          moved = TRUE
        }
        if (!quietly) cat("\n")
      }
      else if (j %in% pa(graph, i)) {
        if (!quietly) cat(i,"<|-",j,criterion,"= ")
        # try.graph = graph
        # try.graph$edges$directed = setdiff(graph$edges$directed, list(c(j,i)))
        try.graph = removeEdges(graph, list(directed=list(c(j,i))))
        f1 = fitADMG(dat, try.graph, tol=tol, sparse=TRUE, quietly=TRUE, r = r)
        try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
        
        if (!quietly) cat(try.stat)
        
        if(try.stat < current.stat - 1e-3) {
          current.graph=try.graph; current.stat=try.stat
          if (!quietly) cat(" *")
          moved = TRUE
        }
        if (!quietly) cat("\n")
      }
      
      
      if (!(j %in% sib(graph, i)) && !(i %in% pa(graph, j)) && !(j %in% pa(graph, i))) {
        if (!quietly) cat(i,"<->",j,criterion,"= ")
        # try.graph = graph
        # try.graph$edges$bidirected = c(graph$edges$bidirected, list(c(i,j)))
        try.graph = addEdges(graph, list(bidirected=list(c(i,j))))
        
        f1 = fitADMG(dat, try.graph, tol=tol, sparse=TRUE, quietly=TRUE, r = r)
        try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
        
        if (!quietly) cat(try.stat)
        
        if(try.stat < current.stat - 1e-3) {
          current.graph=try.graph; current.stat=try.stat
          if (!quietly) cat(" *")
          moved = TRUE
        }
        if (!quietly) cat("\n")
      }
      else if (j %in% sib(graph, i)) {
        if (!quietly) cat(i,"<|>",j,criterion,"= ")
        # try.graph = graph
        # try.graph$edges$bidirected = setdiff(graph$edges$bidirected, list(c(i,j), c(j,i)))
        try.graph = removeEdges(graph, list(directed=list(c(i,j))))
        
        f1 = fitADMG(dat, try.graph, tol=tol, sparse=TRUE, quietly=TRUE, r = r)
        try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
        
        if (!quietly) cat(try.stat)
        
        if(try.stat < current.stat - 1e-3) {
          current.graph=try.graph; current.stat=try.stat
          if (!quietly) cat(" *")
          moved = TRUE
        }
        if (!quietly) cat("\n")
      }
    }
    
    list(graph = current.graph, stat = current.stat, moved = moved)
  }

.searchADMG2 <-
  function(graph, dat, stat, criterion = "AIC", quietly = FALSE, r = TRUE, tol = 1e-2){
    
    graph <- withEdgeList(graph)
    
    n = length(graph$v)
    moved = FALSE
    
    best.stat <- stat
    best.graph <- graph
    
    
    if(!quietly) print("trying to add edges")
    
    for (i in seq_len(n-1)) for (j in seq_len(n-i)+i) {
      
      # try adding a directed edge from i to j
      if (!(j %in% anc(graph, i)) && !(i %in% pa(graph, j)) && !(i %in% sib(graph, j))) {
        
        # try.graph = graph
        # try.graph$edges$directed = c(graph$edges$directed, list(c(i,j)))
        try.graph = addEdges(graph, list(directed=list(c(i,j))))
        
        if(is_arid(try.graph)){
          
          #print("try.graph")
          #print(try.graph)
          
          f1 = fitADMG(dat, try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
          try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
          
          if(try.stat < best.stat){
            best.stat <- try.stat
            best.graph <- try.graph
            moved <- TRUE
          }
        }
      }
      
      # try adding a bidirected edge from i to j
      if (!(i %in% pa(graph, j)) && !(i %in% sib(graph, j))) {
        
        # try.graph = graph
        # try.graph$edges$bidirected = c(graph$edges$bidirected, list(c(i,j)))
        try.graph = addEdges(graph, list(bidirected=list(c(i,j))))
        
        
        if(is_arid(try.graph)){
          
          #print("try.graph")
          #print(try.graph)
          
          f1 = fitADMG(dat, try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
          try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
          
          if(try.stat < best.stat){
            best.stat <- try.stat
            best.graph <- try.graph
            moved <- TRUE
          }
        }
      }
    }
    
    if(!quietly) print("trying to remove/replace <-> edges")
    
    for(edge in graph$edges$bidirected){
      
      i <- edge[1]
      j <- edge[2]
      
      try.graph = graph
      
      # try removing i <-> j
      try.graph$edges$bidirected <- Filter(function(e) ! setequal(e, edge), try.graph$edges$bidirected)
      
      #print("try.graph")
      #print(try.graph)
      
      f1 = fitADMG(dat, try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
      try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
      
      if(try.stat < best.stat){
        best.stat <- try.stat
        best.graph <- try.graph
        moved <- TRUE
      }
      
      # try replacing i <-> j with i -> j
      if (!(j %in% anc(try.graph, i)) && !(i %in% pa(try.graph, j)) && !(i %in% sib(try.graph, j))) {
        
        add.try.graph = try.graph
        
        add.try.graph$edges$directed = c(try.graph$edges$directed, list(c(i,j)))
        
        if(is_arid(add.try.graph)){
          
          #print("add.try.graph")
          #print(add.try.graph)
          
          f1 = fitADMG(dat, add.try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
          try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
          
          if(try.stat < best.stat){
            best.stat <- try.stat
            best.graph <- try.graph
            moved <- TRUE
          }
        }
      }
      
      # try replacing i <-> j with j -> i
      if (!(i %in% anc(try.graph, j)) && !(j %in% pa(try.graph, i)) && !(j %in% sib(try.graph, i))) {
        
        add.try.graph = try.graph
        
        add.try.graph$edges$directed = c(try.graph$edges$directed, list(c(j,i)))
        
        if(is_arid(add.try.graph)){
          
          #print("add.try.graph")
          #print(add.try.graph)
          
          f1 = fitADMG(dat, add.try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
          try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
          
          if(try.stat < best.stat){
            best.stat <- try.stat
            best.graph <- add.try.graph
            moved <- TRUE
          }
        }
      }
    }
    
    if(!quietly) print("trying to remove/replace --> edges")
    
    for(edge in graph$edges$directed){
      
      i <- edge[1]
      j <- edge[2]
      
      try.graph = graph
      
      # try removing i -> j
      try.graph$edges$directed <- Filter(function(e) ! setequal(e, edge), try.graph$edges$directed)
      
      #print("try.graph")
      #print(try.graph)
      
      f1 = fitADMG(dat, try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
      try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
      
      if(try.stat < best.stat){
        best.stat <- try.stat
        best.graph <- try.graph
        moved <- TRUE
      }
      
      # try replacing i -> j with i <-> j
      add.try.graph = try.graph
      
      add.try.graph$edges$bidirected = c(try.graph$edges$bidirected, list(c(i,j)))
      
      if(is_arid(add.try.graph)){
        
        #print("add.try.graph")
        #print(add.try.graph)
        
        f1 = fitADMG(dat, add.try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
        try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
        
        if(try.stat < best.stat){
          best.stat <- try.stat
          best.graph <- try.graph
          moved <- TRUE
        }
      }
      
      # try replacing i -> j with j <- i
      if (!(i %in% anc(try.graph, j)) && !(j %in% pa(try.graph, i)) && !(j %in% sib(try.graph, i))) {
        
        add.try.graph = try.graph
        
        add.try.graph$edges$directed = c(try.graph$edges$directed, list(c(j,i)))
        
        if(is_arid(add.try.graph)){
          
          #print("add.try.graph")
          #print(add.try.graph)
          
          f1 = fitADMG(dat, add.try.graph, tol=tol, sparse=FALSE, quietly=TRUE, r = r)
          try.stat = switch(criterion, AIC=summary(f1, fisher=FALSE)$AIC, BIC=summary(f1, fisher=FALSE)$BIC)
          
          if(try.stat < best.stat){
            best.stat <- try.stat
            best.graph <- add.try.graph
            moved <- TRUE
          }
        }
      }
    }
    
    list(graph = best.graph, stat = best.stat, moved = moved)
  }
