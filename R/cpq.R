
# backend for conditional probability queries.
conditional.probability.query = function(fitted, event, evidence, method,
                                         extra, probability = TRUE, cluster = NULL, debug = FALSE) {
  
  if (method == "ls") {
    
    if (!is.null(cluster)) {
      
      # get the number of slaves.
      s = nSlaves(cluster)
      # divide the number of particles among the slaves.
      batch = n = ceiling(extra$n / s)
      
      if (probability) {
        
        results = parSapply(cluster, seq(s),
                            function(x) {
                              
                              logic.sampling(fitted = fitted, event = event,
                                             evidence = evidence, n = n, batch = batch, debug = debug)
                              
                            })
        
        return(mean(results))
        
      }#THEN
      else {
        
        results = parLapply(cluster, seq(s),
                            function(x) {
                              
                              logic.distribution(fitted = fitted, nodes = event,
                                                 evidence = evidence, n = n, batch = batch, debug = debug)
                              
                            })
        
        return(do.call(rbind, results))
        
      }#ELSE
      
    }#THEN
    else {
      
      if (probability) {
        
        logic.sampling(fitted = fitted, event = event, evidence = evidence,
                       n = extra$n, batch = extra$batch, debug = debug)
        
      }#THEN
      else {
        
        logic.distribution(fitted = fitted, nodes = event, evidence = evidence,
                           n = extra$n, batch = extra$batch, debug = debug)
        
      }#ELSE
      
    }#ELSE
    
  }#THEN
  else if (method == "exact") {
    
    if (!is.fitted.discrete(fitted))
      stop("exact inference does only support discrete bns.")
    
    if (probability) {
      
      p = exact.prob(fitted = fitted, event = event, evidence = evidence, debug = debug)
      
      return(p)
      
    }#THEN
    else {
      
      probs = exact.cpt(fitted, event, evidence, debug)
      for(col in 1:(ncol(probs)-1))
        probs[, col] = factor(probs[, col])
      p.col = ncol(probs)
      cpt = table(probs[, -p.col], dnn=names(probs)[-p.col])
      cpt[1:length(cpt)] = probs[, p.col]
      cpt = prop.table(margin.table(cpt, 1:length(event)))
      
      return(cpt)
      
    }
    
  }#THEN
  
}#CONDITIONAL.PROBABILITY.QUERY

exact.cpt = function(fitted, event, evidence, debug = FALSE){
  
  extract.names = function(call) {
    names = character(0)
    for (x in as.list(call)[-1]) {
      if (is.call(x)) {
        names = c(names, extract.names(x))
      }#THEN
      else if (is.name(x)) {
        names = c(names, as.character(x))
      }#ELSE
    }#FOR
    return(names)
  }#EXTRACT.NAMES
  
  if(is.call(event) | is.logical(event)) {
    nodes = extract.names(event)
  }#THEN
  else {
    nodes = event
  }
  
  # Joint distribution of the target and conditional nodes knowing their
  # parents, their parent's parents, etc.
  nodes = intersect(unique(c(nodes, extract.names(evidence))), names(bn))
  to.check = nodes
  while(length(to.check) > 0) {
    for (node in to.check) {
      parents = bn[[node]]$parents
      to.check = setdiff(to.check, node)
      to.check = union(to.check, setdiff(parents, nodes))
      nodes = union(nodes, parents)
    }
  }
  nbnodes = length(nodes)
  
  cpt.table = expand.grid(lapply(nodes, function(x) {
    dimnames(bn[[x]]$prob)[[1]]
  }))
  names(cpt.table) = nodes
  
  m = nrow(cpt.table)
  
  # evaluate the expression defining the evidence.
  if (identical(evidence, TRUE))
    r = rep(TRUE, m)
  else
    r = eval(evidence, cpt.table, parent.frame())
  
  # double check that this is a logical vector.
  if (!is.logical(r))
    stop("evidence must evaluate to a logical vector.")
  # double check that it is of the right length.
  if (length(r) != m)
    stop("logical vector for evidence is of length ", length(r), " instead of ", m, ".")
  
  cpt.table = cpt.table[r, , drop=FALSE]
  
  for (node in nodes) {
    cpt = as.data.frame(bn[[node]]$prob)
    if(ncol(cpt) < 3) # BUG: whith 0 parents the variable name is missing
      names(cpt)[1] = node
    names(cpt)[ncol(cpt)] = paste("p", node, sep=".")
    cpt.table = merge(cpt.table, cpt, by=names(cpt)[-length(cpt)])
  }
  
  joint.table = cbind(cpt.table[, nodes, drop=FALSE], rep(1, nrow(cpt.table)))
  names(joint.table) = c(names(joint.table)[1:nbnodes], "p")
  for(i in 1:nbnodes) {
    joint.table[, nbnodes + 1] = joint.table[, nbnodes + 1] * cpt.table[, nbnodes + i]
  }
  
  return(joint.table)
  
}#EXACT.CPT

exact.prob = function(fitted, event, evidence, debug = FALSE){
  
  cpt = exact.cpt(fitted, event, evidence, debug)
  
  m = nrow(cpt)
  
  # evaluate the expression defining the event.
  if (identical(event, TRUE))
    r = rep(TRUE, m)
  else
    r = eval(event, cpt, parent.frame())
  
  # double check that this is a logical vector.
  if (!is.logical(r))
    stop("event must evaluate to a logical vector.")
  # double check that it is of the right length.
  if (length(r) != m)
    stop("logical vector for event is of length ", length(r), " instead of ", m, ".")
  
  numerator = sum(cpt[r, ncol(cpt)])
  denominator = sum(cpt[, ncol(cpt)])
  
  if(denominator == 0)
    return(0)
  
  return(numerator/denominator)
  
}#EXACT.PROB

exact.dist = function(fitted, nodes, evidence, n, batch, debug = FALSE) {
  
  # TODO
  stop("Not supported yet.")
  
}#EXACT.DIST

logic.sampling = function(fitted, event, evidence, n, batch, debug = FALSE) {
  
  cpxe = cpe = 0L
  filtered = logical(n)
  matching = logical(n)
  r = logical(n)
  
  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch
  
  for (m in c(rep(batch, nbatches), last.one)) {
    
    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL
    
    # generate random data from the bayesian network.
    if (m > 0) {
      
      if (is.fitted.discrete(fitted))
        generated.data = rbn.discrete(x = fitted, n = m)
      else
        generated.data = rbn.continuous(x = fitted, n = m)
      
    }#THEN
    else
      break
    
    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")
    
    # evaluate the expression defining the evidence.
    if (identical(evidence, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(evidence, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("evidence must evaluate to a logical vector.")
    # double check that it is of the right length.
    if (length(r) != m)
      stop("logical vector for evidence is of length ", length(r),
           " instead of ", m, ".")
    # filter out the samples not matching the evidence we assume.
    filtered = r & !is.na(r)
    
    # update the global counters.
    cpe = cpe + length(which(filtered))
    
    if (debug) {
      
      lwfilter = length(which(filtered))
      if (!identical(evidence, TRUE))
        cat("  > evidence matches ", lwfilter, " samples out of ", m,
            " (p = ", lwfilter/m, ").\n", sep = "")
      else
        cat("  > evidence matches ", m, " samples out of ", m,
            " (p = 1).\n", sep = "")
      
    }#THEN
    
    # evaluate the expression defining the event.
    if (identical(event, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(event, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("event must evaluate to a logical vector.")
    # double check that it is of the right length.
    if (length(r) != m)
      stop("logical vector for event is of length ", length(r),
           " instead of ", m, ".")
    # filter out the samples not matching the event we are looking for.
    matching = filtered & r & !is.na(r)
    
    # update the global counters.
    cpxe = cpxe + length(which(matching))
    
    if (debug) {
      
      lwmatch = length(which(matching))
      lwratio = ifelse(lwfilter == 0, 0, lwmatch/lwfilter)
      if (!identical(event, TRUE))
        cat("  > event matches ", lwmatch, " samples out of ", lwfilter,
            " (p = ", lwratio, ").\n", sep = "")
      else
        cat("  > event matches ", lwfilter, " samples out of ", lwfilter,
            " (p = 1).\n", sep = "")
      
    }#THEN
    
  }#FOR
  
  # prevent divide-by-zero errors.
  result = ifelse(cpe == 0, 0, cpxe / cpe)
  
  if (debug && (nbatches > 1)) {
    
    cat("* generated a grand total of", n, "samples.\n")
    cat("  > event matches ", cpxe, " samples out of ", cpe,
        " (p = ", result, ").\n", sep = "")
    
  }#THEN
  
  return(result)
  
}#LOGIC.SAMPLING

logic.distribution = function(fitted, nodes, evidence, n, batch, debug = FALSE) {
  
  filtered = logical(n)
  result = NULL
  
  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch
  
  for (m in c(rep(batch, nbatches), last.one)) {
    
    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL
    
    # generate random data from the bayesian network.
    if (m > 0) {
      
      if (is.fitted.discrete(fitted))
        generated.data = rbn.discrete(x = fitted, n = m)
      else
        generated.data = rbn.continuous(x = fitted, n = m)
      
    }#THEN
    else
      break
    
    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")
    
    # evaluate the expression defining the evidence.
    r = eval(evidence, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("evidence must evaluate to a logical vector.")
    # double check that it is of the right length.
    if ((length(r) != 1) && (length(r) != m))
      stop("logical vector for evidence is of length ", length(r),
           " instead of ", m, ".")
    # filter out the samples not matching the evidence we assume.
    filtered = r & !is.na(r)
    
    if (debug) {
      
      lwfilter = length(which(filtered))
      if (!identical(evidence, TRUE))
        cat("  > evidence matches ", lwfilter, " samples out of ", m,
            " (p = ", lwfilter/m, ").\n", sep = "")
      else
        cat("  > evidence matches ", m, " samples out of ", m,
            " (p = 1).\n", sep = "")
      
    }#THEN
    
    # update the return value.
    result = rbind(result, generated.data[filtered, nodes, drop = FALSE])
    
  }#FOR
  
  # reset the row names.
  rownames(result) = NULL
  
  if (debug && (nbatches > 1)) 
    cat("* generated a grand total of", n, "samples.\n")
  
  return(result)
  
}#LOGIC.DISTRIBUTION
