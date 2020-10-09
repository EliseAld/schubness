#####
# MISSING NETDIST FUNCTIONS
#####
transfmat <- function(x) {
  Adj <- x
  n <- ncol(Adj)
  tag <- "undir"
  if (!isSymmetric(x,check.attributes=FALSE, check.names=FALSE)){
    zero <- matrix(0, nrow=n, ncol=n)
    tmp <- rbind(cbind(zero,t(Adj)),cbind(Adj,zero))
    Adj <- tmp
    tag <- "dir"
    n <- ncol(tmp)
  }
  if (any(Adj > 1) || any(Adj < 0)){
    warning("Edge weight should be >= 0 and <= 1, scaling has been automatically applied!", call.=FALSE)
    Adj <- (Adj - min(Adj)) / (max(Adj) - min(Adj))
  }
  return(list(adj = Adj, tag = tag, N=n))
}
Lap <- function(x){
  L <- -x
  diag(L) <- apply(x,2,sum)
  return(L)
}
ipsen_minus_one <- function(g,n){
  return(sqrt(
    1/(pi*g) +
      1/(2*g*(atan(sqrt(n)/g)+pi/2)**2)*(pi/2+ (sqrt(n)/g)/(1+(sqrt(n)/g)**2)+atan(sqrt(n)/g))
    -4*(pi-(g/sqrt(n))*log(1/(1+(sqrt(n)/g)**2))+atan(sqrt(n)/g))/
      (pi*g*(4+(sqrt(n)/g)**2)*(atan(sqrt(n)/g)+pi/2)))-1)
}
ipsen_minus_one_dir  <- function(g,n){
  return(ZZ(g)^2*MM(0,g)+WW(n,g)^2*MM(n-2,g)+WW(n,g)^2*MM(n,g)+WWp(n,g)^2*MM(2*n-2,g)
         -2*ZZ(g)*WW(n,g)*LL(0,n-2,g)-2*ZZ(g)*WW(n,g)*LL(0,n,g)-2*ZZ(g)*WWp(n,g)*LL(0,2*n-2,g)
         +2*WW(n,g)*WW(n,g)*LL(n-2,n,g) +2*WW(n,g)*WWp(n,g)*LL(n-2,2*n-2,g) +2*WW(n,g)*WWp(n,g)*LL(n,2*n-2,g)-1)
}
ipsen <- function(object, ga=NULL){
  if (is.null(ga)){
    if (object$tag == "undir"){
      optgamma <- uniroot(ipsen_minus_one,
                          c(0.01,1),
                          n=object$N,
                          maxiter=100000,
                          tol=.Machine$double.eps)$root
    }
    else {
      optgamma <- uniroot(ipsen_minus_one_dir,
                          c(0.01,1),
                          n=object$N,
                          maxiter=100000,
                          tol=.Machine$double.eps)$root
    }
  }
  else {
    optgamma <- ga
  }
  laplist <- object$L
  n.cores <- NULL
  if (!is.na(match("n.cores",names(list())))) {
    n.cores <- list()[["n.cores"]]
  }
  verbose <- FALSE
  if (!is.na(match("verbose", names(list())))) {
    verbose <- list()[["verbose"]]
  }
  if(detectCores() >= 2 && (is.null(n.cores) || n.cores>1)){
    if (is.null(n.cores) || n.cores >= detectCores()){
      if (length(laplist) < detectCores()){
        n.cores <- length(laplist)
      }
      else {
        n.cores <- detectCores() - 1
      }
    }
    cl <- makeCluster(n.cores)
    clusterEvalQ(cl,
                 {K <- function(mygamma,given_omega) {
      return(1/integrate(lorentz,lower=0,upper=Inf,mygamma=mygamma,given_omega=given_omega, stop.on.error = FALSE)$value)
      }
      rho <- function(omega, mygamma, ll) {
        ll[[2]]*lorentz(omega,mygamma,ll[[1]])
      }
      lorentz <- function(omega,mygamma,given_omega){
        l <-0
        for (i in 2:length(given_omega)) {
          l = l + mygamma/( (omega-given_omega[i])**2+mygamma**2)                          }
        return(l)
      }
    })
    if (verbose)
      cat("Start computing eigenvalues with multiple cores\n")
    ll <- clusterApply(cl,laplist,function(x,mygamma=optgamma){
      myomega <- sqrt(abs(round(spec(x),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    })
    stopCluster(cl)
  } 
  else {
    ll <- lapply(1:length(laplist),function(x,mygamma,laplist){
      if (verbose)
        cat("Done",x,"/",length(laplist),"\n")
      aa <- laplist[[x]]
      myomega <- sqrt(abs(round(spec(aa),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    }, mygamma=optgamma, laplist=laplist)
  }
  mydistfun <- function(a,b, optgamma){
    integrand <- function(omega, mygamma, given_omega_G, given_omega_H){
      (rho(omega, optgamma,a)-rho(omega,optgamma,b))**2
    }
    tmp <- sqrt(integrate(integrand,lower=0,upper=Inf,mygamma=optgamma,given_omega_G=a[[1]],given_omega_H=b[[1]], stop.on.error=FALSE,rel.tol=.Machine$double.eps,subdivisions=1e4)$value)
    return(tmp)
  }
  if (verbose)
    cat("Start computing mutual distances\n")
  if (length(laplist) == 2){
    dist <- mydistfun(ll[[1]], ll[[2]], optgamma=optgamma)
    names(dist) <- "IM"
  } else {
    idx <- combn(length(ll),2)
    tmpdist <- sapply(1:dim(idx)[2], function(x,ll,optgamma, idx){
      if (verbose)
        cat("D(",idx[1,x],",", idx[2,x],")\n")
      mydistfun(ll[[idx[1,x]]], ll[[idx[2,x]]], optgamma)
    }, ll=ll, optgamma=optgamma, idx=idx)
    dist <- matrix(NA,ncol=length(ll), nrow=length(ll))
    dist[t(idx)] <- dist[t(idx)[,c(2,1)]] <- tmpdist
    diag(dist) <- 0
  }
  return(dist)
}
him <- function(object, ga=NULL, components=TRUE, ltag=FALSE, rho=1){
  ipd <- ipsen(object$LAP, ga)
  had <- hamming(object$ADJ)
  gloc <- sqrt(1/(1+rho)) * sqrt(had**2+ rho*(ipd**2))
  if (length(gloc)==1)
    names(gloc) <- "HIM"
  if(components==TRUE){
    if (ltag){
      dist <- list(H=had,IM=ipd,HIM=gloc)
    } else {
      dist <- c(had, ipd, gloc)
      names(dist) <- c("H","IM","HIM")
    }
  } else {
    dist <- gloc
    if (!ltag)
      names(dist) <- "HIM"
  }
  return(dist)
}

#####
# MISSING DYNO FUNCTIONS
#####
add_cluster_graph2<-function (dataset, milestone_network, grouping = NULL, explicit_splits = FALSE) {
  testthat::expect_true(is_data_wrapper(dataset))
  if (is.null(grouping)) {
    testthat::expect_true(is_wrapper_with_grouping(dataset))
  }
  else {
    dataset <- dataset %>% add_grouping(grouping)
  }
  grouping <- get_grouping(dataset)
  grouping <- grouping[!is.na(grouping)]
  milestone_ids <- unique(c(milestone_network$to, milestone_network$from))
  check_milestone_network2(milestone_ids, milestone_network, allow_self_loops = T)
  if (explicit_splits) {
    milestone_network <- cluster_graph_add_explicit_splits(milestone_network)
    milestone_ids <- unique(c(milestone_network$to, milestone_network$from))
  }
  both_directions <- bind_rows(milestone_network %>% select(from, 
                                                            to) %>% mutate(label = from, percentage = 0), milestone_network %>% 
                                 select(from, to) %>% mutate(label = to, percentage = 1))
  progressions <- tibble(cell_id = names(grouping), label = grouping) %>% 
    left_join(both_directions, by = "label") %>% group_by(cell_id) %>% 
    arrange(desc(percentage)) %>% slice(1) %>% ungroup() %>% 
    select(-label)
  add_trajectory(dataset = dataset, milestone_ids = milestone_ids, 
                 milestone_network = milestone_network, divergence_regions = NULL, 
                 progressions = progressions, allow_self_loops = T)
}
check_milestone_network2 <- function(milestone_ids, milestone_network, allow_self_loops = T) {
  assert_that(is.data.frame(milestone_network),
              ncol(milestone_network) == 4,
              setequal(colnames(milestone_network), c("from", "to", "length", "directed")))
  milestone_network <- milestone_network %>% select(from, to, length, directed)
  
  assert_that(is.character(milestone_network$from),
              is.character(milestone_network$to),
              is.numeric(milestone_network$length),
              is.logical(milestone_network$directed),
              all(milestone_network$from %in% milestone_ids),
              all(milestone_network$to %in% milestone_ids),
              !any(duplicated(milestone_network %>% select(from, to))))
  if (!allow_self_loops) {
    assert_that(!any((milestone_network$from == milestone_network$to) & milestone_network$length > 0))
  }
  check1 <- milestone_network
  if (allow_self_loops) check1 <- check1 %>% filter(from != to)
  check <-
    inner_join(
      check1 %>% transmute(from, to, left = "left"),
      check1 %>% transmute(from = to, to = from, right = "right"),
      by = c("from", "to")
    )
  assert_that(nrow(check) == 0, msg = "Milestone network should not contain A->B B->A edges")
  
  milestone_network
}
calculate_metrics2 <- function(model1,model2,metric_ids,expression_source=model1$expression) {
  if (!all(sapply(seq_along(metric_ids), function(i) !is.function(metric_ids[[i]]) || 
                  !is.null(names(metric_ids)[[i]])))) {
    stop("All custom metrics (functions) must be named!")
  }
  valid_metrics <- dyneval::metrics$metric_id
  character_metrics <- as.character(keep(metric_ids, is.character))
  if (!all(character_metrics %in% valid_metrics)) {
    stop("Invalid metrics: ", glue::glue_collapse(setdiff(character_metrics, 
                                                          valid_metrics), ", "))
  }
  summary_list <- list()
  model1 <- dynwrap::simplify_trajectory(model1, allow_self_loops = T)
  if (!is.null(model2)) {
    model2 <- dynwrap::simplify_trajectory(model2, allow_self_loops = T)
  }
  if ("correlation" %in% metric_ids) {
    testthat::expect_true(dynwrap::is_wrapper_with_waypoint_cells(model1))
    testthat::expect_true(is.null(model2) || dynwrap::is_wrapper_with_waypoint_cells(model2))
    if (!is.null(model2)) {
      testthat::expect_true(all(model2$cell_ids %in% model1$cell_ids))
      model2$cell_ids <- model1$cell_ids
      waypoints <- unique(c(model1$waypoint_cells, model2$waypoint_cells))
      time0 <- Sys.time()
      model1$geodesic_dist <- dynwrap::calculate_geodesic_distances(model1, 
                                                                    waypoints)
      model2$geodesic_dist <- dynwrap::calculate_geodesic_distances(model2, 
                                                                    waypoints)
      time1 <- Sys.time()
      summary_list$time_waypointedgeodesic <- as.numeric(difftime(time1, 
                                                                  time0, units = "sec"))
    }
    if (!is.null(model2)) {
      model1$geodesic_dist[is.infinite(model1$geodesic_dist)] <- .Machine$double.xmax
      model2$geodesic_dist[is.infinite(model2$geodesic_dist)] <- .Machine$double.xmax
      testthat::expect_equal(rownames(model1$geodesic_dist), 
                             rownames(model2$geodesic_dist))
      testthat::expect_equal(colnames(model1$geodesic_dist), 
                             colnames(model2$geodesic_dist))
      time0 <- Sys.time()
      if (length(unique(c(model2$geodesic_dist))) == 1 || 
          length(unique(c(model1$geodesic_dist))) == 1) {
        summary_list$correlation <- 0
      }
      else {
        summary_list$correlation <- cor(model1$geodesic_dist %>% 
                                          as.vector, model2$geodesic_dist %>% as.vector, 
                                        method = "spearman") %>% max(0)
      }
      time1 <- Sys.time()
      summary_list$time_correlation <- as.numeric(difftime(time1, 
                                                           time0, units = "sec"))
    }
    else {
      summary_list <- c(summary_list, list(correlation = 0))
    }
  }
  if ("edge_flip" %in% metric_ids) {
    if (!is.null(model2)) {
      net1 <- model2$milestone_network
      net2 <- model1$milestone_network
      time0 <- Sys.time()
      summary_list$edge_flip <- calculate_edge_flip(net1, 
                                                    net2)
      time1 <- Sys.time()
      summary_list$time_edge_flip <- as.numeric(difftime(time1, 
                                                         time0, units = "sec"))
    }
    else {
      summary_list$edge_flip <- 0
    }
  }
  if ("him" %in% metric_ids) {
    if (!is.null(model2)) {
      net1 <- model2$milestone_network
      net2 <- model1$milestone_network
      time0 <- Sys.time()
      summary_list$him <- calculate_him(net1, net2)
      time1 <- Sys.time()
      summary_list$time_him <- as.numeric(difftime(time1, 
                                                   time0, units = "sec"))
    }
    else {
      summary_list$him <- 0
    }
  }
  if ("isomorphic" %in% metric_ids) {
    if (!is.null(model2)) {
      graph1 <- model2$milestone_network %>% igraph::graph_from_data_frame()
      graph2 <- model1$milestone_network %>% igraph::graph_from_data_frame()
      time0 <- Sys.time()
      summary_list$isomorphic <- as.numeric(igraph::isomorphic(graph1, 
                                                               graph2))
      time1 <- Sys.time()
      summary_list$time_isomorphic <- as.numeric(difftime(time1, 
                                                          time0, units = "sec"))
    }
    else {
      summary_list$isomorphic <- 0
    }
  }
  if (any(c("rf_mse", "rf_rsq", "rf_nmse", "lm_mse", "lm_rsq", 
            "lm_nmse") %in% metric_ids)) {
    time0 <- Sys.time()
    position_predict <- calculate_position_predict(model1, 
                                                   model2)
    time1 <- Sys.time()
    summary_list$time_pp <- as.numeric(difftime(time1, time0, 
                                                units = "sec"))
    summary_list <- c(summary_list, position_predict$summary[intersect(metric_ids, 
                                                                       names(position_predict$summary))])
  }
  if (any(c("featureimp_cor", "featureimp_wcor") %in% metric_ids)) {
    time0 <- Sys.time()
    featureimp <- calculate_featureimp_cor(model1, model2, 
                                           expression_source = expression_source)
    time1 <- Sys.time()
    summary_list$time_featureimp <- as.numeric(difftime(time1, 
                                                        time0, units = "sec"))
    summary_list$featureimp_cor <- featureimp$featureimp_cor
    summary_list$featureimp_wcor <- featureimp$featureimp_wcor
  }
  if (any(c("featureimp_ks", "featureimp_wilcox") %in% metric_ids)) {
    time0 <- Sys.time()
    featureimp <- calculate_featureimp_enrichment(model1, 
                                                  model2, expression_source = expression_source)
    time1 <- Sys.time()
    summary_list$time_featureimp_enrichment <- as.numeric(difftime(time1, 
                                                                   time0, units = "sec"))
    summary_list$featureimp_ks <- featureimp$featureimp_ks
    summary_list$featureimp_wilcox <- featureimp$featureimp_wilcox
  }
  if (any(c("recovery_branches", "relevance_branches", "F1_branches", 
            "recovery_milestones", "relevance_milestones", "F1_milestones") %in% 
          metric_ids)) {
    model1_simplified <- dynwrap::simplify_trajectory(model1, 
                                                      allow_self_loops = TRUE)
    if (!is.null(model2)) {
      model2_simplified <- dynwrap::simplify_trajectory(model2, 
                                                        allow_self_loops = TRUE)
    }
    else {
      model2_simplified <- NULL
    }
    if (any(c("recovery_branches", "relevance_branches", 
              "F1_branches") %in% metric_ids)) {
      time0 <- Sys.time()
      mapping_branches <- calculate_mapping_branches2(model1_simplified, 
                                                     model2_simplified)
      time1 <- Sys.time()
      summary_list$time_mapping_branches <- as.numeric(difftime(time1, 
                                                                time0, units = "sec"))
      summary_list <- c(summary_list, mapping_branches)
    }
    if (any(c("recovery_milestones", "relevance_milestones", 
              "F1_milestones") %in% metric_ids)) {
      time0 <- Sys.time()
      mapping_milestones <- calculate_mapping_milestones2(model1_simplified, 
                                                         model2_simplified)
      time1 <- Sys.time()
      summary_list$time_mapping_milestones <- as.numeric(difftime(time1, 
                                                                  time0, units = "sec"))
      summary_list <- c(summary_list, mapping_milestones)
    }
  }
  for (i in seq_along(metric_ids)) {
    f <- metric_ids[[i]]
    fn <- names(metric_ids)[[i]]
    if (is.function(f)) {
      if (!is.null(model2)) {
        time0 <- Sys.time()
        output <- f(model1, model2)
        time1 <- Sys.time()
        summary_list[[paste0("time_", fn)]] <- as.numeric(difftime(time1, 
                                                                   time0, units = "sec"))
        if (length(output) != 1) {
          stop("Metric ", sQuote(fn), " should return exactly 1 numeric score.")
        }
      }
      else {
        output <- 0
      }
      names(output) <- fn
      summary_list[names(output)] <- output
    }
  }
  summary <- as_tibble(summary_list)
  summary
}
calculate_mapping_branches2 <- function(dataset, prediction, simplify = TRUE) {
  mapping <- calculate_mapping2(dataset, prediction, grouping = "branches", simplify = simplify)
  names(mapping) <- paste0(names(mapping), "_branches")
  mapping
}
calculate_mapping_milestones2 <- function(dataset, prediction, simplify = TRUE) {
  mapping <- calculate_mapping(dataset, prediction, "milestones", simplify = simplify)
  names(mapping) <- paste0(names(mapping), "_milestones")
  mapping
}
calculate_mapping2 <- function(dataset, prediction, grouping = c("branches", "milestones"), simplify = TRUE) {
  grouping <- match.arg(grouping)
  
  if (is.null(prediction) || is.null(dataset)) {
    lst(
      recovery = 0,
      relevance = 0,
      F1 = 0
    )
  } else {
    if (simplify) {
      dataset <- dynwrap::simplify_trajectory(dataset, allow_self_loops = TRUE)
      prediction <- dynwrap::simplify_trajectory(prediction, allow_self_loops = TRUE)
    }
    
    if (grouping == "branches") {
      groups_dataset <- dataset %>% dynwrap::group_onto_trajectory_edges()
      groups_prediction <- prediction %>% dynwrap::group_onto_trajectory_edges()
    } else if (grouping == "milestones") {
      groups_dataset <- dataset %>% dynwrap::group_onto_nearest_milestones()
      groups_prediction <- prediction %>% dynwrap::group_onto_nearest_milestones()
    }
    
    groups_dataset <- groups_dataset %>% as.character() %>% enframe("cell_id", "group_dataset")
    groups_dataset_levels <- unique(groups_dataset$group_dataset) %>% na.omit()
    groups_prediction <- groups_prediction %>% as.character() %>% enframe("cell_id", "group_prediction")
    groups_prediction_levels <- unique(groups_prediction$group_prediction) %>% na.omit()
    
    groups <- full_join(groups_dataset, groups_prediction, "cell_id")
    
    # calculate the size of the intersections and of each group separately
    intersections <-
      groups %>%
      filter(!is.na(group_dataset), !is.na(group_prediction)) %>%
      group_by(group_dataset, group_prediction) %>%
      summarise(intersection = n()) %>%
      ungroup() %>%
      mutate(
        group_dataset = factor(group_dataset, levels = groups_dataset_levels),
        group_prediction = factor(group_prediction, levels = groups_prediction_levels)
      ) %>%
      complete(
        group_dataset, group_prediction,
        fill = list(intersection = 0)
      ) %>%
      mutate_if(is.factor, as.character)
    
    n_dataset <-
      groups %>%
      filter(!is.na(group_dataset)) %>%
      group_by(group_dataset) %>%
      summarise(n_dataset = n())
    
    n_prediction <-
      groups %>%
      filter(!is.na(group_prediction)) %>%
      group_by(group_prediction) %>%
      summarise(n_prediction = n())
    
    # now join and calculate the jaccard
    jaccards <- intersections %>%
      left_join(n_dataset, "group_dataset") %>%
      left_join(n_prediction, "group_prediction") %>%
      mutate(jaccard = intersection / (n_dataset + n_prediction - intersection))
    
    # calculate the recovery and relevance
    recoveries <- jaccards %>%
      group_by(group_dataset) %>%
      arrange(-jaccard) %>%
      slice(1) %>%
      ungroup()
    
    relevances <-
      jaccards %>%
      group_by(group_prediction) %>%
      arrange(-jaccard) %>%
      slice(1) %>%
      ungroup()
    
    # calculate the final scores
    zero_if_na <- function(x) if(is.na(x)) {0} else {x}
    
    lst(
      recovery = mean(recoveries$jaccard) %>% zero_if_na(),
      relevance = mean(relevances$jaccard) %>% zero_if_na(),
      F1 = calculate_harmonic_mean2(recovery, relevance)
    )
  }
}
calculate_harmonic_mean2 <- function(..., weights = NULL) {
  x <- process_combination_input2(...)
  if (is.null(weights)) {
    ncol(x) / rowSums(1/x)
  } else {
    sum(weights) / rowSums(process_weights(weights, nrow(x))/x)
  }
}
process_combination_input2 <- function(...) {
  dots <- list(...)
  if (length(dots) > 1 && all(map_lgl(dots, is.numeric))) {
    do.call(cbind, dots)
  } else if (is.list(..1) && all(map_lgl(..1, is.numeric))) {
    do.call(cbind, ..1)
  } else if (is.matrix(..1) && is.numeric(..1)) {
    ..1
  } else if (is.numeric(..1)) {
    do.call(cbind, as.list(..1))
  } else {
    stop("Invalid input")
  }
}
calculate_him2 <- function(net1, net2, simplify=T, d="HIM") {
  adjacencies <- get_matched_adjacencies(net1, net2, simplify = T)
  if (max(adjacencies[[2]]) == 0 || max(adjacencies[[1]]) == 0) {
    return(0)
  }
  
  g1 <- transfmat(adjacencies[[1]]/sum(adjacencies[[1]]))
  g2 <- transfmat(adjacencies[[2]]/sum(adjacencies[[2]]))
  myadj <- list("HIM",G=list(g1$adj,g2$adj),N=g1$N,tag=g1$tag)
  mylap <- list(L=list(Lap(g1$adj), Lap(g2$adj)),N=g1$N, tag=g1$tag)
  netdist <- him(list(ADJ=myadj,LAP=mylap),  ga=0.1,  components=T, ltag=FALSE)['HIM']
  netdist[netdist < 0] <- 0
  return(1 - netdist)
}
get_matched_adjacencies <- function(net1, net2, simplify = TRUE) {
  if (simplify) {
    directed1 <- any(net1$directed)
    directed2 <- any(net2$directed)
    net1 <- net1 %>%
      rename(weight = length) %>%
      filter(!(from == to & weight == 0)) %>% # remove self loop edges with length 0
      igraph::graph_from_data_frame(directed = F) %>%
      dynwrap::simplify_igraph_network() %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed1) %>%
      insert_two_nodes_into_selfloop() %>%
      change_single_edge_into_double() %>%
      insert_one_node_into_duplicate_edges()
    net2 <- net2 %>%
      rename(weight = length) %>%
      filter(!(from == to & weight == 0)) %>% # remove self loop edges with length 0
      igraph::graph_from_data_frame(directed = F) %>%
      dynwrap::simplify_igraph_network() %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed2) %>%
      insert_two_nodes_into_selfloop() %>%
      change_single_edge_into_double() %>%
      insert_one_node_into_duplicate_edges()
  }
  
  adj1 <- get_adjacency_lengths(net1)
  adj2 <- get_adjacency_lengths(net2)
  
  # make the adjacency matrices have the same dimensions
  if (nrow(adj1) > nrow(adj2)) {
    adj2 <- complete_matrix(adj2, nrow(adj1), fill = 0)
  } else {
    adj1 <- complete_matrix(adj1, nrow(adj2), fill = 0)
  }
  
  lst(adj1, adj2)
}

