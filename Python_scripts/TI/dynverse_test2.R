library(SeuratDisk)
library(Seurat)
library(dynwrap)
library(dyno)
library(dyneval)
library(tidyverse)
library(assertthat)
library(ggpubr)
library(pbapply)
library(BBmisc)
library(Matrix)
library(FedData)
#https://zenodo.org/record/1443566
path <- "/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/dynverse/"
#dataset <- datasets[311:337]
#dataRDS <- pblapply(dataset,
#                    function(x) readRDS(paste0(path,'/1_rds/',x)))
#Order <- pblapply(dataRDS,
#                  function(x) {Order=x$prior_information$groups_id;
#                  rownames(Order)<-Order$cell_id;
#                  Order<-Order[x$cell_ids,];
#                  return(Order$group_id)})
#data <- pblapply(seq(length(dataset)),
#                   function(x) dataRDS[[x]]$counts)
#pblapply(seq(length(data)),
#         function(x) {write.csv(data[[x]],
#                                file=paste0(path,
#                                            strsplit(dataset[x],"_")[[1]][1],"/csv/",
#                                            strsplit(dataset[x],"[.]")[[1]][1],
#                                            "_data.csv"));
#                      write.csv(Order[[x]],
#                                file=paste0(path,
#                                            strsplit(dataset[x],"_")[[1]][1],"/csv/",
#                                            strsplit(dataset[x],"[.]")[[1]][1],
#                                            "_order.csv"))})

source("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/dynverse_test_functions.R")
methods <- c("nothing","mp_normal","ls","ls_nicdm","dsl")
datasets_choice = unname(sapply(c('dyntoy_looping_1.rds', 'silver_planaria-full_plass.rds', 'silver_mouse-cell-atlas-combination-2.rds',
                                  'dyntoy_linear_2.rds', 'gold_hematopoiesis-gates_olsson.rds', 'silver_planaria-pair-10_plass.rds',
                                  'gold_psc-astrocyte-maturation-glia_sloan.rds', 'dyntoy_linear_3.rds', 'prosstt_multifurcating_6.rds',
                                  'silver_mouse-cell-atlas-combination-3.rds', 'dyngen_2.rds', 'silver_planaria-combination-1_plass.rds',
                                  'gold_germline-human-female-weeks_li.rds', 'dyntoy_looping_2.rds',
                                  'silver_fetal-liver-fetal-hematopoiesis_mca.rds', 'gold_stimulated-dendritic-cells-LPS_shalek.rds',
                                  'gold_germline-human-female_guo.rds', 'dyntoy_disconnected_10.rds',
                                  'silver_mouse-cell-atlas-combination-1.rds', 'dyntoy_binary_tree_7.rds',
                                  'silver_planaria-combination-13_plass.rds', 'dyntoy_linear_1.rds', 'silver_planaria-pair-7_plass.rds',
                                  'dyntoy_binary_tree_6.rds', 'dyngen_1.rds', 'silver_hepatoblast-differentiation_yang.rds',
                                  'splatter_multifurcating_10.rds', 'silver_placenta-trophoblast-differentiation-invasive_mca.rds',
                                  'dyntoy_looping_3.rds', 'dyntoy_looping_7.rds',
                                  'silver_oligodendrocyte-differentiation-subclusters_marques.rds',
                                  'silver_trophoblast-stem-cell-trophoblast-differentiation_mca.rds',
                                  'silver_mouse-cell-atlas-combination-4.rds', 'dyntoy_binary_tree_2.rds',
                                  'silver_distal-lung-epithelium_treutlein.rds', 'dyngen_10.rds',
                                  'silver_planaria-combination-7_plass.rds', 'silver_trophectoderm-monkey_nakamura.rds',
                                  'dyntoy_binary_tree_3.rds', 'silver_mouse-cell-atlas-combination-5.rds',
                                  'silver_epiblast-monkey_nakamura.rds', 'silver_planaria-neuron-differentiation_plass.rds',
                                  'gold_stimulated-dendritic-cells-PAM_shalek.rds', 'silver_epidermis-hair-spatial_joost.rds',
                                  'dyntoy_looping_6.rds', 'silver_olfactory-projection-neurons_horns.rds',
                                  'silver_mouse-cell-atlas-combination-7.rds', 'dyntoy_binary_tree_1.rds', 'dyntoy_linear_7.rds',
                                  'silver_oligodendrocyte-differentiation-clusters_marques.rds', 'dyngen_12.rds', 'dyntoy_linear_6.rds',
                                  'prosstt_multifurcating_3.rds', 'silver_mouse-cell-atlas-combination-6.rds',
                                  'silver_planaria-pair-1_plass.rds', 'dyntoy_bifurcating_10.rds',
                                  'silver_bone-marrow-mesenchyme-erythrocyte-differentiation_mca.rds',
                                  'gold_mESC-differentiation_hayashi.rds', 'gold_developing--cells_schlitzer.rds',
                                  'silver_planaria-parenchyme-differentiation_plass.rds', 'dyntoy_diverging_converging_10.rds',
                                  'dyntoy_multifurcating_10.rds', 'dyngen_61.rds', 'prosstt_bifurcating_1.rds', 'dyngen_60.rds',
                                  'dyngen_110.rds', 'dyntoy_linear_10.rds', 'silver_germline-human-female_li.rds',
                                  'silver_cell-cycle_leng.rds', 'prosstt_tree_1.rds', 'silver_planaria-pair-11_plass.rds',
                                  'prosstt_linear_6.rds', 'splatter_bifurcating_7.rds', 'splatter_bifurcating_5.rds',
                                  'dyntoy_converging_6.rds', 'splatter_bifurcating_10.rds', 'dyntoy_tree_10.rds',
                                  'silver_planaria-combination-12_plass.rds', 'silver_planaria-pair-6_plass.rds',
                                  'dyntoy_converging_10.rds', 'prosstt_tree_3.rds', 'dyntoy_cyclic_7.rds', 'dyngen_62.rds',
                                  'prosstt_bifurcating_2.rds', 'prosstt_bifurcating_3.rds',
                                  'silver_olfactory-projection-neurons-DA1_horns.rds', 'dyntoy_cyclic_6.rds', 'dyntoy_converging_7.rds',
                                  'splatter_bifurcating_4.rds', 'dyntoy_converging_3.rds', 'silver_planaria-combination-6_plass.rds',
                                  'prosstt_linear_1.rds', 'prosstt_tree_6.rds', 'silver_thymus-t-cell-differentiation_mca.rds',
                                  'dyntoy_cyclic_2.rds', 'prosstt_bifurcating_7.rds', 'prosstt_bifurcating_6.rds', 'dyntoy_cyclic_3.rds',
                                  'dyngen_72.rds', 'dyngen_102.rds', 'splatter_binary_tree_10.rds', 'dyntoy_converging_2.rds',
                                  'dyntoy_looping_10.rds', 'dyntoy_connected_10.rds', 'splatter_bifurcating_1.rds',
                                  'gold_pancreatic-beta-cell-maturation_zhang.rds', 'splatter_bifurcating_3.rds',
                                  'splatter_tree_10.rds', 'dyngen_70.rds', 'dyngen_100.rds',
                                  'dyntoy_cyclic_1.rds', 'prosstt_bifurcating_4.rds', 'silver_planaria-combination-14_plass.rds',
                                  'prosstt_bifurcating_5.rds', 'dyngen_65.rds', 'dyngen_101.rds', 'dyngen_71.rds',
                                  'dyntoy_binary_tree_10.rds', 'prosstt_linear_3.rds', 'dyntoy_converging_1.rds',
                                  'gold_stimulated-dendritic-cells-PIC_shalek.rds', 'splatter_bifurcating_2.rds',
                                  'gold_human-embryos_petropoulos.rds', 'silver_ICM-monkey_nakamura.rds',
                                  'silver_germline-human-male_li.rds', 'silver_planaria-pair-12_plass.rds',
                                  'dyntoy_diverging_with_loops_1.rds', 'silver_neonatal-inner-ear-SC-HC_burns.rds',
                                  'silver_planaria-combination-3_plass.rds', 'prosstt_binary_tree_1.rds', 'dyntoy_multifurcating_3.rds',
                                  'dyntoy_tree_2.rds', 'dyngen_40.rds', 'prosstt_bifurcating_8.rds', 'dyngen_82.rds', 'dyntoy_tree_3.rds',
                                  'dyntoy_multifurcating_2.rds', 'dyntoy_cyclic_10.rds', 'silver_kidney-collecting-duct-clusters_park.rds',
                                  'splatter_binary_tree_7.rds', 'silver_neonatal-inner-ear-all_burns.rds', 'silver_planaria-pair-8_plass.rds',
                                  'splatter_binary_tree_5.rds', 'dyntoy_diverging_with_loops_2.rds',
                                  'silver_planaria-epidermis-differentiation_plass.rds',
                                  'silver_kidney-distal-convoluted-tubule_mca.rds', 'dyntoy_tree_1.rds',
                                  'silver_epidermis-hair-uHF_joost.rds', 'dyngen_80.rds', 'dyngen_95.rds', 'dyngen_42.rds',
                                  'prosstt_binary_tree_3.rds', 'silver_neonatal-rib-cartilage_mca.rds',
                                  'gold_macrophage-salmonella_saliba.rds', 'silver_fibroblast-reprogramming_treutlein.rds',
                                  'silver_hematopoiesis-clusters_olsson.rds', 'silver_planaria-pharynx-differentiation_plass.rds',
                                  'silver_olfactory-projection-neurons-DC3_VA1d_horns.rds', 'dyntoy_diverging_with_loops_3.rds',
                                  'silver_planaria-pair-5_plass.rds', 'splatter_binary_tree_4.rds', 'silver_planaria-combination-11_plass.rds',
                                  'dyntoy_diverging_with_loops_7.rds', 'silver_blastocyst-monkey_nakamura.rds', 'dyngen_52.rds', 'dyngen_91.rds',
                                  'dyngen_90.rds', 'prosstt_binary_tree_6.rds', 'silver_dentate-gyrus-neurogenesis_hochgerner.rds',
                                  'silver_planaria-pair-14_plass.rds', 'silver_planaria-muscle-differentiation_plass.rds',
                                  'dyntoy_diverging_with_loops_6.rds', 'splatter_binary_tree_1.rds',
                                  'silver_embryonic-mesenchyme-stromal-cell-cxcl14-cxcl12-axis_mca.rds',
                                  'silver_planaria-combination-5_plass.rds', 'gold_mesoderm-development_loh.rds',
                                  'splatter_binary_tree_3.rds', 'silver_planaria-pair-3_plass.rds', 'dyntoy_multifurcating_6.rds',
                                  'dyntoy_tree_7.rds', 'dyngen_92.rds', 'dyntoy_tree_6.rds', 'dyntoy_multifurcating_7.rds', 'dyngen_50.rds',
                                  'silver_kidney-collecting-duct-subclusters_park.rds', 'silver_kidney-bursh-border-to-s1_mca.rds',
                                  'splatter_binary_tree_2.rds', 'silver_planaria-combination-8_plass.rds',
                                  'silver_planaria-combination-2_plass.rds', 'splatter_tree_7.rds', 'dyntoy_diverging_with_loops_10.rds',
                                  'silver_mouse-cell-atlas-combination-10.rds', 'silver_planaria-pair-13_plass.rds', 'splatter_linear_5.rds',
                                  'silver_neonatal-inner-ear-TEC-SC_burns.rds', 'dyngen_22.rds', 'splatter_linear_4.rds',
                                  'silver_cortical-interneuron-differentiation_frazer.rds', 'gold_germline-human-male_guo.rds',
                                  'silver_planaria-pair-9_plass.rds', 'gold_germline-human-both_guo.rds', 'splatter_multifurcating_1.rds',
                                  'splatter_multifurcating_3.rds', 'dyntoy_diverging_converging_6.rds', 'splatter_tree_4.rds',
                                  'silver_epidermis-hair-IFE_joost.rds', 'dyntoy_bifurcating_7.rds', 'gold_myoblast-differentiation_trapnell.rds',
                                  'silver_mouse-cell-atlas-combination-8.rds', 'dyntoy_disconnected_7.rds', 'dyngen_20.rds', 'dyntoy_connected_7.rds',
                                  'dyntoy_connected_6.rds', 'silver_mammary-gland-involution-endothelial-cell-aqp1-gradient_mca.rds',
                                  'dyntoy_disconnected_6.rds', 'splatter_linear_7.rds', 'silver_planaria-pair-4_plass.rds', 'dyntoy_bifurcating_6.rds',
                                  'silver_planaria-combination-10_plass.rds', 'splatter_tree_5.rds', 'dyntoy_diverging_converging_7.rds',
                                  'silver_placenta-trophoblast-differentiation_mca.rds', 'splatter_multifurcating_2.rds',
                                  'gold_pancreatic-alpha-cell-maturation_zhang.rds', 'dyntoy_diverging_converging_3.rds', 'splatter_tree_1.rds',
                                  'dyntoy_bifurcating_2.rds', 'splatter_linear_3.rds', 'splatter_linear_10.rds', 'dyngen_31.rds',
                                  'dyntoy_disconnected_2.rds', 'dyntoy_connected_2.rds',
                                  'dyngen_30.rds', 'splatter_linear_2.rds', 'silver_planaria-combination-4_plass.rds',
                                  'silver_embronic-mesenchyme-neuron-differentiation_mca.rds', 'dyntoy_diverging_converging_2.rds',
                                  'splatter_multifurcating_7.rds', 'gold_NKT-differentiation_engel.rds', 'silver_planaria-pair-2_plass.rds',
                                  'splatter_multifurcating_5.rds', 'splatter_tree_2.rds',
                                  'silver_planaria-phagocyte-differentiation_plass.rds', 'dyngen_32.rds',
                                  'dyntoy_connected_1.rds', 'splatter_linear_1.rds', 'silver_planaria-combination-9_plass.rds',
                                  'silver_neonatal-inner-ear-TEC-HSC_burns.rds', 'gold_psc-astrocyte-maturation-neuron_sloan.rds',
                                  'splatter_tree_3.rds', 'dyntoy_diverging_converging_1.rds', 'splatter_multifurcating_4.rds'),
                                function(x) strsplit(x,".rds")[[1]]))
choice="dyntoy"
dataset_choice <- unname(sapply(sort(datasets_choice[grep(choice,datasets_choice)]),
                                function(x) paste0(strsplit(x,"_")[[1]][1],"/h5/",x)))
n_dim_choice <- c(25,50,100,500)
metric_choice <- c(#"him",
  "correlation","F1_branches","featureimp_wcor")
#pblapply(dataset_choice,
#         function(dts) lapply(n_dim_choice,
#                              function(dim) lapply(methods,
#                                                   function(mtd) {if (file.exists(paste0(path,dts,"_",dim,"dims_",mtd[1],"2.h5ad"))) {
#                                                                     Convert(paste0(path,dts,"_",dim,"dims_",mtd,"2.h5ad"),
#                                                                             dest = "h5seurat", overwrite = TRUE, verbose = F)}})))
#sapply(dataset_choice,function(x) {seurat <- LoadH5Seurat(paste0(path,x,"_25dims_nothing2.h5seurat"),verbose=F);
#tmp<-ifelse(ncol(seurat)==length(seurat@misc$Order),"",paste0("remove",x));
#tmp<-tmp[tmp!=""];
#print(unname(tmp))}
#) # remove datasets with less order than cells
make_paga_traj <- function(dataset_id, n_dim) {
  if (file.exists(paste0(path,dataset_id, "_",n_dim,"dims_",methods[1],"2.h5seurat"))) {
    seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id,"_",n_dim,"dims_", x,"2.h5seurat"),
                                              verbose=F))
    seurat <- lapply(seq(length(methods)),
                     function(x) {seurat[[x]]@misc$paga <- seurat[[x]]@misc$paga$connectivities_tree;
                     return(seurat[[x]])})
    paga <- lapply(seurat,
                   function(x) x@misc$paga)
    names(seurat) <- methods
    names(paga) <- methods
    n_col <- lapply(paga,
                    function(pg) sapply(seq(ncol(pg)),
                                        function(x) pg@p[x+1] - pg@p[x]))
    from <- lapply(n_col,
                   function(y) unlist(sapply(seq(y), function(x) rep(as.character(x),y[x]))))
    to <- lapply(paga,
                 function(pg) as.character(pg@i+1))
    missing_cluster_ids <- lapply(seq(length(seurat)),
                                  function(x) setdiff(as.character(as.numeric(seurat[[x]]$leiden)),
                                                      c(from[[x]],to[[x]])))
    from <- lapply(seq(length(from)),
                   function(x) c(from[[x]],missing_cluster_ids[[x]]))
    to <- lapply(seq(length(to)),
                 function(x) c(to[[x]],missing_cluster_ids[[x]]))
    length <- lapply(seq(length(methods)),
                     function(pg) c(paga[[pg]]@x,rep(0,length(missing_cluster_ids[[pg]]))))
    directed <- lapply(length,
                       function(x) rep(T, length(x)))
    milestone_network <- lapply(seq(length(methods)),
                                function(x) tibble("from"=from[[x]],
                                                   "to"=to[[x]],
                                                   "length"=length[[x]],
                                                   "directed"=directed[[x]]))
    dataset <- wrap_expression(counts = t(seurat$nothing@assays$RNA@counts),
                               expression = t(seurat$nothing@assays$RNA@data))
    traj <- lapply(seq(length(methods)),
                   function(x) {dataset %>% add_grouping(as.character(as.numeric(seurat[[x]]$leiden))) %>% add_cluster_graph2(milestone_network[[x]])})
    traj <- lapply(seq(length(methods)),
                   function(x) traj[[x]] %>% add_dimred(dyndimred::dimred_mds,
                                                        expression_source = dataset$expression))
    return(list("seurat"=seurat,
                "trajectory"=traj,
                "milestone_network"=milestone_network,
                "dataset_id"=dataset_id,
                "n_dim"=n_dim))
  }
  else {
  }
}
plot_paga_traj <- function(paga_traj_output) {
  idx_to_plot=seq(length(methods))
  if (paga_traj_output$dataset_id == 'GSE45719' & paga_traj_output$n_dim==50) {
    idx_to_plot=c(1,2,3,5)
  }
  pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/",
                  paga_traj_output$dataset_id, "/",
                  paga_traj_output$n_dim,"_dims_paga_graph.pdf" ))
  ggarrange(plotlist=lapply(idx_to_plot,
                            function(x) plot_dimred(paga_traj_output$trajectory[[x]], 
                                                    expression_source = t(paga_traj_output$seurat[[x]]@assays$RNA@data))),
            labels=methods)
  dev.off()
}
plot_compare_paga_traj <- function(paga_traj_output) {
  dataset <- wrap_expression(counts = t(paga_traj_output$seurat$nothing@assays$RNA@counts),
                             expression = t(paga_traj_output$seurat$nothing@assays$RNA@data))
  pseudotime <- paga_traj_output$seurat$nothing@misc$Order
  names(pseudotime) <- colnames(paga_traj_output$seurat$nothing)
  GroundTruth <- dataset %>% add_linear_trajectory(pseudotime) %>% add_dimred(dyndimred::dimred_mds,
                                                                              expression_source = dataset$expression) %>% add_cell_waypoints()
  model <- paga_traj_output$trajectory
  names(model) <- methods
  model <- map(model, add_cell_waypoints)
  # take same metrics as in paper https://www-nature-com.proxy.insermbiblio.inist.fr/articles/s41587-019-0071-9#Sec9
  metric_ids <- dyneval::metrics %>% filter(category != "average") %>% filter(perfect==1) %>% filter(metric_id %in% c(#"him",
    "correlation","F1_branches","featureimp_wcor")) %>% pull(metric_id)
  metrics <- lapply(model,
                    function(x) calculate_metrics2(GroundTruth,x,metric_ids))
  metrics_df <- data.frame("metric_value"=c(sapply(metrics,function(x) return(t(x[,metric_ids])))),
                           "metric_id"=metric_ids,
                           "method_id"=rep(methods,each=length(metric_ids)))
  metrics_df <- bind_rows(metrics_df, data.frame("metric_value"=colMeans(sapply(metrics,function(x) return(t(x[,metric_ids])))),
                                                 "metric_id"="Mean",
                                                 "method_id"=methods))
  metrics_df$metric_id <- factor(metrics_df$metric_id, levels = c(metric_ids,"Mean"))
  #pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/data/forTI/",
  #                paga_traj_output$dataset_id, "/",
  #                paga_traj_output$n_dim,"_dims_paga_heatmap_comp.pdf" ))
  ggplot(metrics_df, aes(metric_id, method_id, fill = metric_value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme(axis.text.x = element_text(angle = 90))
  #dev.off()
}
score_paga_traj <- function(paga_traj_output) {
  dataset <- wrap_expression(counts = t(paga_traj_output$seurat$nothing@assays$RNA@counts),
                             expression = t(paga_traj_output$seurat$nothing@assays$RNA@data))
  pseudotime <- paga_traj_output$seurat$nothing@misc$Order
  if (length(grep("_",pseudotime))>0) {
    pseudotime <- sapply(pseudotime, function(x) strsplit(x,"_")[[1]][2])
  }
  if (length(grep("M",pseudotime))>0|length(grep("W",pseudotime))>0|length(grep("F",pseudotime))>0) {
    pseudotime <- strsplit(pseudotime,"")
    pseudotime <- sapply(pseudotime,
                         function(x) {tmp=x[!is.na(unlist(lapply(x,as.numeric)))];
                         as.numeric(paste(tmp,collapse=""))})
  }
  if (sum(is.na(pseudotime))>0) {
    pseudotime <- paga_traj_output$seurat$nothing@misc$Order
    IDs <- unique(pseudotime)
    pseudotime <- rapply(as.list(seq(length(IDs))),
                         function(x) gsub(IDs[x],x,pseudotime))
    pseudotime <- as.numeric(pseudotime[!is.na(as.numeric(pseudotime))])
  }
  names(pseudotime) <- colnames(paga_traj_output$seurat$nothing)
  GroundTruth <- dataset %>% add_linear_trajectory(pseudotime) %>% add_dimred(dyndimred::dimred_mds,
                                                                              expression_source = dataset$expression) %>% add_cell_waypoints()
  model <- paga_traj_output$trajectory
  names(model) <- methods
  model <- lapply(model,
                  function(x) add_cell_waypoints(x))
  # take same metrics as in paper https://www-nature-com.proxy.insermbiblio.inist.fr/articles/s41587-019-0071-9#Sec9
  metric_ids <- dyneval::metrics %>% filter(category != "average") %>% filter(perfect==1) %>% filter(metric_id %in% metric_choice) %>% pull(metric_id)
  metrics <- lapply(model,
                    function(x) calculate_metrics2(GroundTruth,x,metric_ids))
  metrics <- lapply(metrics,
                    function(x) x[,metric_ids])
  names(metrics) <- methods
  return(metrics)
}
compare_byScore <- function(scores) {
  # first normalize scores across methods for each dataset + apply unit proba density function
  score_norm <- lapply(scores,
                       function(x) {tmp<-normalize(x);
                       tmp<-apply(tmp,2,dnorm);
                       rownames(tmp)<-methods;
                       return(tmp)})
  # arithmetic mean over datasets
  avg_score <- sapply(methods,
                      function(x) {tmp<-sapply(score_norm,
                                               function(y) y[x,]);
                      rowMeans(tmp)})
  data.frame(t(avg_score))
}
compare_allScore <- function(scores) {
  avg_score <- compare_byScore(scores)
  # geometric mean
  data.frame("OverallScore"=apply(avg_score,
                                  1,
                                  function(x) exp(mean(log(x[x!=0])))))
}
make_traj_and_score <- function(dataset_id, n_dim) {
  if (file.exists(paste0(path,dataset_id, "_",n_dim,"dims_",methods[1],"2.h5seurat"))) {
    paga_traj_output <- make_paga_traj(dataset_id, n_dim)
    score <- score_paga_traj(paga_traj_output)
    score <- data.frame(do.call(rbind,score))
    rownames(score)<-methods
    return(score)
  }
}
cluster_nb <- function(dataset_id, n_dim) {
  if (file.exists(paste0(path,dataset_id,"_",n_dim,"dims_",methods[1],"2.h5seurat"))) {
    seurat <- lapply(methods,
                     function(x) LoadH5Seurat(paste0(path,dataset_id,"_",n_dim,"dims_", x,"2.h5seurat"),
                                              verbose=F))
    return("cluster_nb"=sapply(seurat,function(x) length(unique(x$leiden))))
  }
  else {
  }
}


#paga_traj_output <- pblapply(n_dim_choice,
#                       function(dims) pblapply(dataset_choice[1:15],
#                                            function(set) make_paga_traj(set,dims)))
#paga_traj_output <- lapply(paga_traj_output,
#                           function(x1) {idx = rep(T,length(x1));
#                           idx <- sapply(seq(length(x1)),
#                                    function(x2) idx[x2]=ifelse(is.null(x1[[x2]]),F,T));
#                           x1<-x1[idx]}) # remove empty spots
#scores_result <- pblapply(paga_traj_output,
#                          function(x) pblapply(x, score_paga_traj))
cluster_nb <- pblapply(n_dim_choice,
                       function(dims) pblapply(dataset_choice,
                                               function(set) cluster_nb(set,dims)))
scores_result <- pblapply(n_dim_choice,
                          function(dims) {tmp<-pblapply(dataset_choice,
                                                        function(set) make_traj_and_score(set,dims));
                          names(tmp) <- dataset_choice[1:length(tmp)];
                          return(tmp)})
names(scores_result) <- n_dim_choice
scores_result <- lapply(scores_result,
                        function(x1) {idx = rep(T,length(x1));
                        idx <- sapply(seq(length(x1)),
                                      function(x2) idx[x2]=ifelse(is.null(x1[[x2]]),F,T));
                        x1<-x1[idx]}) # remove empty spots
scores_result <- scores_result[lengths(scores_result)!=0]
scores <- lapply(scores_result,
                 function(x) {tmp<-cbind(compare_byScore(x),compare_allScore(x));
                 data.frame("Score"=unlist(tmp),
                            "Metric"=rep(colnames(tmp), each= nrow(tmp)),
                            "Hub_reduction"=rep(rownames(tmp), times=ncol(tmp)))})
pdf(file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/bubbleplot_score_",choice,"_sup.pdf"),
    width=12,height=8)#sup refers to the nb of clusters chosen or not
ggarrange(plotlist=lapply(scores,
                          function(x) ggplot(x, aes(Metric, Hub_reduction)) +
                            geom_point(aes(size=Score, color = Score)) +
                            geom_text(aes(label=round(Score,3))) +
                            scale_size(range = c(0,20)) +
                            scale_color_gradient2(low="navyblue",
                                                  mid="mediumvioletred",
                                                  high="yellow",
                                                  midpoint=0.3) +
                            theme(axis.text.x = element_text(angle = 90))),
          labels=paste0(n_dim_choice," PCs"),
          common.legend = T)
dev.off()
saveRDS(scores_result,
        file=paste0("/Users/elise/Desktop/GitHub/Hubness_sc/Python_scripts/TI/score_",choice,"_sup.rds"))

