library(magrittr)
library(ggalt)
library(ggfortify)
library(ggtext)
# library(ggVennDiagram)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

plot_var_explained <- function (object, x = "view", y = "factor", split_by = NA, plot_total = FALSE, 
          factors = "all", min_r2 = 0, max_r2 = NULL, legend = TRUE, 
          use_cache = TRUE, ...) 
{

  if (length(unique(c(x, y, split_by))) != 3) {
    stop(paste0("Please ensure x, y, and split_by arguments are different.\n", 
                "  Possible values are `view`, `group`, and `factor`."))
  }
  if (is.na(split_by)) 
    split_by <- setdiff(c("view", "factor", "group"), c(x, 
                                                        y, split_by))
  if ((use_cache) & .hasSlot(object, "cache") && ("variance_explained" %in% 
                                                  names(object@cache))) {
    r2_list <- object@cache$variance_explained
  }
  else {
    r2_list <- calculate_variance_explained(object, factors = factors, 
                                            ...)
  }
  r2_mk <- r2_list$r2_per_factor
  r2_mk_df <- melt(lapply(r2_mk, function(x) melt(as.matrix(x), 
                                                  varnames = c("factor", "view"))), id.vars = c("factor", 
                                                                                                "view", "value"))
  colnames(r2_mk_df)[ncol(r2_mk_df)] <- "group"
  if ((length(factors) == 1) && (factors[1] == "all")) {
    factors <- factors_names(object)
  }
  else {
    if (is.numeric(factors)) {
      factors <- factors_names(object)[factors]
    }
    else {
      stopifnot(all(factors %in% factors_names(object)))
    }
    r2_mk_df <- r2_mk_df[r2_mk_df$factor %in% factors, ]
  }
  r2_mk_df$factor <- factor(r2_mk_df$factor, levels = factors)
  r2_mk_df$group <- factor(r2_mk_df$group, levels = groups_names(object))
  r2_mk_df$view <- factor(r2_mk_df$view, levels = views_names(object))
  groups <- names(r2_list$r2_total)
  views <- colnames(r2_list$r2_per_factor[[1]])
  if (!is.null(min_r2)) 
    r2_mk_df$value[r2_mk_df$value < min_r2] <- 0.001
  min_r2 = 0
  if (!is.null(max_r2)) {
    r2_mk_df$value[r2_mk_df$value > max_r2] <- max_r2
  }
  else {
    max_r2 = max(r2_mk_df$value)
  }
  
  title <- round(head(model@cache$variance_explained$r2_total[[1]]), 1)
  title <- paste0(paste0(c(names(object@data[1]), names(object@data[2])), c(": ", ": "), 
                         as.character(title), c("%", "%")), collapse = ', ')
  print(title)
  
  p1 <- ggplot(r2_mk_df, aes(x = .data[[x]], y = .data[[y]])) + 
    geom_tile(aes(fill = .data$value), color = "black") + 
    geom_text(aes(.data[[x]], .data[[y]], label=round(.data$value, 1)), colour = "red", check_overlap = TRUE) +
    facet_wrap(as.formula(sprintf("~%s", split_by)), nrow = 1) + 
    labs(x = "", y = "", title = "") + scale_fill_gradientn(colors = c("gray97", 
                                                                       "darkblue"), guide = "colorbar", limits = c(min_r2, max_r2)) + 
    guides(fill = guide_colorbar("Var. (%)")) + theme(axis.text.x = element_text(size = rel(1), 
                                                                                 color = "black"), axis.text.y = element_text(size = rel(1.1), 
                                                                                                                              color = "black"), axis.line = element_blank(), axis.ticks = element_blank(), 
                                                      panel.background = element_blank(), strip.background = element_blank(), 
                                                      strip.text = element_text(size = rel(1))) +
    ggtitle(title)
  
  if (isFALSE(legend)) 
    p1 <- p1 + theme(legend.position = "none")
  if (length(unique(r2_mk_df[, split_by])) == 1) 
    p1 <- p1 + theme(strip.text = element_blank())
  if (plot_total) {
    r2_m_df <- melt(lapply(r2_list$r2_total, function(x) lapply(x, 
                                                                function(z) z)), varnames = c("view", "group"), value.name = "R2")
    colnames(r2_m_df)[(ncol(r2_m_df) - 1):ncol(r2_m_df)] <- c("view", 
                                                              "group")
    r2_m_df$group <- factor(r2_m_df$group, levels = MOFA2::groups_names(object))
    r2_m_df$view <- factor(r2_m_df$view, levels = views_names(object))
    min_lim_bplt <- min(0, r2_m_df$R2)
    max_lim_bplt <- max(r2_m_df$R2)
    p2 <- ggplot(r2_m_df, aes(x = .data[[x]], y = .data$R2)) + 
      geom_bar(stat = "identity", fill = "deepskyblue4", 
               color = "black", width = 0.9) + facet_wrap(as.formula(sprintf("~%s", 
                                                                             split_by)), nrow = 1) + xlab("") + ylab("Variance explained (%)") + 
      scale_y_continuous(limits = c(min_lim_bplt, max_lim_bplt), 
                         expand = c(0.005, 0.005)) + theme(axis.ticks.x = element_blank(), 
                                                           axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), 
                                                           axis.title.y = element_text(color = "black"), axis.line = element_line(color = "black"), 
                                                           panel.background = element_blank(), strip.background = element_blank(), 
                                                           strip.text = element_text())
    if (length(unique(r2_m_df[, split_by])) == 1) 
      p2 <- p2 + theme(strip.text = element_blank())
    plot_list <- list(p1, p2)
  }
  else {
    plot_list <- p1
  }
  return(plot_list)
}



make_pathway_pcas <- function(pw_list, bestf_list, idmap, w = 2500, kegg.gs, e, mycolor, myshape, axis1 = 1, axis2 = 2){ #f
  for (i in 1:length(pw_list)){
    for (j in 1:3){
      pw <- kegg.gs[grepl(pw_list[[i]], names(kegg.gs))]	
      name <- paste0(names(bestf_list[j]), " - ", names(pw_list[i]), " pathway")
      print(name)
      mydf <- bestf_list[[j]]  %>%
        filter(type == "pep") %>% 
        left_join(idmap, by = c("ID" = "anonID")) %>%
        filter(!is.na(entrezid), entrezid %in% pw[[1]]) %>%
        select(matches("inspire")) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Patientcode")
      
      p <- plot_pca(df_topca = mydf, e = e, mycolor, myshape, mytitle = name, axis1 = axis1, axis2 = axis2)
      
      ggsave(paste0("plots/pca/pca_per_pathway/", name, ".png"), p, width = 11, height = 8, units = "cm")
      

    }
  }
  
}

filter_data_for_pw <- function(fc_df, idmap, pw){
  fc_df$entrezid <- as.character(fc_df$entrezid)
  idmap$entrezid <- as.character(idmap$entrezid)
  fc_df %>%
    filter(entrezid %in% pw[[1]]) %>% 
    left_join(idmap, by = "entrezid")  %>% select(-entrezid)

}


parse_dfs_to_heatmapformat <- function(prot_df, pep_df, idmap, pw){
  prot_df_pw <- prot_df %>%
    filter(!is.na(ID), ID %in% pw[[1]]) %>% 
    left_join(idmap, by = c("ID" = "entrezid")) %>% mutate(type = "Pr")
  pep_df_pw <- pep_df %>%
    filter(!is.na(ID), ID %in% pw[[1]]) %>% 
    left_join(idmap, by = c("ID" = "entrezid")) %>% mutate(type = "")
  
  alldf <- rbind(prot_df_pw, pep_df_pw) %>% select(-ID)
  return(alldf)
  
}


# make_heatmaps <- function(heatmap_mat, name, ha, idmap, w = 2500, e, rows_clustered = T, 
#                           cols_clustered = F, coloring = "normal", row_annotation = F, res.path){ 
#   
#   
#   
#   plot_heatmap(df = heatmap_mat, ha = ha, name = name, e = e, w = w, coloring = coloring, 
#                rows_clustered = rows_clustered, cols_clustered = cols_clustered, row_annotation = row_annotation,
#                res.path = res.path)
# }

make_pathway_tilemaps <- function(pw_list, fc_df, idmap, manual_size = F, w = NULL, h = NULL, kegg.gs, rows_clustered = T, 
                          cols_clustered = F, coloring = "normal", res.path = "results/pathway_tilemaps/"){ 
  # Makes an FC tilemap for each pathway.
  # Input:
  # coloring: "normal" or "limited_scale"
  # fc_df: columns should be: comparison (values e.g. "T vs C"; arranged factor), entrezid, data_from (values e.g. "Pep", 'Ki', 'MS'; arranged factor), and value (LFC).
  # idmap: entrezid to uniprot name, contains all data
  
  
  
  if (!dir.exists(res.path)){
    dir.create(res.path, recursive = T)
  }
  
  # for each pathway:
  for (i in 1:length(pw_list)){
    # get pw elements

    pw <- kegg.gs[grepl(pw_list[[i]], names(kegg.gs))]	
    name <- paste0(names(pw_list[i]))
    print(name)
    #browser()
    pws_not_found <- c()
    if (length(pw) == 0) {
      print(paste("Pathway not found: ", name))
      pws_not_found <- c(pws_not_found, name)
      next  # Skip to the next iteration
    } else {
      # filter data for pw elements
      mymat <- filter_data_for_pw(fc_df = fc_df, idmap = idmap, pw = pw)
      plot_pathway_tilemap(df = mymat, name = name, manual_size = F,  w = NULL, h = NULL, coloring = coloring, 
                           rows_clustered = rows_clustered, cols_clustered = cols_clustered,
                           res.path = res.path)
    }
    print("Pathways not found")
    print(pws_not_found)
  }
}



parse_bestf_to_pcaformat <- function(bestf, idmap, pw){ 
  bestf %>%
    filter(type == "pep") %>% 
    left_join(idmap, by = c("ID" = "anonID")) %>%
    filter(!is.na(entrezid), entrezid %in% pw[[1]]) %>%
    select(-ID) %>%
    rename("ID" = "ID.y") %>%
    select(-FactorNum, -type, -coef, -PepProtein_UniprotID, -PepProtein_SeqSimilarity, -Sequence, -entrezid, -PepProtein_UniprotName) %>%
    distinct() %>%
    arrange(ID) %>%
    column_to_rownames("ID")%>%
    as.matrix()
}



plot_heatmap <- function(df, ha = NULL, value_limits = NULL, name,legend_name, w = 3000, h = 2000, rows_clustered = F, cols_clustered = F, 
                         show_row_names = F, coloring = "normal", res.path){
  # df: matrix with optional row and colnames
  # ha: heatmap annotation df with first column containing the colnames of X, 2nd col is annotation
  # value_limits: if not null, clips the df according to min and max value. example: c(0,1)
  
  parse_ha <- function(ha, X){
    # parses the ha input to a vector with names of the same colnames as X, values as the annotation. 
    # Orders the vector in the same way as X
    first_col <- colnames(ha)[1]
    second_col <- colnames(ha)[2]
    col_names <- colnames(X)
    ha_vec <- setNames(ha[[second_col]][match(col_names, ha[[first_col]])], col_names)
    return(ha_vec)
  }
  
  
  
  if(!is.null(value_limits)){
    df <- as.data.frame(df) %>%
      rownames_to_column("row_id") %>%
      pivot_longer(-row_id, names_to = "col_id", values_to = "value") %>%
      mutate(value = pmax(value_limits[1], pmin(value, value_limits[2]))) %>%
      pivot_wider(names_from = col_id, values_from = value) %>%
      column_to_rownames("row_id") %>% as.matrix()
  }
  
  if (!dir.exists(res.path)){
    dir.create(res.path, recursive = T)
  }

  # coloring
  if (coloring == "diff_normal"){
    colorfunc <- circlize::colorRamp2(seq(min(df), max(df), length = 3), c("green", "black", "red"))
  } else if (coloring == 'diff_limited'){
    colorfunc <- circlize::colorRamp2(seq(-1, 1, length = 3), c("green", "black", "red"))
  } else if (coloring == "viridis"){
    colorfunc <- viridis::viridis(100)
  } else if (coloring == "onlyneg"){
    colorfunc <- circlize::colorRamp2(seq(-7, 0, length = 2), c("green", "black"))
  } else if (coloring == "onlypos"){
    colorfunc <- circlize::colorRamp2(seq(0, (max(df)), length = 2), c("blue", "red"))
  } else if (coloring == 'discrete'){
    colorfunc <- c("DOWN" = "green", "0" = "black", "UP" = "red") 
  }

  # Kinase group annotation and colors
  ha <- NULL
  if (!is.null(ha)) {
    ha_vec <- parse_ha(ha, X)
    kinase_groups_unique <- unique(ha_vec)
    kinase_group_colors <- setNames(RColorBrewer::brewer.pal(min(length(kinase_groups_unique), 12), "Paired")[seq_along(kinase_groups_unique)], kinase_groups_unique)
    ha <- ComplexHeatmap::HeatmapAnnotation(Kinase_Group = ha_vec, 
                                            col = list(Kinase_Group = kinase_group_colors))
  }

  png(filename = paste0(res.path, "/Heatmap_", name, ".png"), 
      width = w, height = h, units = "px", res = 150)

  heatmap_args <- list(
    df,
    name = legend_name, 
    cluster_columns = cols_clustered,
    cluster_rows = rows_clustered,
    show_row_names = show_row_names,
    show_column_names = TRUE,
    column_names_side = "top",
    row_km = 0,
    column_km = 0,
    column_title_rot = 0, 
    column_title = name,
    col = colorfunc,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8)
  )

  if (!is.null(ha)) {
    heatmap_args$top_annotation <- ha
    p <- do.call(Heatmap, heatmap_args)
  } else {
    p <- do.call(Heatmap, heatmap_args)
  }

  p <- draw(p)
  dev.off()
}


plot_pathway_tilemap <- function(df, name, manual_size, w, h, rows_clustered = F, 
                                 cols_clustered = F, coloring = "normal",
                                 res.path){
  

  if (coloring == "limited_scale"){
    p25_neg <- stats::quantile(df$value[df$value < 0], probs = 0.25, na.rm = TRUE)
    p75_pos <- stats::quantile(df$value[df$value > 0], probs = 0.75, na.rm = TRUE)
    
    limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
    color_limits <- c(-limit_val, limit_val)
    
    df$value <- pmax(pmin(df$value, color_limits[2]), color_limits[1])
  }  
  
  

  p <- ggplot(df, aes(x = data_from, y = uniprotname, fill = value)) + 
    geom_tile() +
    facet_wrap(df$comparison, nrow = 1) + 
    labs(x = "Comparison", y = "", title = name) + 
    theme(axis.text.x = element_text(size = 15, angle = 35, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 20),
          strip.text = element_text(size = 15),
          title = element_text(size = 20)) +
    scale_fill_gradient2(low = "green",
                         mid = "black",
                         high = "red")
  
  if (manual_size == F){
    # set width
    n_data_from <- length(unique(df$data_from))
    
    if(n_data_from == 1){
      w = 30
    } else if (n_data_from == 2){
      w = 40
    } else if (n_data_from == 3){
      w = 50
    } else{
      w = 60
    }
    
    
    if (nrow(df) < 400){
      h = 60
    } else if (nrow(df) >= 400  & nrow(df) <= 500){
      h = 70
    } else if (nrow(df) >500 & nrow(df <= 750)){
      h = 80
    } else if (nrow(df) > 750){
      h = 90
    }
  }
  
  ggsave(paste0(res.path, "/", name, ".png"), p, width = w, height = h, units = "cm" )
  
}


plot_rna_pep_hist <- function(rnav, pepv, titles, data_from){
  ncol_fig <- length(rnav)
  
  jpeg(paste0("results/plots/hist_prot_pep_normalizations_", data_from, ".jpg"), width = 1200, height = 600, units = "px")
  
  par(mfrow= c(2,ncol_fig))
  #par(mar = c(2.5, 2.5, 2.5, 2.5))
  for (i in 1:ncol_fig){
    hist(rnav[[i]], 20, main = titles[i], xlab = 'Protein')
  } 
  for (j in 1:ncol_fig){
    hist(pepv[[j]], 20, main = NULL, xlab = "pep")
  }
  
  dev.off()
}




############
### MOFA ###
############





plot_factor_pca <- function(model, mycolor, myshape, f1, f2, perfactor_var, w = 14, h = 14){
  
  a <- plot_factors_pg(model, 
                    factors = c(f1, f2),
                    color_by = mycolor,
                    shape_by = myshape,
                    #color_name = "sample_type",
                    dot_size = 2.5,
                    scale = F,
                    stroke = NA) + 
    xlab(paste0("Factor ", f1, " (R2 = ", round(sum(perfactor_var[paste0('Factor', f1),]), 1), "%)")) + 
    ylab(paste0("Factor ", f2, " (R2 = ", round(sum(perfactor_var[paste0('Factor', f2),]), 1), "%)")) +
    theme_gray() +
    theme(axis.text.x= element_text(colour = 'black'),
          axis.text.y = element_text(colour = 'black')) +
    ggtitle(paste0(modelname))
  
  
  
  myfilename = paste0("./results/plots/mofamodel_", modelname, 
                      "/factorplot_color_", mycolor, 
                      "_shape_", myshape,
                      "_f", f1, "-f", f2,
                      ".png")
  ggsave(myfilename, a, width = w, height = h, units = "cm")
  a
}

plot_factors_pg <- function (object, factors = c(1, 2), groups = "all", show_missing = TRUE, 
          scale = FALSE, color_by = NULL, shape_by = NULL, color_name = NULL, 
          shape_name = NULL, dot_size = 2, alpha = 1, legend = TRUE, 
          stroke = NULL, return_data = FALSE) {
  
  
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (length(unique(factors)) == 1) {
    .args <- as.list(match.call()[-1])
    .args <- .args[names(.args) != "factors"]
    return(do.call(plot_factor, c(.args, list(factors = unique(factors)))))
  }
  else if (length(factors) > 2) {
    .args <- as.list(match.call()[-1])
    p <- do.call(.plot_multiple_factors, .args)
    return(p)
  }
  
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name)) 
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name)) 
    shape_name <- shape_by
  factors <- .check_and_get_factors(object, factors)
  Z <- get_factors(object, factors = factors, groups = groups, 
                   as.data.frame = TRUE)

  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  Z <- Z[complete.cases(Z), ]
  df <- merge(Z, color_by, by = "sample")
  df <- merge(df, shape_by, by = "sample")
  df$shape_by <- as.character(df$shape_by)
  if (isFALSE(show_missing)) 
    df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  df <- spread(df, key = "factor", value = "value")
  df <- df[, c(colnames(df)[seq_len(4)], factors)]
  df <- set_colnames(df, c(colnames(df)[seq_len(4)], "x", "y"))
  
  if (scale) {
    df$x <- df$x/max(abs(df$x))
    df$y <- df$y/max(abs(df$y))
  }
  if (return_data) 
    return(df)
  if (is.null(stroke)) {
    stroke <- .select_stroke(N = length(unique(df$sample)))
  }

  p <- ggplot(df, aes(x = .data$x, y = .data$y, color = .data$color_by, 
                      shape = .data$shape_by)) + 
    geom_point(size = dot_size, alpha = alpha) +
    #geom_encircle(aes(group = .data$shape_by, fill = .data$shape_by), alpha = 0.1, expand = 0, s_shape = 1) +
    labs(x = factors[1], y = factors[2]) + 
    #scale_shape_manual(values=c(16, 17)) + 
    scale_shape_manual(values=1:nlevels(df %>% pull(shape_by) %>% as.factor())) +
    theme(axis.text = element_text(size = rel(0.8), color = "black"), 
          axis.title = element_text(size = rel(1.1), color = "black"), 
          axis.line = element_line(color = "black", linewidth = 0.5), 
          axis.ticks = element_line(color = "black", linewidth = 0.5))
  
  #browser()
  #p <- .add_legend(p, df, legend, color_name, shape_name)
  if (!is.null(color_name)) {
    p <- p + labs(color = color_name)
  }
  if (!is.null(shape_name)) {
    p <- p + labs(shape = shape_name)
  }
  return(p)
}

.plot_multiple_factors <- function(object, factors = "all", show_missing = TRUE, dot_size = 1,
                                   color_by = NULL, color_name = "", shape_by = NULL, shape_name = "") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Collect relevant data
  Z <- get_factors(object, factors=factors, as.data.frame=TRUE)
  
  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Remove samples with missing factor values
  Z <- Z[complete.cases(Z),]
  
  # Merge factor values with color and shape information
  df <- merge(Z, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  
  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # Spread over factors
  df <- tidyr::spread(df, key="factor", value="value")
  
  # Prepare the legend
  p <- ggplot(df, aes_string(x=factors[1], y=factors[2], color="color_by", shape="shape_by")) +
    geom_point() +
    scale_shape_manual(values=1:nlevels(df %>% pull(shape_by) %>% as.factor()))
    theme(
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_text(size=rel(1.2))
    )
  if (length(unique(df$color))>1) { p <- p + labs(color=color_name) } else { p <- p + guides(color=FALSE) + scale_color_manual(values="black") }
  if (is.numeric(df$color)) p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n=5, name="RdYlBu")))(10)) 
  if (length(unique(df$shape))>1) { p <- p + labs(shape=shape_name) } else { p <- p + guides(shape = FALSE) }
  if (length(unique(df$color))>1 || length(unique(df$shape))>1) { legend <- GGally::grab_legend(p) } else { legend <- NULL }
  
  
  # Generate plot
  p <- GGally::ggpairs(df, 
                       columns = factors,
                       lower = list(continuous=GGally::wrap("points", size=dot_size)), 
                       diag = list(continuous='densityDiag'), 
                       upper = list(continuous=GGally::wrap("points", size=dot_size)), 
                       mapping = aes_string(color="color_by", shape="shape_by"), 
                       title = "", 
                       legend = legend
  )
  p <- p + theme_bw() + theme(
    # axis.line = element_line(color="black", size=rel(1.0)),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
  
  return(p)
}

.set_shapeby <- function (object, shape_by) 
{
  if (is.null(shape_by)) {
    shape_by <- rep("1", sum(object@dimensions[["N"]]))
  }
  else if (shape_by[1] == "group") {
    shape_by <- samples_metadata(object)$group
  }
  else if ((length(shape_by) == 1) && is.character(shape_by) & 
           (shape_by %in% colnames(samples_metadata(object)))) {
    shape_by <- samples_metadata(object)[, shape_by]
  }
  else if ((length(shape_by) == 1) && is.character(shape_by) && 
           (shape_by[1] %in% unlist(features_names(object)))) {
    viewidx <- which(sapply(features_names(object), function(x) shape_by %in% 
                              x))
    foo <- list(shape_by)
    names(foo) <- names(viewidx)
    shape_by <- lapply(get_data(object, features = foo), 
                       function(l) Reduce(cbind, l))[[1]][1, ]
  }
  else if (is(shape_by, "data.frame")) {
    stopifnot(all(colnames(shape_by) %in% c("sample", "color")))
    stopifnot(all(unique(shape_by$sample) %in% unlist(samples_names(object))))
  }
  else if (length(shape_by) > 1) {
    stopifnot(length(shape_by) == sum(object@dimensions[["N"]]))
  }
  else {
    stop("'shape_by' was specified but it was not recognised, please read the documentation")
  }
  if (!is(shape_by, "data.frame")) {
    df = data.frame(sample = samples_metadata(object)$sample, 
                    shape_by = as.factor(shape_by), stringsAsFactors = FALSE)
  }
  return(df)
}

.set_colorby <- function (object, color_by) 
{
  if (is.null(color_by)) {
    color_by <- rep("1", sum(object@dimensions[["N"]]))
  }
  else if (color_by[1] == "group") {
    color_by <- samples_metadata(object)$group
  }
  else if ((length(color_by) == 1) && (is.character(color_by) | 
                                       is.factor(color_by)) & (color_by[1] %in% colnames(samples_metadata(object)))) {
    color_by <- samples_metadata(object)[, color_by]
  }
  else if ((length(color_by) == 1) && is.character(color_by) && 
           (color_by[1] %in% unlist(features_names(object)))) {
    viewidx <- which(sapply(features_names(object), function(x) color_by %in% 
                              x))
    foo <- list(color_by)
    names(foo) <- names(viewidx)
    color_by <- lapply(get_data(object, features = foo), 
                       function(l) Reduce(cbind, l))[[1]][1, ]
  }
  else if ((length(color_by) == 1) && is.character(color_by) && 
           (color_by[1] %in% colnames(get_factors(object)[[1]]))) {
    color_by <- do.call(rbind, get_factors(object))[, color_by]
  }
  else if (is(color_by, "data.frame")) {
    stopifnot(all(colnames(color_by) %in% c("sample", "color")))
    stopifnot(all(unique(color_by$sample) %in% unlist(samples_names(object))))
  }
  else if (length(color_by) > 1) {
    stopifnot(length(color_by) == sum(get_dimensions(object)$N))
  }
  else {
    stop("'color_by' was specified but it was not recognised, please read the documentation")
  }
  if (!is(color_by, "data.frame")) {
    df <- data.frame(sample = samples_metadata(object)$sample, 
                     color_by = color_by, stringsAsFactors = FALSE)
  }
  if (length(unique(df$color_by)) < 5) 
    df$color_by <- as.factor(df$color_by)
  return(df)
}

.check_and_get_factors <- function (object, factors) 
{
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(factors)))
  if (is.numeric(factors)) {
    stopifnot(all(factors <= object@dimensions$K))
    factors_names(object)[factors]
  }
  else {
    if (paste0(factors, collapse = "") == "all") {
      factors_names(object)
    }
    else {
      stopifnot(all(factors %in% factors_names(object)))
      factors
    }
  }
}


.select_stroke <- function (N) 
{
  if (N <= 1000) {
    stroke <- 0.5
  }
  else if (N > 1000 & N <= 10000) {
    stroke <- 0.2
  }
  else {
    stroke <- 0.05
  }
}


.add_legend <- function (p, df, legend, color_name, shape_name) 
{
  if (is.numeric(df$color_by)) {
    p <- p + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 5, 
                                                                           name = "RdYlBu")))(10)) + labs(fill = color_name)
  }
  else {
    if (length(unique(df$color_by)) > 1) {
      p <- p + guides(fill = guide_legend(override.aes = list(shape = 21))) + 
        labs(fill = color_name)
    }
    else {
      p <- p + guides(fill = "none", color = "none") + 
        scale_color_manual(values = "black") + scale_fill_manual(values = "gray60")
    }
  }
  if (length(unique(df$shape_by)) > 1) {
    p <- p + scale_shape_manual(values = c(21, 23, 24, 25)[1:length(unique(df$shape_by))]) + 
      guides(shape = guide_legend(override.aes = list(fill = "black"))) + 
      labs(shape = shape_name)
  }
  else {
    p <- p + scale_shape_manual(values = c(21)) + guides(shape = "none")
  }
  if (legend) {
    p <- p + guides(color = guide_legend(override.aes = list(fill = "white"))) + 
      theme(legend.text = element_text(size = rel(0.8)), 
            legend.title = element_text(size = rel(0.8)), 
            legend.key = element_rect(fill = "white", color = "white"))
  }
  else {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

###########
### PCA ###
###########

plot_pca_continous_color <- function(df_topca, e, joinby, mycolor, myshape, 
                     mytitle, axis1, axis2){
  
  pca.res <- prcomp(df_topca[,-1], center = T)
  df_e <- left_join(df_topca, e, by = joinby) %>% select(all_of(c(mycolor, myshape)))
  
  pca_data <- as.data.frame(pca.res$x)
  pca_data_e <- cbind(pca_data, df_e)
  #browser()
  scoreplot <- ggplot(pca_data_e,
                      aes(x = .data[[axis1]], y = .data[[axis2]], color = .data[[mycolor]], shape = .data[[myshape]])) + 
    geom_point(size = 2.5) + 
    scale_color_gradient(low = "#32CD32", high = "red") +
    ggtitle(mytitle) +
    theme_classic() +
    theme(axis.text.x= element_text(colour = 'black'),
          axis.text.y = element_text(colour = 'black'),
          plot.title = element_textbox_simple())
  
  if (!is.null(myshape)){
    scoreplot <- scoreplot + scale_shape_manual(values=1:nlevels(pca_data_e %>% pull(myshape) %>% as.factor()))
  }

  scoreplot
  
}

plot_pca2 <- function(df_topca, e, joinby, mycolor, myshape = NULL,
                      mytitle, toscale = FALSE, sel_pcs = 1:5, mypalette) {

  # Set rownames
  df_topca <- df_topca %>% column_to_rownames(joinby)
  
  # PCA
  pca.res <- prcomp(df_topca, center = TRUE, scale. = toscale)
  
  # Get variation for the figure
  pca_var <- pca.res$sdev^2
  pca_var_exp <- round(100 * pca_var / sum(pca_var), 1)
  
  new_names <- paste0("PC", 1:length(pca_var_exp), " (", pca_var_exp, "%)")
  colnames(pca.res$x) <- new_names
  pc_labels <- new_names[sel_pcs] 
  
  # Join scores with metadata
  pcaScores <- pca.res$x %>% 
    as.data.frame() %>% 
    mutate(join_id = rownames(.)) %>% 
    left_join(e, by = c('join_id' = joinby))
  #browser()
  # make shape/color columns factors if needed
  if (!is.null(myshape)) {
    pcaScores[[myshape]] <- as.factor(pcaScores[[myshape]])
  }

  pcaScores[[mycolor]] <- as.factor(pcaScores[[mycolor]])
  
  # Build aesthetics string for both color and shape
  mapping <- if (!is.null(myshape)) {
    ggplot2::aes_string(color = mycolor, shape = myshape)
  } else {
    ggplot2::aes_string(color = mycolor)
  }

  # Create plot
  scoreplot <- GGally::ggpairs(
    data = pcaScores,
    columns = pc_labels,
    mapping = mapping,
    upper = list(continuous = "blank")
  ) + 
    scale_color_brewer(palette = mypalette) + 
    ggtitle(mytitle) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # Dummy plot for legend
  if (!is.null(myshape)) {
    dummy <- ggplot(pcaScores, aes(
      x = !!sym(pc_labels[1]),
      y = !!sym(pc_labels[2]),
      color = !!sym(mycolor),
      shape = !!sym(myshape)
    )) +
      geom_point() +
      scale_color_brewer(palette = mypalette) + # ensure legend matches plot
      theme_void() +
      theme(legend.position = "right")
  } else {
    dummy <- ggplot(pcaScores, aes(
      x = !!sym(pc_labels[1]),
      y = !!sym(pc_labels[2]),
      color = !!sym(mycolor)
    )) +
      geom_point() +
      scale_color_brewer(palette = mypalette) + # ensure legend matches plot
      theme_void() +
      theme(legend.position = "right")
  }

  legend <- cowplot::get_legend(dummy)
  plot_grob <- GGally::ggmatrix_gtable(scoreplot)

  # Add a new column to the right for the legend
  combined <- gtable::gtable_add_cols(plot_grob, widths = unit(1.5, "in"))

  # Add legend to all rows, centered vertically
  combined <- gtable::gtable_add_grob(
    combined,
    grobs = legend,
    t = 1, l = ncol(combined), b = nrow(combined)
  )

  return(combined)
  
}



plot_pcas2 <- function(df_list, colors, shapes, e, joinby, 
                       namesfrom, valuesfrom,
                       savetodir = "results/plots/PCA",
                       w = 12, h = 10, discretecolor = T,
                       oneplot = F, toscale = F, sel_pcs = 1:5,
                       axis1 = 1, axis2 = 2, mypalette = "Paired"){
  # mycolor should be a vector of colors
  # oneplot means to output only one pca plot with specified axes

  p_list <- as.list(rep(1, length(df_list)))


  lapply(1:length(colors), function(i){
    colorx <- colors[i]
    shapex <- shapes[i]
    # If shapex is NA or NULL, set to NULL for downstream logic
    if (is.na(shapex) || is.null(shapex)) shapex <- NULL
    for (i in 1:length(df_list)){
      topca <- df_list[[i]] %>% 
        pivot_wider(names_from = all_of(namesfrom), values_from = all_of(valuesfrom), values_fill = 0)
      print(paste0("Working on: ", names(df_list[i]), ", color: ", colorx, ", shape: ", shapex))
      
      if(discretecolor){
        if (oneplot){
          combs <- combn(sel_pcs, 2, simplify = FALSE)
          # browser()
          lapply(combs, function(x){
            # make dir
            axisdir <- paste0(savetodir, "/ax", x[1], "-", x[2])
            if (!dir.exists(axisdir)){
              dir.create(axisdir, recursive = T)
            }
            myfilename = paste0(savetodir, "/ax", x[1], "-", x[2], "/", names(df_list[i]), "_PCA_color-", colorx,
                                "-shape-", shapex, ".png")
            myfilename_bip = paste0(savetodir, "/ax", x[1], "-", x[2], "/", names(df_list[i]), "_PCA_color-", colorx,
                                "-shape-", shapex, "_biplot.png")
            
            ps <- plot_pca_individual(df_topca = topca, e = e, joinby = joinby, mycolor = colorx, myshape = shapex, 
                                     mytitle = names(df_list[i]), axis1 = x[1], axis2 = x[2], mypalette = mypalette)
            ggsave(myfilename, ps$s, width = w, height = h, units = "cm")
            ggsave(myfilename_bip, ps$b, width = w, height = h, units = "cm")
            p_list[[i]] <- p
          })
        } else {
          # multiple plots with GGally
          myfilename = paste0(savetodir, "/", names(df_list[i]), "_PCA_color-", colorx, "_pcs", paste(sel_pcs, collapse = "-"), 
                              "-shape-", shapex, "_scale", toscale, ".png")
          mytitle = paste0( names(df_list[i]), " - ", colorx, " - ", shapex)
          p <- plot_pca2(df_topca = topca, e = e, joinby = joinby, mycolor = colorx, myshape = shapex,
                         mytitle = mytitle, toscale = toscale, sel_pcs = sel_pcs, mypalette = mypalette)
          ggsave(myfilename, p, width = w, height = h, units = "cm")
          p_list[[i]] <- p
        }
      } else{
        # Continous color pca
        myfilename = paste0(savetodir, "/ax", axis1, "-", axis2, "/", names(df_list[i]), "_PCA_color-", colorx,
                            "-shape-", shapex, ".png")
        p <- plot_pca_continous_color(df_topca = topca, e = e, joinby = joinby, mycolor = colorx, myshape = shapex, 
                                      mytitle = names(df_list[i]), axis1 = axis1, axis2 = axis2)
        ggsave(myfilename, p, width = w, height = h, units = "cm")
        p_list[[i]] <- p
      }
      
      
      
      
      
      names(p_list) <- paste0(names(df_list), "_p")
      p_list
    }
  })
}



plot_pca_individual <- function(df_topca, e, joinby, mycolor, myshape,
                     mytitle, axis1, axis2, frame = F, mypalette){
  #browser()
  pca.res <- prcomp(df_topca %>% column_to_rownames(joinby), center = T)
  biplot <- factoextra::fviz_pca_biplot(pca.res, repel = TRUE,
                                        col.var = "#2E9FDF", col.ind = "#696969",
                                        select.var = list(contrib = 20),
                                        axes = c(axis1, axis2)) +
    theme_bw()
  
  
  df_e <- left_join(df_topca, e, by = joinby)
  
  scoreplot <- autoplot(pca.res, data = df_e,
                        color = mycolor, size = 2.5, shape = myshape,
                        x = axis1, y = axis2, scale = 0,
                        frame = frame, frame.colour = myshape, frame.alpha = 0.05) +
    #scale_x_continuous(limits = c(-0.3, 1)) +
    #xlab("Dimension 1 (variance = 12.4%)") +
    scale_color_brewer(palette = mypalette) + 
    ggtitle(mytitle) +
    theme_classic() +
    theme(axis.text.x= element_text(colour = 'black'),
          axis.text.y= element_text(colour = 'black'),
          plot.title = element_textbox_simple())

  # Only add shape scale if myshape is not NULL and is a valid column in e
  if (!is.null(myshape) && myshape %in% colnames(e)) {
    scoreplot <- scoreplot + scale_shape_manual(values=1:nlevels(df_e %>% pull(myshape) %>% as.factor()))
  }

  return(list(s = scoreplot, b = biplot))

}

plot_pcas <- function(df_list, mycolor, myshape, e, joinby,
                      namesfrom, valuesfrom,
                      axis1 = 1, axis2 = 2, savetodir = "results/plots/PCA",
                      w = 12, h = 10, discretecolor = T){
  p_list <- as.list(rep(1, length(df_list)))

  if (!dir.exists(savetodir)){
    dir.create(savetodir, recursive = T)
  }

  for (i in 1:length(df_list)){
    myfilename = paste0(savetodir, "/", names(df_list[i]), "_PCA_color-", mycolor,
                        "-shape-", myshape, "-ax-", axis1, "-", axis2, ".png")
    topca <- df_list[[i]] %>%
      pivot_wider(names_from = all_of(namesfrom), values_from = all_of(valuesfrom), values_fill = 0)
    print(paste0("Working on: ", names(df_list[i])))
    if(discretecolor){
      p <- plot_pca(df_topca = topca, e = e, joinby = joinby, mycolor = mycolor, myshape = myshape,
                    mytitle = names(df_list[i]), axis1 = axis1, axis2 = axis2)
    } else{
      p <- plot_pca_continous_color(df_topca = topca, e = e, joinby = joinby, mycolor = mycolor, myshape = myshape,
                    mytitle = names(df_list[i]), axis1 = axis1, axis2 = axis2)
    }

    p_list[[i]] <- p
    ggsave(myfilename, p, width = w, height = h, units = "cm")
  }
  names(p_list) <- paste0(names(df_list), "_p")
  p_list
}



# 
# do_venn <- function(list_of_sets, savepath, mytitle, w = 11, h = 8){
#   p <- ggVennDiagram(list_of_sets,
#                      label = "count",
#                      label_alpha = 0 ) +
#     scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
#     scale_x_continuous(expand = expansion(mult = .5)) +
#     ggtitle(mytitle)
# 
#   filename = paste0(savepath, "/Venn ", paste0(names(list_of_sets), collapse = "-"), "-", mytitle, '.jpg')
#   
#   ggsave(filename, p, width = 11, height = 8, units = "cm")
# }


############################
### Kinase Pathway Plots ###
############################

# Helper function to create combined heatmap for overlapping pathways
create_combined_heatmap <- function(heatmap_data_list, w_combined, h, save_folder){
  # Find overlapping and union pathways
  all_pathways <- lapply(heatmap_data_list, function(x) rownames(x$matrix))
  common_pathways <- Reduce(intersect, all_pathways)
  union_pathways <- Reduce(union, all_pathways)
  
  # Calculate global color limits across all comparisons
  all_logfc_vals <- unlist(lapply(heatmap_data_list, function(x) {
    matrix_subset <- x$matrix[rownames(x$matrix) %in% common_pathways, , drop = FALSE]
    as.numeric(matrix_subset)
  }))
  all_logfc_vals <- all_logfc_vals[!is.na(all_logfc_vals) & is.finite(all_logfc_vals)]
  
  if (length(all_logfc_vals) > 0) {
    p25_neg <- quantile(all_logfc_vals[all_logfc_vals < 0], probs = 0.25, na.rm = TRUE)
    p75_pos <- quantile(all_logfc_vals[all_logfc_vals > 0], probs = 0.75, na.rm = TRUE)
    limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
    
    if (!is.finite(limit_val) || limit_val == 0) {
      limit_val <- max(abs(all_logfc_vals), na.rm = TRUE)
    }
    
    global_color_limits <- c(-limit_val, limit_val)
  } else {
    global_color_limits <- c(-1, 1)
  }
  
  col_fun <- colorRamp2(c(global_color_limits[1], 0, global_color_limits[2]), 
                        c("#0000CC", "#BDBDBD", "#C62828"))
  
  
  
  if (length(common_pathways) > 0) {
    
    # Create individual heatmaps for each comparison (common pathways only)
    individual_heatmaps <- list()
    
    for (comp_name in names(heatmap_data_list)) {
      # Get matrix subset for common pathways only
      comp_matrix <- heatmap_data_list[[comp_name]]$matrix
      comp_matrix_subset <- comp_matrix[rownames(comp_matrix) %in% common_pathways, , drop = FALSE]

      comp_matrix_subset <- comp_matrix_subset[
        !apply(comp_matrix_subset == 0, 1, all),
        !apply(comp_matrix_subset == 0, 2, all),
        drop = FALSE
      ]
      
      # Ensure same pathway order across all comparisons
      comp_matrix_subset <- comp_matrix_subset[common_pathways, , drop = FALSE]
      comp_matrix_subset[is.na(comp_matrix_subset)] <- 0
      
      # Create heatmap for this comparison
      ht_individual <- Heatmap(
        comp_matrix_subset,
        name = "LogFC",
        col = col_fun,
        cluster_rows = F,
        cluster_columns = FALSE,  # Order columns by name instead
        column_order = order(colnames(comp_matrix_subset)),
        row_order = order(rownames(comp_matrix_subset)),    
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = comp_name,
        column_title_gp = gpar(fontsize = 8, fontface = "bold"),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(8, "cm"),
        rect_gp = gpar(col = "black", lwd = 0.5),
        show_heatmap_legend = FALSE  # Only show legend on first heatmap
      )
      
      individual_heatmaps[[comp_name]] <- ht_individual
    }
    
    # Show legend only on the first heatmap
    individual_heatmaps[[1]]@heatmap_param$show_heatmap_legend <- TRUE
    individual_heatmaps[[1]]@heatmap_param$heatmap_legend_param <- list(
      title = "LogFC", 
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
    
    # Combine heatmaps horizontally
    ht_combined <- Reduce(`+`, individual_heatmaps)
    
    # Calculate dimensions based on number of comparisons
    n_comparisons <- length(heatmap_data_list)
    n_pathways_combined <- length(common_pathways)
    
    # Adjust width based on number of comparisons
    if (!is.null(w_combined)){
      if (n_comparisons == 2) {
        w_combined <- 25
      } else if (n_comparisons == 3) {
        w_combined <- 35
      } else {
        w_combined <- 45
      }
    }
    
    
    
    # Save combined heatmap
    combined_filename <- paste0(save_folder, "/Kinase_Pathway_heatmap_overlapping_comparisons.png")
    png(combined_filename, 
        width = w_combined, 
        height = h, 
        units = "cm",
        res = 300)
    draw(ht_combined)
    dev.off()
    
    cat("Saved combined heatmap:", combined_filename, "\n")
    cat("Found", length(common_pathways), "overlapping pathways across", 
        length(heatmap_data_list), "comparisons\n")
    
  } else {
    cat("No overlapping pathways found across comparisons\n")
  }
}

# Function to create combined heatmap using union of pathways (for exactly 2 comparisons)
create_union_combined_heatmap <- function(heatmap_data_list, w_combined, h, save_folder) {
  # Only works with exactly 2 comparisons
  if (length(heatmap_data_list) != 2) {
    cat("Union combined heatmap only works with exactly 2 comparisons. Found", 
        length(heatmap_data_list), "comparisons.\n")
    return()
  }
  
  # Get pathway sets for each comparison
  comp_names <- names(heatmap_data_list)
  pathways_comp1 <- rownames(heatmap_data_list[[1]]$matrix)
  pathways_comp2 <- rownames(heatmap_data_list[[2]]$matrix)
  
  # Find common and unique pathways
  common_pathways <- intersect(pathways_comp1, pathways_comp2)
  unique_comp1 <- setdiff(pathways_comp1, pathways_comp2)
  unique_comp2 <- setdiff(pathways_comp2, pathways_comp1)
  
  # Create ordered pathway list: common, unique_comp1, unique_comp2
  ordered_pathways <- c(common_pathways, unique_comp1, unique_comp2)
  
  # Create pathway category factor for row splitting
  pathway_categories <- c(
    rep("1", length(common_pathways)),
    rep("2", length(unique_comp1)),
    rep("3", length(unique_comp2))
  )
  names(pathway_categories) <- ordered_pathways
  pathway_split <- factor(pathway_categories, levels = c("1", "2", "3"))
  
  cat("Found", length(common_pathways), "common pathways,", 
      length(unique_comp1), "unique to", comp_names[1], ",", 
      length(unique_comp2), "unique to", comp_names[2], "\n")
  
  # Calculate global color limits across all comparisons (including all pathways)
  all_logfc_vals <- unlist(lapply(heatmap_data_list, function(x) {
    as.numeric(x$matrix)
  }))
  all_logfc_vals <- all_logfc_vals[!is.na(all_logfc_vals) & is.finite(all_logfc_vals)]
  
  if (length(all_logfc_vals) > 0) {
    p25_neg <- quantile(all_logfc_vals[all_logfc_vals < 0], probs = 0.25, na.rm = TRUE)
    p75_pos <- quantile(all_logfc_vals[all_logfc_vals > 0], probs = 0.75, na.rm = TRUE)
    limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
    
    if (!is.finite(limit_val) || limit_val == 0) {
      limit_val <- max(abs(all_logfc_vals), na.rm = TRUE)
    }
    
    global_color_limits <- c(-limit_val, limit_val)
  } else {
    global_color_limits <- c(-1, 1)
  }
  
  col_fun <- colorRamp2(c(global_color_limits[1], 0, global_color_limits[2]), 
                        c("#0000CC", "#BDBDBD", "#C62828"))
  
  # Create matrices for each comparison with union of pathways
  individual_heatmaps <- list()
  
  for (i in 1:2) {
    comp_name <- comp_names[i]
    comp_matrix <- heatmap_data_list[[comp_name]]$matrix
    
    # Create a new matrix with all pathways, filled with zeros for missing pathways
    union_matrix <- matrix(0, nrow = length(ordered_pathways), ncol = ncol(comp_matrix),
                          dimnames = list(ordered_pathways, colnames(comp_matrix)))
    
    # Fill in the actual values for pathways that exist in this comparison
    existing_pathways <- intersect(ordered_pathways, rownames(comp_matrix))
    union_matrix[existing_pathways, ] <- comp_matrix[existing_pathways, , drop = FALSE]
    
    # Replace NA with 0
    union_matrix[is.na(union_matrix)] <- 0
    
    # Create heatmap for this comparison
    ht_individual <- Heatmap(
      union_matrix,
      name = "LogFC",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,  # Order columns by name instead
      column_order = order(colnames(union_matrix)),
      row_order = order(rownames(union_matrix)),
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_split = pathway_split,  # Split rows by pathway category and cluster within splits
      column_title = comp_name,
      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      row_names_max_width = unit(8, "cm"),
      rect_gp = gpar(col = "black", lwd = 0.5),
      show_heatmap_legend = FALSE  # Only show legend on first heatmap
    )
    
    individual_heatmaps[[comp_name]] <- ht_individual
  }
  
  # Show legend only on the first heatmap
  individual_heatmaps[[1]]@heatmap_param$show_heatmap_legend <- TRUE
  individual_heatmaps[[1]]@heatmap_param$heatmap_legend_param <- list(
    title = "LogFC", 
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  )
  
  # Combine heatmaps horizontally
  ht_combined <- Reduce(`+`, individual_heatmaps)
  
  # Set width for 2 comparisons
  if (is.null(w_combined)) {
    w_combined <- 35  # Default for union heatmap with 2 comparisons
  }
  
  # Save combined union heatmap
  union_filename <- paste0(save_folder, "/Kinase_Pathway_heatmap_union_comparisons.png")
  png(union_filename, 
      width = w_combined, 
      height = h, 
      units = "cm",
      res = 300)
  draw(ht_combined)
  dev.off()
  
  cat("Saved union combined heatmap:", union_filename, "\n")
}

# Function to create kinase-pathway heatmaps from network analysis results
plot_kinase_pathway_heatmaps <- function(result_folder, save_folder = NULL, w = 20, h = 14, w_combined = 25) {
  
  if (is.null(save_folder)) {
    save_folder <- file.path(result_folder, "kinase_pathway_heatmaps")
  }
  
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }
  
  # Find all pathways files (excluding *_all.csv files)
  pathways_files <- list.files(result_folder, 
                              pattern = "pathways_.*\\.csv$", 
                              full.names = TRUE, recursive = TRUE)
  
  # Exclude pathways_*_all.csv files
  pathways_files <- pathways_files[!grepl("pathways_.*_all\\.csv$", pathways_files)]
  
  if (length(pathways_files) == 0) {
    stop("No pathways_*.csv files found in the specified folder (excluding *_all.csv files)")
  }
  
  # Extract spec_cutoff and comparison names from filenames
  extract_file_info <- function(filepath) {
    basename_file <- basename(filepath)
    # Remove pathways_ prefix and .csv suffix
    full_comparison <- gsub("^pathways_(.*)\\.csv$", "\\1", basename_file)
    
    # Extract spec cutoff
    spec_match <- regmatches(full_comparison, regexpr("_spec[0-9.]+$", full_comparison))
    if (length(spec_match) > 0) {
      spec_cutoff <- gsub("_spec", "", spec_match)
      comparison_name <- gsub("_spec[0-9.]+$", "", full_comparison)
    } else {
      spec_cutoff <- "unknown"
      comparison_name <- full_comparison
    }
    
    return(list(spec_cutoff = spec_cutoff, comparison = comparison_name, full_comparison = full_comparison))
  }
  
  # Get file information for all files
  file_info_list <- lapply(pathways_files, extract_file_info)
  names(file_info_list) <- pathways_files
  
  # Group files by spec_cutoff
  spec_cutoffs <- unique(sapply(file_info_list, function(x) x$spec_cutoff))
  cat("Found spec_cutoffs:", paste(spec_cutoffs, collapse = ", "), "\n")
  
  # Process each spec_cutoff separately
  for (current_spec_cutoff in spec_cutoffs) {
    cat("\n=== Processing spec_cutoff:", current_spec_cutoff, "===\n")
    
    # Get files for this spec_cutoff
    files_for_cutoff <- names(file_info_list)[sapply(file_info_list, function(x) x$spec_cutoff == current_spec_cutoff)]
    
    if (length(files_for_cutoff) == 0) {
      cat("No files found for spec_cutoff:", current_spec_cutoff, "\n")
      next
    }
    
    # Create spec_cutoff-specific save folder
    spec_save_folder <- file.path(save_folder, paste0("spec_", current_spec_cutoff))
    if (!dir.exists(spec_save_folder)) {
      dir.create(spec_save_folder, recursive = TRUE)
    }
    
    # Extract comparison names for this spec_cutoff
    comparison_names <- sapply(files_for_cutoff, function(f) file_info_list[[f]]$comparison)
    names(files_for_cutoff) <- comparison_names
    
    # Process each comparison for this spec_cutoff
    heatmap_data_list <- list()
    
    for (i in seq_along(files_for_cutoff)) {
      comparison <- names(files_for_cutoff)[i]
      pw_file <- files_for_cutoff[i]
      
      cat("Processing comparison:", comparison, "(spec_cutoff:", current_spec_cutoff, ")\n")
      
      # Find corresponding nodes file (without pathways)
      # Convert from pathways_comparison_spec*.csv to nodes_comparison_spec*.csv
      nodes_file_basename <- gsub("^pathways_", "nodes_", basename(pw_file))
      nodes_file <- file.path(dirname(pw_file), nodes_file_basename)
      
      if (!file.exists(nodes_file)) {
        warning("Could not find corresponding nodes file for: ", comparison, " (spec_cutoff: ", current_spec_cutoff, ")")
        next
      }

      # Read the files
      tryCatch({
        pathways_data <- read_csv(pw_file, show_col_types = FALSE)
        nodes <- read_csv(nodes_file, show_col_types = FALSE)
        
        # Process pathways data: use Group_Pw and Genes columns, separate genes by ";"
        nodes_pw_expanded <- pathways_data %>%
          select(Group_Pw, Genes) %>%
          separate_longer_delim(Genes, delim = ";") %>%
          filter(!is.na(Genes), Genes != "", !is.na(Group_Pw), Group_Pw != "") %>%
          rename(pathway = Group_Pw, Protein = Genes) %>%
          left_join(nodes %>% select(Protein, type, LogFC), 
                   by = "Protein") %>%
          filter(type %in% c("Kinase", "Artificial")) %>%
          select(Protein, pathway, LogFC) %>%
          distinct()
        
        if (nrow(nodes_pw_expanded) == 0) {
          warning("No kinase data found for comparison: ", comparison, " (spec_cutoff: ", current_spec_cutoff, ")")
          next
        }
        
        # Create matrix for heatmap
        heatmap_matrix <- nodes_pw_expanded %>%
          pivot_wider(names_from = Protein, values_from = LogFC, 
                     values_fn = mean, values_fill = NA) %>%
          column_to_rownames("pathway") %>%
          as.matrix()
        
        # Remove rows and columns with all NA
        heatmap_matrix <- heatmap_matrix[!apply(is.na(heatmap_matrix), 1, all), 
                                       !apply(is.na(heatmap_matrix), 2, all), drop = FALSE]
        heatmap_matrix[is.na(heatmap_matrix)] <- 0
        
        if (ncol(heatmap_matrix) == 0 || nrow(heatmap_matrix) == 0) {
          warning("Empty matrix for comparison: ", comparison, " (spec_cutoff: ", current_spec_cutoff, ")")
          next
        }
        
        # Store for later combined analysis
        heatmap_data_list[[comparison]] <- list(
          matrix = heatmap_matrix,
          data = nodes_pw_expanded
        )
        
        # Apply the same color logic as visualize_network_pg
        logfc_vals <- as.numeric(heatmap_matrix)
        logfc_vals <- logfc_vals[!is.na(logfc_vals) & is.finite(logfc_vals)]
        
        if (length(logfc_vals) > 0) {
          p25_neg <- quantile(logfc_vals[logfc_vals < 0], probs = 0.25, na.rm = TRUE)
          p75_pos <- quantile(logfc_vals[logfc_vals > 0], probs = 0.75, na.rm = TRUE)
          limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
          
          if (!is.finite(limit_val) || limit_val == 0) {
            limit_val <- max(abs(logfc_vals), na.rm = TRUE)
          }
          
          color_limits <- c(-limit_val, limit_val)
        } else {
          color_limits <- c(-1, 1)
        }
        
        # Color function (blue -> gray -> red)
        col_fun <- colorRamp2(c(color_limits[1], 0, color_limits[2]), 
                             c("#0000CC", "#BDBDBD", "#C62828"))
        
        # Calculate reasonable dimensions
        n_kinases <- ncol(heatmap_matrix)
        n_pathways <- nrow(heatmap_matrix)
        
        # Base width and height with reasonable limits
        if (!is.null(w)){
          w_adj <- max(20, max(8, n_kinases * 0.4))  # Max reasonable width
        } else {
          w_adj <- w
        }
        if (!is.null(h)){
          h_adj <- max(14, max(4, n_pathways * 0.25 + 2))  # Extra space for row names
        } else {
          h_adj <- h
        }
        
        
        # Create heatmap
        ht <- Heatmap(
          heatmap_matrix,
          name = "LogFC",
          col = col_fun,
          cluster_rows = TRUE,  # cluster pathways
          cluster_columns = TRUE,  # Cluster kinases
          show_row_names = TRUE,
          show_column_names = TRUE,
          column_title = paste("Kinase-Pathway Heatmap:", comparison, "(spec:", current_spec_cutoff, ")"),
          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          row_names_gp = gpar(fontsize = 6),  # Smaller font for long pathway names
          column_names_gp = gpar(fontsize = 8),
          row_names_max_width = unit(8, "cm"),  # More space for row names
          rect_gp = gpar(col = "black", lwd = 0.5),
          heatmap_legend_param = list(title = "LogFC", 
                                     title_gp = gpar(fontsize = 10),
                                     labels_gp = gpar(fontsize = 8))
        )
        
        # Save individual heatmap
        filename <- paste0(spec_save_folder, "/Kinase_Pathway_heatmap_", 
                          gsub("[^A-Za-z0-9_-]", "_", comparison), ".png")
        png(filename, 
            width = w_adj, 
            height = h_adj,
            units = "cm",
            res = 300)
        draw(ht)
        dev.off()
        
        cat("Saved heatmap:", filename, "\n")
        
      }, error = function(e) {
        warning("Error processing ", comparison, " (spec_cutoff: ", current_spec_cutoff, "): ", e$message)
      })
    }
    
    # Create combined heatmaps for this spec_cutoff
    if (length(heatmap_data_list) >= 2 && length(heatmap_data_list) <= 4) {
      
      create_combined_heatmap(heatmap_data_list = heatmap_data_list, w_combined = w_combined, h = h, 
                              save_folder = spec_save_folder)
      
      # Create union combined heatmap for exactly 2 comparisons
      if (length(heatmap_data_list) == 2) {
        create_union_combined_heatmap(heatmap_data_list = heatmap_data_list, w_combined = w_combined, h = h, 
                                     save_folder = spec_save_folder)
      }
      
    } else if (length(heatmap_data_list) > 4) {
      cat("Too many comparisons (", length(heatmap_data_list), ") for combined heatmap in spec_cutoff ", current_spec_cutoff, ". Maximum 4 allowed. Skipping combined heatmap.\n")
    } else if (length(heatmap_data_list) == 1) {
      cat("Only one comparison found for spec_cutoff ", current_spec_cutoff, ". Skipping combined heatmap.\n")
    }
    
    cat("Heatmap generation completed for spec_cutoff:", current_spec_cutoff, ". Files saved in:", spec_save_folder, "\n")
  }
  
  cat("All heatmap generation completed. Files saved in:", save_folder, "\n")
}

