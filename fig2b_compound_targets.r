library( tidyverse )
library( seriation )   # For optimal leaf reordering
library(synapser)
library(here)
library(qs)
library(magrittr)
library(data.table)
library(powerjoin)
library(seriation)

pathData <- "~/data"

wd <- here("fig2")
dir.create(wd, showWarnings = FALSE)

synapser::synLogin()
syn <- synExtra::synDownloader(pathData, .cache = TRUE)

rdbu_colormap <- RColorBrewer::brewer.pal(5, "RdBu")

cmap_nonlinear_colormap <- list(
  colors = c(
    rdbu_colormap[1:3],
    rdbu_colormap[3:5]
  ) %>%
    rev(),
  values = scales::rescale(c(-100, -90, -80, 80, 90, 100), from = c(-100, 100))
)

## Load all results
R_all <- syn("syn26468923") %>%
  # file.path(pathData, "clue_results_combined.rds") %>%
  qread()

R_cells <- R_all %>%
  filter(result_type == "pert", score_level == "cell") %>%
  pluck("data", 1) %>%
  as_tibble()

R <- R_all %>%
  # filter( result_type == "pert", score_level == "cell" ) %>%
  # pluck( "data", 1 ) %>%
  # filter(cell_id == "MCF7")
  filter( result_type == "pert", score_level == "summary" ) %>%
  pluck( "data", 1 )

condition_conc_vars <- c("cells", "drug_id", "lspci_id", "stim", "stim_conc", "time")

M <- syn("syn25292310") %>%
  # syn("syn22000707") %>%
  qread() %>%
  unnest(meta)

cmap_gene_sets <- syn("syn25314203.4") %>%
  qread()

dge_gene_sets <- syn("syn25303778") %>%
  qread()

R_cmap_mcf7 <- syn("syn30410740") %>%
  read_csv() %>%
  filter(cell_id == "MCF7")

R_pcl <- R_all %>%
  filter(result_type == "pcl", score_level == "summary") %>%
  pluck("data", 1)

gene_set_meta <- bind_rows(
  l1000 = cmap_gene_sets %>%
    select(
      cell_aggregate_method, replicate_method,
      lspci_id, drug_conc = pert_dose, time = pert_time, cells = cell_id,
      cutoff, gene_set_id
    ),
  dge = dge_gene_sets %>%
    select(
      concentration_method, replicate_method,
      lspci_id, drug_conc, cells, stim, stim_conc,
      time, gene_set_id
    ),
  .id = "source"
) %>%
  mutate(cells = coalesce(cells, "summary"))

cmap_meta <- syn("syn21547097") %>%
  read_csv(
    col_types = cols(
      cmap_name = col_character(),
      compound_aliases = col_character(),
      target = col_character(),
      moa = col_character()
    )
  )

compound_names <- syn("syn26260344") %>%
  read_csv() %>%
  select(lspci_id, name) %>%
  drop_na() %>%
  bind_rows(
    anti_join(cmap_meta, ., by = "lspci_id") %>%
      select(name = pert_iname, lspci_id) %>%
      drop_na(name)
  ) %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()

M_cell_info <- bind_rows(
  select(cmap_gene_sets, gene_set_id, lspci_id, cells = cell_id),
  select(dge_gene_sets, gene_set_id, lspci_id, cells)
) %>%
  distinct()


gene_sets_with_cmap_results <- gene_set_meta %>%
  filter(
    gene_set_id %in% R_all$data[[4]]$gene_set
  )

## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2 <- gene_set_meta %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    # replicate_method == "replicates_aggregated" | is.na(replicate_method),
    # concentration_method == "concentration_aggregated" | is.na(concentration_method),
    is.na(stim)
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  left_join(
    cmap_meta %>%
      distinct(pert_id, lspci_id_t = lspci_id),
    by = "pert_id"
  ) %>% {
    bind_rows(
      filter(., source == "l1000"),
      # Actually didn't perform cell aggregation for dge data... need to do here
      filter(., source == "dge") %>%
        group_by(source, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_t) %>%
        summarize(
          tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
            {.[order(abs(.))[2]]},
          # tau = mean(tau),
          cells_q = paste(sort(unique(cells_q)), collapse = "_"),
          .groups = "drop"
        ) %>%
        mutate(cell_aggregate_method = "cells_aggregated")
    )
  } %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  group_by(source, cell_aggregate_method, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_q, cells_t, cutoff) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    # tau = mean(tau),
    .groups = "drop"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  ) %>%
  left_join(
    compound_names %>%
      rename(name_t = name),
    by = c("lspci_id_t" = "lspci_id")
  )

cluster_mat <- function(df, target_col = name_t) {
  ## Perform hierarchical clustering on drugT profiles (columns in the final plot)
  ## Use the DGE slice because it is cleaner and less saturated
  R3_mat <- df %>%
    select(name_q, {{target_col}}, tau) %>%
    pivot_wider(names_from = {{target_col}}, values_from = tau) %>%
    column_to_rownames("name_q") %>%
    as.matrix()
  # Just cluster
  DM_rows <- dist(R3_mat)
  clust_rows <- hclust(DM_rows, method = "average") %>%
    reorder(DM_rows, method = "OLO")
  DM_cols <- dist(t(R3_mat))
  clust_cols <- hclust(DM_cols, method = "average") %>%
    reorder(DM_cols, method = "OLO")
  lvl_rows <- clust_rows$labels[clust_rows$order]
  lvl_cols <- clust_cols$labels[clust_cols$order]

  ## Fix the order via factor levels
  R4 <- df %>% mutate(name_q = factor(name_q, lvl_rows),
                      {{target_col}} := factor({{target_col}}, rev(lvl_cols)))
  R4
}

complete_df <- function(df) {
  # Complete missing observations at z-scores that yielded insufficient
  # genes for Clue with NA
  R4_completed <- bind_rows(
    df %>%
      filter(source == "dge"),
    df %>%
      filter(source == "l1000") %>%
      complete(nesting(lspci_id_q, lspci_id_t, name_q, name_t, source), cutoff)
  )
  R4_completed
}


# Start off with data where replicates, cells, and concentrations aggregated and
# using a z-threshold of 0.7
R3 <- R2 %>%
  filter(
    cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    replicate_method == "replicates_aggregated" | is.na(replicate_method),
    concentration_method == "concentration_aggregated" | is.na(concentration_method),
    source == "dge" | (source == "l1000" & cutoff == 0.7)
  ) %>%
  mutate(
    name_q = str_trunc(name_q, 16, ellipsis = "…"),
    name_t = str_trunc(name_t, 16, ellipsis = "…")
  )


# Select which drugs to show

R4 <- R3 %>%
  filter(
    lspci_id_t %in% {
      R3 %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(abs(tau))) %>%
        slice_head(n = 5) %>%
        pull("lspci_id_t")
    }
  )

R4_clustered <- R4 %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()
R4_clustered <- R4 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_clustered %>%
  arrange(name_q, desc(abs(tau)))


## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
  theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
        axis.text.y = etxt(12), axis.title = etxt(14),
        legend.text = etxt(12), legend.title = etxt(14),
        axis.ticks = element_blank(), strip.text = etxt(14))
}

## Plotting a heatmap of clue hits
fplot <- function(X) {
  library(ggtext)
  X_ <- X %>%
    arrange(name_q) %>%
    mutate(
      name_q_color = paste0(
        "<span style=\"color:",
        if_else(
          as.character(name_q) %in% as.character(name_t),
          "#ff0000", "#000000"
        ),
        "\">",
        name_q,
        "</span>"
      ),
      across(name_q_color, fct_inorder)
    )
  # browser()
  ggplot( X_, aes(x=name_t, y=name_q_color, fill=tau) ) +
    theme_minimal() + theme_bold() +
    geom_tile() +
    geom_tile(data=filter(X_, as.character(name_q)==as.character(name_t)), color="black", size=1) +
    scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
    # exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, guide = FALSE, limits = c(-100, 100)) +
    xlab( "CMap Target" ) +
    theme(
      axis.text.y = ggtext::element_markdown()
    )
}

composite_plot <- function(X) {
  ## DGE plot
  gg1 <- fplot( filter(X, source == "dge") ) +
    scale_x_discrete(position = "top") +
    # scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(hjust=0, vjust=0.5),
          plot.margin = margin(r = 0.25, l = 0.25, unit = "in")) +
    ylab( "3' DGE Query" )

  ## L1000 plot
  gg2 <- fplot( filter(X, source == "l1000") ) +
    ylab( "L1000 Query" ) +
    # scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
    # exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, name = "Tau", limits = c(-100, 100)) +
    theme(legend.position = "bottom",
          plot.margin = margin(r = 0.25, l = 0.25, t = .125, unit = "in"))

  ## Create the composite plot
  egg::ggarrange(
    gg1 + coord_equal(),
    gg2 + coord_equal(),
    heights=c(7.5,7), widths = c(7), draw=FALSE
  )
}

p <- composite_plot(R4_clustered)
ggsave(
  file.path(wd, "fig2b_compound_targets_top5_both.pdf"),
  p, width = 22, height = 26
)

# Only positive Tau

# Select which drugs to show
R4 <- R3 %>%
  filter(
    lspci_id_t %in% {
      R3 %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("lspci_id_t")
    }
  )

R4_clustered <- R4 %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()
R4_clustered <- R4 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_clustered %>%
  arrange(name_q, desc(abs(tau)))

p <- composite_plot(R4_clustered)
ggsave(
  file.path(wd, "fig2b_compound_targets_top5_both_pos.pdf"),
  p, width = 22, height = 26
)


# Only targets of DGE queries
R4 <- R3 %>%
  filter(
    lspci_id_t %in% {
      R3 %>%
        filter(source == "dge") %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(abs(tau))) %>%
        slice_head(n = 5) %>%
        pull("lspci_id_t")
    }
  )

R4_clustered <- R4 %>%
  filter(source == "dge") %>%
  cluster_mat()
R4_clustered <- R4 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_clustered$name_t))
  ) %>%
  complete_df()

R4_clustered %>%
  arrange(name_q, desc(abs(tau)))

p <- composite_plot(R4_clustered)
ggsave(
  file.path(wd, "fig2b_compound_targets_top5_dge.pdf"),
  p, width = 15, height = 24
)


# Only targets of DGE queries, only pos
R4 <- R3 %>%
  filter(
    lspci_id_t %in% {
      R3 %>%
        filter(source == "dge") %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("lspci_id_t")
    }
  )

R4_clustered <- R4 %>%
  filter(source == "dge") %>%
  cluster_mat()
R4_clustered <- R4 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_clustered$name_t))
  ) %>%
  complete_df()

R4_clustered %>%
  arrange(name_q, desc(abs(tau)))

p <- composite_plot(R4_clustered)
ggsave(
  file.path(wd, "fig2b_compound_targets_top5_dge_pos.pdf"),
  p, width = 15, height = 24
)
ß

library(ggbeeswarm)

self_similarity_beeswarm_data <- R3 %>%
  # filter(is.na(cutoff) | cutoff == 0.7) %>%
  mutate(
    self_similarity = if_else(
      lspci_id_q == lspci_id_t,
      "self_similarity",
      "cross_similarity"
    ) %>%
      fct_relevel("cross_similarity")
  ) %>%
  arrange(self_similarity)

# This is NOT self-similarity, it's the whole query-target similarity
# ALL targets, violing version
self_similarity_beeswarm_split <- self_similarity_beeswarm_data %>%
  mutate(
    across(
      source,
      \(x) recode(x, "dge" = "DGE", "l1000" = "L1000")
    )
  ) %>%
  ggplot(
    aes(source, tau, color = self_similarity, size = self_similarity)
  ) +
  geom_violin(
    data = ~filter(.x, self_similarity == "cross_similarity"),
    fill = "black",
    color = alpha("black", .8),
    scale = "width"
  ) +
  geom_quasirandom(
    data = ~filter(.x, self_similarity == "self_similarity"),
    shape = 16,
    # size = 1,
    # color = "NA",
    method = "quasirandom",
    bandwidth = 0.2,
    width = 0.45,
    alpha = .5
  ) +
  # geom_beeswarm(
  #     data = ~filter(.x, self_similarity == "cross_similarity"),
  #     shape = 21,
  #     color = "NA",
  #     priority = "random"
  # ) +
  scale_color_manual(
    values = c(
      self_similarity = "#FF0000",
      cross_similarity = "#00000088"
    ),
    guide = FALSE
  ) +
  scale_size_manual(
    values = c(
      self_similarity = 2,
      cross_similarity = 1.5
    ),
    guide = FALSE
  ) +
  # facet_wrap(~source, nrow = 1) +
  theme_light() +
  theme_bold() +
  scale_x_discrete(position = "top") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = .5)
  ) +
  # scale_y_continuous(position = "right", trans = logit_p_trans(-101, 101)) +
  labs(x = "Query source", y = "Tau")

ggsave(
  file.path(wd, "fig2b_compound_targets_self_similarity_beeswarm_agg_side_by_side.pdf"),
  self_similarity_beeswarm_split,
  width = 2.5, height = 4
)

ggsave(
  file.path(wd, "fig2b_compound_targets_self_similarity_beeswarm_agg_linear.pdf"),
  self_similarity_beeswarm_agg_plot +
    scale_y_continuous(position = "right"),
  width = 1.5, height = 7.5
)


# This is NOT self-similarity, it's the whole query-target similarity
# ONLY top 5 targets, beeswarm version

self_similarity_beeswarm_data <- R4 %>%
  # filter(is.na(cutoff) | cutoff == 0.7) %>%
  mutate(
    self_similarity = if_else(
      lspci_id_q == lspci_id_t,
      "self_similarity",
      "cross_similarity"
    ) %>%
      fct_relevel("cross_similarity")
  ) %>%
  arrange(self_similarity)


self_similarity_beeswarm_split <- self_similarity_beeswarm_data %>%
  mutate(
    across(
      source,
      \(x) recode(x, "dge" = "DGE", "l1000" = "L1000")
    )
  ) %>%
  ggplot(
    aes(source, tau, color = self_similarity, size = self_similarity)
  ) +
  geom_quasirandom(
    shape = 16,
    # size = 1,
    # color = "NA",
    method = "quasirandom",
    bandwidth = 0.2,
    width = 0.45,
    alpha = .5
  ) +
  # geom_beeswarm(
  #     data = ~filter(.x, self_similarity == "cross_similarity"),
  #     shape = 21,
  #     color = "NA",
  #     priority = "random"
  # ) +
  scale_color_manual(
    values = c(
      self_similarity = "#FF0000",
      cross_similarity = "#00000088"
    ),
    guide = FALSE
  ) +
  scale_size_manual(
    values = c(
      self_similarity = 2,
      cross_similarity = 1.5
    ),
    guide = FALSE
  ) +
  # facet_wrap(~source, nrow = 1) +
  theme_light() +
  theme_bold() +
  scale_x_discrete(position = "top") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = .5)
  ) +
  # scale_y_continuous(position = "right", trans = logit_p_trans(-101, 101)) +
  labs(x = "Query source", y = "Tau")

ggsave(
  file.path(wd, "fig2b_compound_targets_self_similarity_beeswarm_agg_side_by_side.pdf"),
  self_similarity_beeswarm_split,
  width = 2.5, height = 4
)

ggsave(
  file.path(wd, "fig2b_compound_targets_self_similarity_beeswarm_agg_linear.pdf"),
  self_similarity_beeswarm_agg_plot +
    scale_y_continuous(position = "right"),
  width = 1.5, height = 7.5
)


self_similarity_stats <- R4 %>%
  group_by(source, lspci_id_q, name_q, cutoff) %>%
  arrange(desc(tau)) %>%
  summarize(
    rank = seq_len(n())[lspci_id_q == lspci_id_t],
    rank_percentile = rank / n(),
    .groups = "drop"
  )

wilcox_self_similarities_col_plot <- wilcox_self_similarities %>%
  mutate(log_p = -log10(p.value)) %>%
  ggplot(aes(x = log_p, y = name_q, fill = source)) +
  geom_col() +
  facet_wrap(~source, nrow = 1) +
  guides(fill = FALSE)

ggsave(
  file.path("self_similarity", "wilcox_cols.pdf"),
  wilcox_self_similarities_col_plot,
  width = 4, height = 8
)

wilcox_self_similarities_scatterplot <- wilcox_self_similarities %>%
  mutate(log_p = -log10(p.value)) %>%
  select(idQ, source, p.value) %>%
  spread(source, p.value) %>%
  ggplot(aes(x = DGE, y = L1000)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal()


R2_pcl <- gene_set_meta %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    # replicate_method == "replicates_aggregated" | is.na(replicate_method),
    # concentration_method == "concentration_aggregated" | is.na(concentration_method),
    !cells %in% c("rencells"),
    is.na(stim) | stim == "control"
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R_pcl %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>% {
    bind_rows(
      filter(., source == "l1000"),
      # Actually didn't perform cell aggregation for dge data... need to do here
      filter(., source == "dge") %>%
        group_by(source, replicate_method, concentration_method, lspci_id_q, pert_id, pert_iname, cells_t) %>%
        summarize(
          tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
            {.[order(abs(.))[2]]},
          # tau = mean(tau),
          cells_q = paste(sort(unique(cells_q)), collapse = "_"),
          .groups = "drop"
        ) %>%
        mutate(cell_aggregate_method = "cells_aggregated")
    )
  } %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  # group_by(source, cell_aggregate_method, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_q, cells_t, cutoff) %>%
  # summarize(
  #   tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
  #     {.[order(abs(.))[2]]},
  #   # tau = mean(tau),
  #   .groups = "drop"
  # ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  )


## Plotting a heatmap of clue hits
fplot_pcl <- function(X) {
  ggplot( X, aes(x=pert_id, y=name_q, fill=tau) ) +
    theme_minimal() + theme_bold() +
    geom_tile() +
    # geom_tile(data=filter(X_, as.character(name_q)==as.character(name_t)), color="black", size=1) +
    # scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
    exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, guide = FALSE, limits = c(-100, 100)) +
    xlab( "CMap Target" ) +
    theme(
      axis.text.y = ggtext::element_markdown()
    )
}

composite_plot_pcl <- function(X) {
  ## DGE plot
  gg1 <- fplot_pcl( filter(X, source == "dge") ) +
    scale_x_discrete(position = "top") +
    # scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=0, vjust=0.5),
          plot.margin = margin(r = 0.25, l = 0.25, unit = "in")) +
    ylab( "3' DGE Query" )

  ## L1000 plot
  gg2 <- fplot_pcl( filter(X, source == "l1000") ) +
    ylab( "L1000 Query" ) +
    # scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
    # exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, name = "Tau", limits = c(-100, 100)) +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          plot.margin = margin(r = 0.25, l = 0.25, t = .125, unit = "in"))

  ## Create the composite plot
  egg::ggarrange(
    gg1 + coord_equal(),
    gg2 + coord_equal(),
    heights=c(7.5,7), widths = c(7), draw=FALSE
  )
}


# Start off with data where replicates, cells, and concentrations aggregated and
# using a z-threshold of 0.7
R3_pcl <- R2_pcl %>%
  filter(
    cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    replicate_method == "replicates_aggregated" | is.na(replicate_method),
    concentration_method == "concentration_aggregated" | is.na(concentration_method),
    source == "dge" | (source == "l1000" & cutoff == 0.7)
  ) %>%
  mutate(
    name_q = str_trunc(name_q, 16, ellipsis = "…")
    # pert_id = str_trunc(pert_id, 16, ellipsis = "…")
  )


# Select which drugs to show

R4_pcl <- R3_pcl %>%
  filter(
    pert_id %in% {
      R3_pcl %>%
        group_by(
          source,
          pert_id
        ) %>%
        # arrange(desc(abs(tau))) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 3) %>%
        pull("pert_id")
    }
  )

R4_clustered_pcl <- R4_pcl %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat(target_col = pert_id)

R4_clustered_pcl <- R4_pcl %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_clustered_pcl$name_q)),
    pert_id = factor(pert_id, levels = levels(R4_clustered_pcl$pert_id))
  ) %>%
  drop_na(tau, name_q)

p <- composite_plot_pcl(R4_clustered_pcl)

ggsave(
  file.path(wd, "fig2b_pcl_targets_top3_both_pos.pdf"),
  p, width = 40, height = 40
)



R4_pcl <- R3_pcl %>%
  filter(
    pert_id %in% {
      R3_pcl %>%
        filter(source == "dge") %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        # arrange(desc(abs(tau))) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 3) %>%
        pull("pert_id")
    }
  )

R4_clustered_pcl <- R4_pcl %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat(target_col = pert_id)

R4_clustered_pcl <- R4_pcl %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_clustered_pcl$name_q)),
    pert_id = factor(pert_id, levels = levels(R4_clustered_pcl$pert_id))
  ) %>%
  drop_na(tau, name_q)

p <- composite_plot_pcl(R4_clustered_pcl)

ggsave(
  file.path(wd, "fig2b_pcl_targets_top3_dge_pos.pdf"),
  p, width = 25, height = 40
)

# ggsave(
#   file.path(wd, "fig2b_compound_targets_top5_both.pdf"),
#   p, width = 22, height = 26
# )

# Only positive Tau

# Select which drugs to show
R4 <- R3 %>%
  filter(
    lspci_id_t %in% {
      R3 %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice(5) %>%
        pull("lspci_id_t")
    }
  )

R4_clustered <- R4 %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()
R4_clustered <- R4 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_clustered %>%
  arrange(name_q, desc(abs(tau)))

p <- composite_plot(R4_clustered)
ggsave(
  file.path(wd, "fig2b_compound_targets_top5_both_pos_.9.pdf"),
  p, width = 22, height = 26
)


## MCF7 query


## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2_mcf7_dge <- dge_gene_sets %>%
  drop_na(lspci_id) %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    replicate_method == "replicates_aggregated",
    concentration_method == "concentration_aggregated",
    cells == "MCF7",
    is.na(stim)
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R_cells %>%
      filter(cell_id == "MCF7") %>%
      rename(cells_t = cell_id) %>%
      mutate(
        name_t = pert_iname
      ),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  left_join(
    cmap_meta %>%
      distinct(pert_id, lspci_id_t = lspci_id),
    by = "pert_id"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  )

R2_mcf7_cmap <- R_cmap_mcf7 %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  ) %>%
  mutate(
    cells_q = "MCF7",
    cells_t = "MCF7",
    name_t = pert_iname
  )

R2_mcf7_combined <- bind_rows(
  l1000 = R2_mcf7_cmap,
  dge = R2_mcf7_dge,
  .id = "source"
) %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  # Also aggregates multiple time points for DGE
  group_by(source, replicate_method, concentration_method, pert_type, lspci_id_q, lspci_id_t, cells_q, cells_t, name_q, name_t) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    # tau = mean(tau),
    .groups = "drop"
  )

# Start off with data where replicates, cells, and concentrations aggregated and
# using a z-threshold of 0.7
R3_mcf7 <- R2_mcf7_combined %>%
  drop_na(lspci_id_t) %>%
  filter(
    pert_type == "trt_cp",
    name_t != "forskolin"
  )
  # mutate(
  #   name_q = str_trunc(name_q, 16, ellipsis = "…"),
  #   name_t = str_trunc(name_t, 16, ellipsis = "…")
  # )


# Select which drugs to show

R4_mcf7 <- R3_mcf7 %>%
  filter(
    lspci_id_t %in% {
      R3_mcf7 %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("lspci_id_t")
    }
  )

R4_mcf7_clustered <- R4_mcf7 %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()

R4_mcf7_clustered <- R4_mcf7 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_mcf7_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_mcf7_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_mcf7_clustered %>%
  arrange(name_q, desc(abs(tau)))


p <- composite_plot(R4_mcf7_clustered)
ggsave(
  file.path(wd, "fig2b_mcf7_queries_compound_targets_top5_both_pos.pdf"),
  p, width = 22, height = 26
)


# This is NOT self-similarity, it's the whole query-target similarity
# ALL targets, violing version
p <- R3_mcf7 %>%
  filter(
    lspci_id_q %in% R4_mcf7_clustered$lspci_id_q
  ) %>%
  mutate(
    across(
      source,
      \(x) recode(x, "dge" = "DGE", "l1000" = "L1000")
    )
  ) %>%
  ggplot(
    aes(source, tau)
  ) +
  geom_violin(
    color = NA,
    fill = "grey10",
    width = 1
    # scale = "width"
  ) +
  # geom_quasirandom(
  #   shape = 16,
  #   # size = 1,
  #   # color = "NA",
  #   method = "quasirandom",
  #   bandwidth = 0.2,
  #   width = 0.45,
  #   alpha = .5
  # ) +
  theme_light() +
  theme_bold() +
  scale_x_discrete(position = "top") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = .5)
  ) +
  # scale_y_continuous(position = "right", trans = logit_p_trans(-101, 101)) +
  labs(x = "Query source", y = "Tau")

ggsave(
  file.path(wd, "fig2b_mcf7_queries_compound_targets_beeswarm.pdf"),
  p,
  width = 2.2, height = 3
)


# Absolute Tau

# Select which drugs to show
R4_mcf7 <- R3_mcf7 %>%
  filter(
    lspci_id_t %in% {
      R3_mcf7 %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(abs(tau))) %>%
        slice_head(n = 5) %>%
        pull("lspci_id_t")
    }
  )

R4_mcf7_clustered <- R4_mcf7 %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()
R4_mcf7_clustered <- R4_mcf7 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_mcf7_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_mcf7_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_mcf7_clustered %>%
  arrange(name_q, desc(abs(tau)))

p <- composite_plot(R4_mcf7_clustered)
ggsave(
  file.path(wd, "fig2b_mcf7_queries_compound_targets_top5_both.pdf"),
  p, width = 22, height = 26
)


# Only targets of DGE queries
R4_mcf7 <- R3_mcf7 %>%
  filter(
    lspci_id_t %in% {
      R3_mcf7 %>%
        filter(source == "dge") %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("lspci_id_t")
    }
  )

R4_mcf7_clustered <- R4_mcf7 %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()
R4_mcf7_clustered <- R4_mcf7 %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_mcf7_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_mcf7_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_mcf7_clustered %>%
  arrange(name_q, desc(abs(tau)))

p <- composite_plot(R4_mcf7_clustered)
ggsave(
  file.path(wd, "fig2b_mcf7_queries_compound_targets_top5_dge_pos.pdf"),
  p, width = 22, height = 26
)



## MCF7 queries PCL


## MCF7 query

## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score

R2_mcf7_pcl <- gene_set_meta %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    # replicate_method == "replicates_aggregated" | is.na(replicate_method),
    # concentration_method == "concentration_aggregated" | is.na(concentration_method),
    replicate_method == "replicates_aggregated" | source == "l1000",
    concentration_method == "concentration_aggregated" | source == "l1000",
    (source == "dge" & !cells %in% c(
      "InMyoFib28085", "rencell"
    )) | (source == "l1000" & cells == "summary"),
    cutoff == .7 | source == "dge",
    is.na(stim) | stim == "control"
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R_pcl %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  #   bind_rows(
  #     filter(., source == "l1000"),
  #     # Actually didn't perform cell aggregation for dge data... need to do here
  #     filter(., source == "dge") %>%
  #       group_by(source, replicate_method, concentration_method, lspci_id_q, pert_id, pert_iname, cells_t) %>%
  #       summarize(
  #         tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
  #           {.[order(abs(.))[2]]},
  #         # tau = mean(tau),
  #         cells_q = paste(sort(unique(cells_q)), collapse = "_"),
  #         .groups = "drop"
  #       ) %>%
  #       mutate(cell_aggregate_method = "cells_aggregated")
  #   )
  # } %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  group_by(source, replicate_method, concentration_method, lspci_id_q, pert_id, pert_iname) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    # tau = mean(tau),
    .groups = "drop"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  )

R2_mcf7_pcl %>%
  group_by(source, name_q, pert_id) %>%
  filter(n() > 1) %>%
  arrange(source, name_q, pert_id)

# Start off with data where replicates, cells, and concentrations aggregated and
# using a z-threshold of 0.7
R3_mcf7_pcl <- R2_mcf7_pcl
  # mutate(
  #   name_q = str_trunc(name_q, 16, ellipsis = "…"),
  #   name_t = str_trunc(name_t, 16, ellipsis = "…")
  # )


# Select which drugs to show

R4_mcf7_pcl <- R3_mcf7_pcl %>%
  filter(
    pert_id %in% {
      R3_mcf7_pcl %>%
        filter(source == "dge") %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("pert_id")
    }
  )

R4_mcf7_pcl_clustered <- R4_mcf7_pcl %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat(target_col = pert_id)

R4_mcf7_pcl_clustered <- R4_mcf7_pcl %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_mcf7_pcl_clustered$name_q)),
    pert_id = factor(pert_id, levels = levels(R4_mcf7_pcl_clustered$pert_id))
  ) %>%
  drop_na(tau, name_q)



p <- composite_plot_pcl(R4_mcf7_pcl_clustered)

ggsave(
  file.path(wd, "fig2b_pcl_targets_mcf7_query_top3_dge_pos.pdf"),
  p, width = 30, height = 30
)




# This is NOT self-similarity, it's the whole query-target similarity
# ALL targets, violing version
p <- R3_mcf7_pcl %>%
  filter(
    lspci_id_q %in% R4_mcf7_pcl_clustered$lspci_id_q
  ) %>%
  mutate(
    across(
      source,
      \(x) recode(x, "dge" = "DGE", "l1000" = "L1000")
    )
  ) %>%
  ggplot(
    aes(source, tau)
  ) +
  geom_violin(
    color = NA,
    fill = "grey10",
    width = 1
    # scale = "width"
  ) +
  # geom_quasirandom(
  #   shape = 16,
  #   # size = 1,
  #   # color = "NA",
  #   method = "quasirandom",
  #   bandwidth = 0.2,
  #   width = 0.45,
  #   alpha = .5
  # ) +
  theme_light() +
  theme_bold() +
  scale_x_discrete(position = "top") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = .5)
  ) +
  # scale_y_continuous(position = "right", trans = logit_p_trans(-101, 101)) +
  labs(x = "Query source", y = "Tau")

ggsave(
  file.path(wd, "fig2b_pcl_targets_mcf7_query_beeswarm.pdf"),
  p,
  width = 2.2, height = 3
)

## Rencell queries pcl target

R2_mcf7_pcl <- gene_set_meta %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    # replicate_method == "replicates_aggregated" | is.na(replicate_method),
    # concentration_method == "concentration_aggregated" | is.na(concentration_method),
    replicate_method == "replicates_aggregated" | source == "l1000",
    concentration_method == "concentration_aggregated" | source == "l1000",
    (source == "dge" & cells %in% c(
      "rencell"
    )) | (source == "l1000" & cells == "summary"),
    cutoff == .7 | source == "dge",
    is.na(stim) | stim == "control"
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R_pcl %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  #   bind_rows(
  #     filter(., source == "l1000"),
  #     # Actually didn't perform cell aggregation for dge data... need to do here
  #     filter(., source == "dge") %>%
  #       group_by(source, replicate_method, concentration_method, lspci_id_q, pert_id, pert_iname, cells_t) %>%
  #       summarize(
  #         tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
  #           {.[order(abs(.))[2]]},
  #         # tau = mean(tau),
  #         cells_q = paste(sort(unique(cells_q)), collapse = "_"),
  #         .groups = "drop"
  #       ) %>%
  #       mutate(cell_aggregate_method = "cells_aggregated")
  #   )
  # } %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  group_by(source, replicate_method, concentration_method, lspci_id_q, pert_id, pert_iname) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    # tau = mean(tau),
    .groups = "drop"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  )

R2_mcf7_pcl %>%
  group_by(source, name_q, pert_id) %>%
  filter(n() > 1) %>%
  arrange(source, name_q, pert_id)

# Start off with data where replicates, cells, and concentrations aggregated and
# using a z-threshold of 0.7
R3_mcf7_pcl <- R2_mcf7_pcl
  # mutate(
  #   name_q = str_trunc(name_q, 16, ellipsis = "…"),
  #   name_t = str_trunc(name_t, 16, ellipsis = "…")
  # )


# Select which drugs to show

R4_mcf7_pcl <- R3_mcf7_pcl %>%
  filter(
    pert_id %in% {
      R3_mcf7_pcl %>%
        filter(source == "dge") %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("pert_id")
    }
  )

R4_mcf7_pcl_clustered <- R4_mcf7_pcl %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat(target_col = pert_id)

R4_mcf7_pcl_clustered <- R4_mcf7_pcl %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_mcf7_pcl_clustered$name_q)),
    pert_id = factor(pert_id, levels = levels(R4_mcf7_pcl_clustered$pert_id))
  ) %>%
  drop_na(tau, name_q)


p <- composite_plot_pcl(R4_mcf7_pcl_clustered)

ggsave(
  file.path(wd, "fig2b_pcl_targets_rencell_query_query_top3_dge_pos.pdf"),
  p, width = 30, height = 30
)


# This is NOT self-similarity, it's the whole query-target similarity
# ALL targets, violing version
p <- R3_mcf7_pcl %>%
  filter(
    lspci_id_q %in% R4_mcf7_pcl_clustered$lspci_id_q
  ) %>%
  mutate(
    across(
      source,
      \(x) recode(x, "dge" = "DGE", "l1000" = "L1000")
    )
  ) %>%
  ggplot(
    aes(source, tau)
  ) +
  geom_violin(
    color = NA,
    fill = "grey10",
    width = 1
    # scale = "width"
  ) +
  # geom_quasirandom(
  #   shape = 16,
  #   # size = 1,
  #   # color = "NA",
  #   method = "quasirandom",
  #   bandwidth = 0.2,
  #   width = 0.45,
  #   alpha = .5
  # ) +
  theme_light() +
  theme_bold() +
  scale_x_discrete(position = "top") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = .5)
  ) +
  # scale_y_continuous(position = "right", trans = logit_p_trans(-101, 101)) +
  labs(x = "Query source", y = "Tau")

ggsave(
  file.path(wd, "fig2b_pcl_targets_rencell_query_beeswarm.pdf"),
  p,
  width = 2.2, height = 3
)



## Rencell query


## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2_rencell_dge <- dge_gene_sets %>%
  drop_na(lspci_id) %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    replicate_method == "replicates_aggregated",
    concentration_method == "concentration_aggregated",
    cells == "rencell",
    is.na(stim)
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R %>%
      filter(cell_id == "summary") %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  left_join(
    cmap_meta %>%
      distinct(pert_id, lspci_id_t = lspci_id),
    by = "pert_id"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  )

R2_rencell_cmap <- cmap_gene_sets %>%
  drop_na(lspci_id) %>%
  filter(
    cell_aggregate_method == "cells_aggregated",
    replicate_method == "replicates_aggregated",
    cutoff == .7
  ) %>%
  rename(
    lspci_id_q = lspci_id
  ) %>%
  mutate(
    cells_q = "summary"
  ) %>%
  select(
    -where(is.list),
    -replicate
  ) %>%
  inner_join(
    R %>%
      filter(cell_id == "summary") %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  left_join(
    cmap_meta %>%
      distinct(pert_id, lspci_id_t = lspci_id),
    by = "pert_id"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  )


R2_rencell_combined <- bind_rows(
  l1000 = R2_rencell_cmap,
  dge = R2_rencell_dge,
  .id = "source"
) %>%
  filter(pert_type == "trt_cp") %>%
  power_left_join(
    compound_names %>%
      rename(name_t_lspci_id = name),
    by = c("lspci_id_t" = "lspci_id"),
    check = check_specs(
      unmatched_keys_left = "warn",
      duplicate_keys_right = "warn"
    ),
    na_matches = "never"
  ) %>%
  mutate(
    name_t = coalesce(
      name_t_lspci_id,
      pert_iname
    )
  ) %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  # Also aggregates multiple time points for DGE
  group_by(source, pert_type, lspci_id_q, lspci_id_t, cells_q, cells_t, name_q, name_t) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    # tau = mean(tau),
    .groups = "drop"
  )


R3_rencell <- R2_rencell_combined %>%
  mutate(
    name_q = str_trunc(name_q, 30, ellipsis = "…"),
    name_t = str_trunc(name_t, 30, ellipsis = "…")
  )


# Select which drugs to show

R4_rencell <- R3_rencell %>%
  filter(
    name_t %in% {
      R3_rencell %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("name_t")
    }
  )

R4_rencell %>%
  drop_na(tau) %>%
  filter(source == "dge")  %>%
  group_by(
    name_q, name_t
  ) %>%
  filter(n() > 1)

R4_rencell_clustered <- R4_rencell %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()

R4_rencell_clustered <- R4_rencell %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_rencell_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_rencell_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_rencell_clustered %>%
  arrange(name_q, desc(abs(tau)))


p <- composite_plot(R4_rencell_clustered)
ggsave(
  file.path(wd, "fig2b_rencell_queries_compound_targets_top5_both_pos.pdf"),
  p, width = 40, height = 26
)

# Absolute Tau

# Select which drugs to show

R4_rencell <- R3_rencell %>%
  filter(
    name_t %in% {
      R3_rencell %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(abs(tau))) %>%
        slice_head(n = 5) %>%
        pull("name_t")
    }
  )

R4_rencell %>%
  drop_na(tau) %>%
  filter(source == "dge")  %>%
  group_by(
    name_q, name_t
  ) %>%
  filter(n() > 1)

R4_rencell_clustered <- R4_rencell %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()

R4_rencell_clustered <- R4_rencell %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_rencell_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_rencell_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_rencell_clustered %>%
  arrange(name_q, desc(abs(tau)))


p <- composite_plot(R4_rencell_clustered)
ggsave(
  file.path(wd, "fig2b_rencell_queries_compound_targets_top5_both.pdf"),
  p, width = 40, height = 26
)



# Only targets of DGE queries

R4_rencell <- R3_rencell %>%
  filter(
    name_t %in% {
      R3_rencell %>%
        filter(source == "dge") %>%
        group_by(
          source,
          lspci_id_q
        ) %>%
        arrange(desc(tau)) %>%
        slice_head(n = 5) %>%
        pull("name_t")
    }
  )

R4_rencell %>%
  drop_na(tau) %>%
  filter(source == "dge")  %>%
  group_by(
    name_q, name_t
  ) %>%
  filter(n() > 1)

R4_rencell_clustered <- R4_rencell %>%
  drop_na(tau) %>%
  filter(source == "dge") %>%
  cluster_mat()

R4_rencell_clustered <- R4_rencell %>%
  mutate(
    name_q = factor(name_q, levels = levels(R4_rencell_clustered$name_q)),
    name_t = factor(name_t, levels = levels(R4_rencell_clustered$name_t))
  ) %>%
  drop_na(tau, name_q)

R4_rencell_clustered %>%
  arrange(name_q, desc(abs(tau)))


p <- composite_plot(R4_rencell_clustered)
ggsave(
  file.path(wd, "fig2b_rencell_queries_compound_targets_top5_dge_pos.pdf"),
  p, width = 40, height = 26
)


# This is NOT self-similarity, it's the whole query-target similarity
# ALL targets, violing version
p <- R3_rencell %>%
  filter(
    lspci_id_q %in% R4_rencell_clustered$lspci_id_q
  ) %>%
  mutate(
    across(
      source,
      \(x) recode(x, "dge" = "DGE", "l1000" = "L1000")
    )
  ) %>%
  ggplot(
    aes(source, tau)
  ) +
  geom_violin(
    color = NA,
    fill = "grey10",
    width = 1
    # scale = "width"
  ) +
  # geom_quasirandom(
  #   shape = 16,
  #   # size = 1,
  #   # color = "NA",
  #   method = "quasirandom",
  #   bandwidth = 0.2,
  #   width = 0.45,
  #   alpha = .5
  # ) +
  theme_light() +
  theme_bold() +
  scale_x_discrete(position = "top") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = .5)
  ) +
  # scale_y_continuous(position = "right", trans = logit_p_trans(-101, 101)) +
  labs(x = "Query source", y = "Tau")

ggsave(
  file.path(wd, "fig2b_rencell_queries_compound_targets_beeswarm.pdf"),
  p,
  width = 2.2, height = 3
)
