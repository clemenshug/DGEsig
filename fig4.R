library(synExtra)
library(tidyverse)
library(cmapR)
library(here)
library(egg)
library(broom)
library(ggrepel)

synapser::synLogin()
syn <- synExtra::synDownloader("data")

wd <- here("fig4")
dir.create(wd, showWarnings = FALSE)

theme_set(theme_bw())

compound_name_map <- synapser::synGet("syn22035396", version = 3) %>%
  chuck("path") %>%
  read_rds() %>%
  filter(fp_name == "morgan_normal") %>%
  chuck("data", 1) %>%
  distinct(
    lspci_id,
    name = str_to_lower(name)
  )

cmap_signatures_profiled <- syn("syn21747571") %>%
  read_rds()

cmap_gene_meta <- syn("syn21547102") %>%
  read_csv()

clue_res_dge <- syn("syn21907139") %>%
  read_rds()

clue_res_l1000 <- syn("syn21907143") %>%
  read_rds()

clue_res_combined <- syn("syn21907166") %>%
  read_rds()

diff_exp_by_conc <- syn("syn21559856") %>%
  read_rds()

diff_exp_linear <- syn("syn21559859") %>%
  read_rds()

pertubation_meta <- syn("syn21547097") %>%
  read_csv()

signature_meta <- syn("syn21547101") %>%
  read_csv()

ensembl_hugo_map <- diff_exp_by_conc %>%
  pull("result") %>%
  map("ensembl_gene_id") %>%
  reduce(union) %>%
  genebabel::query_hgnc("ensembl_gene_id")

etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
  theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
        axis.text.y = etxt(12), axis.title = etxt(14),
        legend.text = etxt(12), legend.title = etxt(14),
        title = etxt(12))
}

calculate_dose_response_l1000 <- function(meta, signatures, z_threshold = 1.645) {
  select_signatures <- signatures %>%
    select(pr_gene_id, one_of(meta[["sig_id"]])) %>%
    filter(
      select(., -pr_gene_id) %>%
        mutate_all(~abs(.x) >= z_threshold) %>%
        purrr::reduce(magrittr::or)
    )
  select_signatures %>%
    gather("sig_id", "z_score", -pr_gene_id) %>%
    inner_join(
      select(meta, sig_id, pert_dose), by = "sig_id"
    ) %>%
    mutate(change_norm = z_score) %>%
    rename(change = z_score, drug_conc = pert_dose) %>%
    group_nest(pr_gene_id)
}

calculate_dose_response_dge <- function(data, p_threshold = 0.05) {
  selected_signatures <- data %>%
    select(condition, result) %>%
    unnest(result) %>%
    # mutate(log2FoldChange = log2FoldChange_MLE) %>%
    group_by(condition) %>%
    mutate(change_norm = scale(log2FoldChange)) %>%
    ungroup() %>%
    group_by(ensembl_gene_id) %>%
    filter(any(na.omit(padj) <= p_threshold)) %>%
    ungroup() %>%
    select(condition, ensembl_gene_id, log2FoldChange, change_norm)
  selected_signatures %>%
    inner_join(
      select(data, condition, drug_conc), by = "condition"
    ) %>%
    rename(change = log2FoldChange) %>%
    drop_na(change) %>%
    group_nest(ensembl_gene_id)
}

calculate_dose_response <- function(drug, cell, time, p_threshold = 0.05, method = "spearman", round_digits = 1) {
  # browser()
  l1000 <- signature_meta %>%
    distinct() %>%
    filter(lspci_id == drug, cell_id == cell, if (!is.null(time)) pert_time == time else TRUE) %>%
    select_if(negate(is.list)) %>%
    calculate_dose_response_l1000(cmap_signatures_profiled, z_threshold = qnorm(1 - p_threshold)) %>%
    inner_join(cmap_gene_meta %>% filter(pr_is_lm == 1) %>% select(pr_gene_id, symbol), by = "pr_gene_id")
  dge <- diff_exp_by_conc %>%
    filter(lspci_id == !!drug, cells == cell, if (!is.null(time)) replace_na(time == !!time, TRUE) else TRUE) %>%
    calculate_dose_response_dge(p_threshold = p_threshold) %>%
    inner_join(select(ensembl_hugo_map, ensembl_gene_id, symbol), by = "ensembl_gene_id")
  list(
    l1000 = l1000,
    dge = dge
  ) %>%
    bind_rows(.id = "method") %>%
    mutate(
      test = map(
        data,
        ~suppressWarnings(cor.test(
          if (round_digits > 0) round(.x[["change"]], digits = round_digits) else x[["change"]],
          .x[["drug_conc"]], method = !!method
        ))
      )
    ) %>%
    group_by(symbol) %>%
    filter(length(unique(method)) == 2, n() == 2) %>%
    ungroup() %>%
    select(symbol, method, test, data)
}

dose_response_pairwise <- function(df) {
  df %>%
    select(symbol, method, test) %>%
    spread(method, test) %>%
    mutate(
      significance = select(., dge, l1000) %>%
        mutate_all(map_lgl, ~.x[["p.value"]] < 0.05) %>%
        {
          case_when(
            .[["dge"]] & .[["l1000"]] ~ "both",
            .[["dge"]] ~ "dge_only",
            .[["l1000"]] ~ "l1000_only",
            TRUE ~ "neither"
          )
        }
    ) %>%
    mutate_at(vars(dge, l1000), map_dbl, "estimate")
}

dose_response_pairwise_fisher <- function(df) {
  tab <- df %>%
    select(dge, l1000) %>%
    gather("method", "correlation") %>%
    mutate_at(vars(correlation), ~if_else(abs(.x) > 0.25, "correlated", "uncorrelated")) %>%
    table()
  test <- fisher.test(tab, conf.int = TRUE,conf.level = .95)
  list(tab, test)
}

dose_response_curve <- function(df) {
  scale_factor <- df %>%
    group_by(method) %>%
    summarize(max_val = max(abs(change))) %>%
    ungroup() %>%
    spread(method, max_val) %>%
    {.[["l1000"]]/.[["dge"]]}
  df %>%
    mutate(
      change_scaled = if_else(method == "dge", change*scale_factor, change)
    ) %>%
    # ggplot(aes(drug_conc, change_norm, color = method)) +
    ggplot(aes(drug_conc, change_scaled, color = method)) +
    geom_hline(yintercept = 0, color = "grey30") +
    geom_point(alpha = 0.5) +
    scale_y_continuous(
      sec.axis = sec_axis(~./scale_factor, name = "DGE")
    ) +
    geom_smooth(aes(fill = method, group = method), method = "lm", alpha = 0.2) +
    scale_x_log10() +
    labs(x = "Dose (uM)", y = "L1000")
}

dose_response_cor_and_curves <- function(df, seed = 42, highlighted_genes = NULL) {
  set.seed(seed)
  dose_cor_plot_data <- dose_response_pairwise(df) %>%
    mutate(
      pos = case_when(
        dge < 0.25 & abs(l1000) < 0.25 ~ "left_middle",
        dge < 0.25 & l1000 < 0.25 ~ "left_bottom",
        dge > 0.25 & abs(l1000) < 0.25 ~ "right_middle",
        dge > 0.25 & l1000 > 0.25 ~ "right_top",
        TRUE ~ "no"
      )
    ) %>%
    group_by(pos) %>%
    mutate(
      selected = if (unique(pos) == "no")
        "no"
      else if (!is.null(highlighted_genes))
        if_else(symbol %in% highlighted_genes, pos, "no")
      else
        sample(c(rep_len("no", n() - 1), unique(pos)), n())
    ) %>%
    ungroup()
  dose_fisher <- dose_response_pairwise_fisher(dose_cor_plot_data)
  cor_plot <- dose_cor_plot_data %>%
    mutate(
      color_selected = if_else(selected == "no", "black", "red") %>%
        factor(levels = c("black", "red"))
    ) %>%
    arrange(color_selected) %>%
    ggplot(aes(dge, l1000, color = color_selected)) +
    geom_hline(yintercept = c(-0.25, 0.25), color = "grey30") +
    geom_vline(xintercept = c(-0.25, 0.25), color = "grey30") +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "yellow", alpha = 0.2)),
      xmin = -1, xmax = -0.25, ymin = -1, ymax = -0.25
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "yellow", alpha = 0.2)),
      xmin = 0.25, xmax = 1, ymin = 0.25, ymax = 1
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "magenta", alpha = 0.2)),
      xmin = -0.25, xmax = 0.25, ymin = 0.25, ymax = 1
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "magenta", alpha = 0.2)),
      xmin = -0.25, xmax = 0.25, ymin = -1, ymax = -0.25
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "cyan", alpha = 0.2)),
      xmin = -1, xmax = -0.25, ymin = -0.25, ymax = 0.25
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "cyan", alpha = 0.2)),
      xmin = 0.25, xmax = 1, ymin = -0.25, ymax = 0.25
    ) +
    coord_cartesian(expand = FALSE) +
    geom_point(aes(size = selected), alpha = 0.8) +
    scale_size_manual(
      values = c("no" = 2, "left_middle" = 3, "left_bottom" = 3, "right_middle" = 3, "right_top" = 3),
      guide = FALSE
    ) +
    geom_text_repel(
      aes(label = symbol),
      data = ~.x %>%
        mutate(symbol = if_else(selected == "no", "", symbol)),
      # box.padding = 0.5,
      # point.padding = 0.5,
      max.iter = 8000
    ) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    scale_color_identity() +
    guides(color = FALSE) +
    labs(x = "Dose-response correlation DGE", y = "Dose-respose correlation L1000")
  dose_plots <- df %>%
    inner_join(
      dose_cor_plot_data %>%
        filter(selected != "no"),
      by = "symbol"
    ) %>%
    select(symbol, pos, method, data) %>%
    unnest(data) %>%
    group_nest(symbol, pos) %>%
    mutate(
      plot = map(data, dose_response_curve) %>%
        map2(symbol, ~.x + labs(title = .y))
    ) %>%
    {set_names(.[["plot"]], .[["pos"]])}
  arrangeGrob(
    grobs = list(
      cor_plot +
        labs(title = "Palbociclib dose-response correlation"),
      dose_plots[["left_middle"]] +
        theme(legend.position = "none"),
      dose_plots[["left_bottom"]] +
        theme(legend.position = "none"),
      dose_plots[["right_middle"]],
      dose_plots[["right_top"]]
    ) %>%
      map(~.x + theme_bold()) %>%
      c(
        list(
          grid::nullGrob(),
          grid::nullGrob()
        )
      ),
    layout_matrix = rbind(
      c(2, 6, 1, 7, 5),
      c(3, 6, 1, 7, 4)
    ),
    widths = unit(c(4, 0.1, 6, 0.1, 5), "in"),
    heights = unit(c(3, 3), "in")
  )
}


palbo_dose_cor <- calculate_dose_response(
  filter(compound_name_map, name == "palbociclib")$lspci_id,
  "MCF7", time = 24
)

picked_genes <- c(
  "FAM20B", "VPS28", "EGR1", "HPRT1"
)

palbo_dose_cor_plot <- dose_response_cor_and_curves(
  palbo_dose_cor, seed = 1, highlighted_genes = picked_genes
)

grid::grid.draw(palbo_dose_cor_plot)

ggsave(
  file.path(wd, "palbo_dose_cor_plot.png"),
  palbo_dose_cor_plot, width = 15.5, height = 6
)

palbo_dose_cor_plot_data <- dose_response_pairwise(palbo_dose_cor)

palbo_dose_fisher <- dose_response_pairwise_fisher(palbo_dose_cor_plot_data)

palbo_dose_cor_plot <- palbo_dose_cor_plot_data %>%
  ggplot(aes(dge, l1000)) +
  geom_hline(yintercept = c(-0.25, 0.25), color = "grey30") +
  geom_vline(xintercept = c(-0.25, 0.25), color = "grey30") +
  geom_point(alpha = 0.8) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(title = "Palbociclib dose-response correlation per gene")

ggsave(
  file.path(wd, "palbo_dose_cor_plot_significance.pdf"),
  palbo_dose_cor_plot,
  width = 6, height = 5
)


palbo_dose_srsf3 <- palbo_dose_cor %>%
  select(symbol, method, data) %>%
  filter(symbol == "SRSF3") %>%
  unnest(data) %>%
  dose_response_curve()

ggsave(
  file.path(wd, "palbo_dose_srsf3.pdf"),
  palbo_dose_srsf3,
  width = 6, height = 4
)