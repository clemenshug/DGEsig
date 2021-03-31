
# library(synExtra)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(here)
library(gt)
library(ComplexUpset)
library(ggvenn)

synapser::synLogin()
syn <- synExtra::synDownloader("~/data/DGE_comp/")

paste_ <- function(...) {
  paste(..., sep = "_")
}

theme_set(theme_light())

wd <- here("fig2")
dir.create(wd, showWarning = FALSE)

meta <- syn("syn22000707") %>%
  read_rds()

pertubation_meta <- syn("syn21547097") %>%
  read_csv(
    col_types = cols(
      lspci_id = col_integer(),
      cmap_name = col_character(),
      compound_aliases = col_character(),
      target = col_character(),
      moa = col_character()
    )
  )

all_gene_sets <- syn("syn22105667") %>%
  read_csv(
    col_types = cols(
      lspci_id = col_integer(),
      cell_aggregate_method = col_character(),
      zscore = col_double(),
      drug_conc = col_double(),
      time = col_double(),
      stim_conc = col_double(),
      stim = col_character(),
      cutoff = col_double(),
      cell_id = col_character(),
      replicate = col_character()
    )
  )

compound_names <- syn("syn22035396") %>%
  read_rds() %>%
  filter(fp_name == "morgan_normal") %>%
  chuck("data", 1) %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()



dge_all_compounds <- meta %>%
  filter(dataset != "fp_transdiff") %>%
  unnest(meta) %>%
  filter(
    # Don't want this cell line from Laura
    !cells %in% c("MDAMB231"),
    # In CDK dataset, only use the four CDK4/6 inhibs, not CDK7
    if_else(
      dataset == "lincs_cdk4_6_7",
      lspci_id %in% c(
        89588, 91464, 78621, 94539
      ) | drug == "control",
      TRUE
    )
  ) %>%
  distinct(dataset, lspci_id) %>%
  drop_na() %>%
  left_join(
    compound_names
  )

good_dge_sets <- all_gene_sets %>%
  filter(source == "dge", padj < 0.1) %>%
  group_by(gene_set_name) %>%
  filter(
    if(sum(direction == "up") >= 10 && sum(direction == "down") >= 10) TRUE else FALSE
  ) %>%
  ungroup()

clue_res_all <- syn("syn21907166") %>%
  read_rds()

clue_res_overlap_query <- clue_res_all %>%
  filter(
    result_type == "pert",
    score_level == "summary"
  ) %>%
  chuck("data", 1)

old_new_lspci_id_map <- clue_res_overlap_query %>%
  distinct(pert_id, lspci_id_old = lspci_id_target) %>%
  drop_na() %>%
  inner_join(
    pertubation_meta %>%
      distinct(pert_id, lspci_id_new = lspci_id) %>%
      drop_na(),
    by = "pert_id"
  )

old_new_lspci_id_map_vec <- with(
  old_new_lspci_id_map,
  set_names(lspci_id_new, lspci_id_old)
)

clue_res_overlap_query_new <- clue_res_overlap_query %>%
  mutate(
    lspci_id_target = old_new_lspci_id_map_vec[as.character(lspci_id_target)],
    lspci_id_query = old_new_lspci_id_map_vec[as.character(lspci_id_query)]
  )

cmap_res_old <- clue_res_overlap_query_new %>%
  filter(query_type == "aggregated", z_score_cutoff %in% c(NA_real_, 0.7)) %>%
  group_by(
    lspci_id_query
  ) %>%
  filter(
    if(any(source == "dge") && any(source == "l1000")) TRUE else FALSE
  ) %>%
  ungroup() %>%
  distinct(lspci_id = lspci_id_query, cmap_query = TRUE) %>%
  drop_na()

old_cmap_upset_data <- meta %>%
  filter(dataset != "fp_transdiff") %>%
  unnest(meta) %>%
  distinct(lspci_id) %>%
  drop_na() %>%
  mutate(
    dge = TRUE
  ) %>%
  left_join(
    pertubation_meta %>%
      filter(dataset != "LINCS_2020") %>%
      distinct(lspci_id) %>%
      drop_na() %>%
      mutate(
        cmap = TRUE
      ),
    by = "lspci_id"
  ) %>%
  left_join(
    pertubation_meta %>%
      distinct(lspci_id) %>%
      drop_na() %>%
      mutate(
        cmap_new = TRUE
      ),
    by = "lspci_id"
  ) %>%
  left_join(
    good_dge_sets %>%
      filter(source == "dge") %>%
      distinct(lspci_id) %>%
      drop_na() %>%
      mutate(
        gene_set = TRUE
      ),
    by = "lspci_id"
  ) %>%
  left_join(
    cmap_res_old,
    by = "lspci_id"
  ) %>%
  mutate(
    across(everything(), replace_na, replace = FALSE)
  )

old_cmap_upset_plot <- upset(
  old_cmap_upset_data,
  intersect = c("dge", "cmap", "cmap_new", "gene_set", "cmap_query"),
  set_sizes = upset_set_size(
    geom = geom_bar()
  ) + 
    geom_text(
      aes(label=stat(count)),
      stat='count',
      color = "white",
      hjust = -0.5
    )
)

cowplot::ggsave2(
  file.path(wd, "fig2a.pdf"),
  old_cmap_upset_plot, width = 6, height = 4
)
