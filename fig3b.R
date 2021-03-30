library(tidyverse)
# library(DESeq2)
library(synapser)
library(here)
library(qs)
library(data.table)

synLogin()

syn <- synExtra::synDownloader("~/data/DGE_comp/")

dge_gene_sets <- syn("syn25303778") %>%
  qread()

cmap_gene_sets <- syn("syn25314203") %>%
  qread()

bulk_gene_sets <- syn("syn23554368") %>%
  read_rds() %>%
  rowwise() %>%
  mutate(
    gene_set_table = deseq %>%
      filter(padj < 0.2) %>%
      list()
  ) %>%
  ungroup()

cmap_perturbation_meta <- syn("syn21547097") %>%
  read_csv(
    col_types = list(
      cmap_name = col_character(),
      compound_aliases = col_character(),
      target = col_character(),
      moa = col_character()
    )
  )

cmap_signature_meta <- syn("syn21547101") %>%
  fread()

lspci_id_name_map <- syn("syn22035396") %>%
  read_rds() %>%
  filter(fp_name == "morgan_normal") %>%
  chuck("data", 1L)

gene_set_comparison <- dge_gene_sets %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "per_replicate",
    is.na(stim),
    !is.na(lspci_id)
  ) %>%
  transmute(
    technique = "dge",
    lspci_id, cells,
    replicate = as.factor(paste(dataset, time, replicate, sep = "_")),
    time,
    gene_set = gene_set_table %>%
      map(select, direction, entrezgene_id)
  ) %>%
  bind_rows(
    cmap_gene_sets %>%
      filter(
        replicate_method == "per_replicate",
        cell_aggregate_method == "per_cell_line",
        cutoff == 0.7
      ) %>%
      transmute(
        technique = "l1000", lspci_id, cells = cell_id, replicate,
        gene_set = as.list(gene_set_table) %>%
          map(select, direction, entrezgene_id)
      ) %>%
      inner_join(
        cmap_signature_meta %>%
          distinct(replicate = sig_id, time = pert_time),
        by = "replicate"
      )
  ) %>%
  arrange(technique, desc(time)) %>%
  group_by(technique, lspci_id, cells, time) %>%
  mutate(
    name = paste0(str_to_upper(technique), " ", time, "h Rep ", rev(seq_len(n()))) %>%
      fct_inorder()
  ) %>%
  ungroup()

gene_set_comparison_venn_data <- gene_set_comparison %>%
  mutate(across(cells, str_to_upper)) %>%
  # mutate(gene_set = map(gene_set, group_nest, direction, .key = "gene_set")) %>%
  unnest(gene_set) %>%
  nest(gene_set = c(direction, entrezgene_id)) %>%
  arrange(lspci_id, cells) %>%
  group_by(lspci_id, cells) %>%
  filter(n() > 1) %>%
  summarize(
    # venn = cur_data() %>%
    venn_data = cur_data() %>%
      list(),
    .groups = "drop"
  )



make_upset <- function(df, ...) {
  df %>%
    transmute(
      name,
      gene_set = map(
        gene_set, head, n = 300
      )
    ) %>%
    unnest(gene_set) %>%
    mutate(dummy = TRUE) %>%
    distinct(entrezgene_id, name, dummy) %>%
    pivot_wider(
      id_cols = entrezgene_id,
      names_from = name,
      values_from = dummy,
      values_fill = FALSE
    ) %>%
    select(-entrezgene_id) %>%
    as.data.frame() %>%
    # UpSetR::upset(
    #   nsets = ncol(.),
    #   keep.order = FALSE,
    #   order.by = "freq",
    #   group.by = "degree",
    #   ...
    # )
    ComplexUpset::upset(intersect = colnames(.), ...)
}

gene_set_comparison_upset_data <- gene_set_comparison_venn_data %>%
  left_join(
    cmap_perturbation_meta %>%
      distinct(lspci_id, pert_iname) %>%
      drop_na(),
    by = "lspci_id"
  ) %>%
  filter(map_lgl(venn_data, ~length(unique(.x[["technique"]])) >= 2))

all_combn <- function(x) {
  map(
    seq_along(x),
    ~combn(x, .x, simplify = FALSE)
  ) %>%
    reduce(c)
}

gene_set_comparison_upset <- gene_set_comparison_upset_data %>%
  # head(n = 3) %>%
  mutate(
    venn = map(
      venn_data,
      ~{
        # browser()
        make_upset(
          .x,
          n_intersections = 30,
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size()
            # 'Intersection ratio' = ComplexUpset::intersection_ratio()
          ),
          queries = c(
            map(
              filter(.x, str_starts(name, "DGE")) %>%
                pull(name) %>%
                all_combn(),
              ~{
                # browser()
                ComplexUpset::upset_query(
                  intersect = as.character(.x),
                  only_components = "intersections_matrix",
                  color = "cornflowerblue",
                  fill = "cornflowerblue"
                )
              }
            ),
            map(
              filter(.x, str_starts(name, "L1000")) %>%
                pull(name) %>%
                all_combn(),
              ~ComplexUpset::upset_query(
                intersect = as.character(.x),
                only_components = "intersections_matrix",
                color = "orangered",
                fill = "orangered"
              )
            )
          ),
          set_sizes = FALSE,
          sort_sets = FALSE
        )
      }
    )
  )

pwalk(
  gene_set_comparison_upset,
  function(pert_iname, venn, cells, ...) {
    withr::with_pdf(
      paste0("fig3b_upset_de_genes_", pert_iname, "_", cells, ".pdf"),
      width = 12, height = 6,
      print(venn)
    )
  }
)
