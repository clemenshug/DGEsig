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

# canonical_table <- syn("syn20835543") %>%
#   read_rds() %>%
#   filter(fp_name == "morgan_normal") %>%
#   chuck("data", 1)

gene_set_comparison <- dge_gene_sets %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "per_replicate",
    is.na(stim),
    !is.na(lspci_id)
  ) %>%
  # Keep latest time point
  # group_by(lspci_id, cells, dataset) %>%
  # arrange(desc(time), .by_group = TRUE) %>%
  # slice_head(1) %>%
  # ungroup() %>%
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

gene_set_comparison_venn <- gene_set_comparison_venn_data %>%
  rowwise() %>%
  mutate(
    venn = with(
      venn_data,
      set_names(
        map(gene_set, "entrezgene_id"),
        paste(technique, replicate, sep = "_")
      )
    ) %>%
      map(na.omit) %>%
      map(head, n = 300) %>%
      possibly(VennDiagram::venn.diagram, NULL)(filename = NULL) %>%
      list()
  ) %>%
  ungroup()

pwalk(
  gene_set_comparison_venn,
  function(pert_iname, venn, ...) {
    ggsave(paste0("venn_de_genes_", pert_iname, ".png"), venn)
  }
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
      paste0("upset_de_genes_", pert_iname, "_", cells, "_2.pdf"),
      width = 12, height = 6,
      print(venn)
    )
    # png(paste0("upset_de_genes_", pert_iname, "_", cells, ".png"))
    # print(venn)
    # dev.off()
    # ggsave(paste0("upset_de_genes_", pert_iname, ".png"), venn)
  }
)


alvocidib_venn_data <- gene_set_comparison_upset_data %>%
  filter(
    pert_iname == "alvocidib",
    cells == "MCF7"
  ) %>%
  chuck("venn_data", 1) %>%
  filter(
    !replicate %in% c("lincs_cdk4_6_7_NA_5", "lincs_cdk4_6_7_NA_6")
  ) %>%
  group_by(technique) %>%
  mutate(
    name = paste0(technique, "_", 1:n()) %>%
      as.factor()
  ) %>%
  ungroup() %>%
  arrange(desc(name))

alvocidib_venn <- alvocidib_venn_data %>%
  make_upset(name = "Gene sets", set_sizes = FALSE, keep_empty_groups = FALSE, min_size = 10, sort_sets = FALSE)

ggsave(
  "upset_de_genes_alvocidib_MCF7_select_conditions.png",
  alvocidib_venn,
  width = 8, height = 5
)

alvocidib_clue_gmt <- alvocidib_venn_data %>%
  unnest(gene_set) %>%
  rename(gene_set = name, gene_id = entrezgene_id) %>%
  clueR::clue_gmt_from_df()

alvocidib_clue_job <- clueR::clue_query_submit(
  alvocidib_clue_gmt[["up"]], alvocidib_clue_gmt[["down"]], use_fast_tool = FALSE,
  name = "alvo"
)

alvocidib_clue_dl <- clueR::clue_query_download(alvocidib_clue_job)

alvocidib_clue_res <- clueR::clue_parse_result(
  alvocidib_clue_dl, score_level = "summary", result_type = "pert", score_type = "tau"
) %>%
  mutate(
    gene_set = recode(
      gene_set,
      !!!set_names(
        sort(unique(as.character(alvocidib_venn_data[["name"]]))),
        paste0("X", 1:length(unique(alvocidib_venn_data[["name"]])))
      )
    )
  )

alvocidib_clue_cor <- alvocidib_clue_res %>%
  filter(pert_type == "trt_cp") %>%
  select(pert_id, gene_set, tau) %>%
  spread(gene_set, tau, fill = 0) %>%
  column_to_rownames("pert_id") %>%
  as.matrix() %>%
  cor(method = "spearman")

alvocidib_clue_cor_hm <- alvocidib_clue_cor %>%
  pheatmap::pheatmap()


ensembl <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", version = "98")

all_ensembl_ids <- bulk_gene_sets %>%
  pull(gene_set_table) %>%
  map("ensembl_gene_id") %>%
  reduce(union)

gene_id_mapping <- biomaRt::getBM(
  c("entrezgene_id", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = all_ensembl_ids,
  mart = ensembl
)

unloadNamespace("biomaRt")
unloadNamespace("AnnotationDBI")

gene_id_mapping_valid <- gene_id_mapping %>%
  drop_na()


bulk_drugs_lspci_id_map <- c(
  palbociclib = 91464,
  abemaciclib = 94539,
  ribociclib = 89588
)

bulk_dge_l1000_upset_data <- bulk_drugs_lspci_id_map %>%
  enframe("drug", "lspci_id") %>%
  mutate(
    venn_data = map2(
      lspci_id, drug,
      ~bind_rows(
        cmap_gene_sets %>%
          filter(
            cell_aggregate_method == "per_cell_line",
            replicate_method == "replicates_aggregated",
            lspci_id == .x,
            cell_id == "MCF7",
            cutoff == 0.7
          ) %>%
          transmute(name = paste("cmap", seq_len(n()), sep = "_"), gene_set = as.list(gene_set_table)),
        dge_gene_sets %>%
          filter(
            concentration_method == "concentration_aggregated",
            replicate_method == "replicates_aggregated",
            cells == "MCF7",
            lspci_id == .x,
            is.na(stim)
            # time == 24
          ) %>%
          transmute(name = paste("dge", seq_len(n()), sep = "_"), gene_set = as.list(gene_set_table)),
        bulk_gene_sets %>%
          filter(
            drug == .y,
            cells == "MCF7",
            time == 24
          ) %>%
          rowwise() %>%
          mutate(
            gene_set_table = gene_set_table %>%
              inner_join(gene_id_mapping_valid, by = "ensembl_gene_id") %>%
              arrange(padj) %>%
              list()
          ) %>%
          ungroup() %>%
          transmute(name = paste("bulk", seq_len(n()), sep = "_"), gene_set = gene_set_table)
      )
    )
  )

bulk_dge_l1000_upset <- bulk_dge_l1000_upset_data %>%
  # head(n = 3) %>%
  mutate(
    venn = map(
      venn_data,
      make_upset
    )
  )

pwalk(
  bulk_dge_l1000_upset,
  function(drug, venn, ...) {
    dir.create(here("bulk_dge_cmap_comparison"), showWarnings = FALSE)
    ggsave(
      here("bulk_dge_cmap_comparison", paste0(drug, "_threshold0.2.pdf")),
      venn,
      width = 7, height = 5
    )
  }
)


bulk_dge_l1000_upset_data <- bulk_drugs_lspci_id_map %>%
  enframe("drug", "lspci_id") %>%
  mutate(
    venn_data = map2(
      lspci_id, drug,
      ~bind_rows(
        cmap_gene_sets %>%
          filter(
            cell_aggregate_method == "per_cell_line",
            replicate_method == "per_replicate",
            lspci_id == .x,
            cell_id == "MCF7",
            cutoff == 0.7
          ) %>%
          transmute(
            name = paste("cmap", seq_len(n()), sep = "_"),
            gene_set = map(
              gene_set_table,
              arrange, desc(abs(zscore))
            )
          ),
        dge_gene_sets %>%
          filter(
            concentration_method == "concentration_aggregated",
            replicate_method == "per_replicate",
            cells == "MCF7",
            lspci_id == .x,
            is.na(stim)
            # time == 24
          ) %>%
          transmute(
            name = paste("dge", seq_len(n()), sep = "_"),
            gene_set = map(
              gene_set_table,
              arrange, padj
            )
          ),
        bulk_gene_sets %>%
          filter(
            drug == .y,
            cells == "MCF7",
            time == 24
          ) %>%
          rowwise() %>%
          mutate(
            gene_set_table = gene_set_table %>%
              inner_join(gene_id_mapping_valid, by = "ensembl_gene_id") %>%
              arrange(padj) %>%
              list()
          ) %>%
          ungroup() %>%
          transmute(name = paste("bulk", seq_len(n()), sep = "_"), gene_set = gene_set_table)
      )
    )
  )

bulk_dge_l1000_upset <- bulk_dge_l1000_upset_data %>%
  # head(n = 3) %>%
  mutate(
    venn = map(
      venn_data,
      make_upset
    )
  )

pwalk(
  bulk_dge_l1000_upset,
  function(drug, venn, ...) {
    dir.create(here("bulk_dge_cmap_comparison"), showWarnings = FALSE)
    ggsave(
      here("bulk_dge_cmap_comparison", paste0(drug, "_per_replicate", ".pdf")),
      venn,
      width = 10, height = 5
    )
  }
)

