library(tidyverse)
library(DESeq2)
library(synapser)
library(data.table)

synLogin()

syn <- synExtra::synDownloader("~/data/DGE_comp/")

deseq_input <- syn("syn21558154") %>%
  read_rds()

sr_drug_repurposing_two_reps <- deseq_input %>%
  filter(dataset == "lincs_cdk4_6_7") %>%
  pull(meta_deseq) %>%
  map("drug_id") %>%
  {lift_dl(c)(.)} %>%
  table()

two_rep_meta <- deseq_input %>%
  filter(dataset == "lincs_cdk4_6_7") %>%
  mutate(
    meta_deseq = map(
      meta_deseq,
      ~.x %>%
        filter(drug_id %in% sr_drug_repurposing_two_reps)
    )
  ) %>%
  mutate(
    counts_deseq = map2(
      counts_deseq, meta_deseq,
      ~.x[, rownames(.y)]
    )
  )

two_rep_deseq <- two_rep_meta %>%
  transmute(
    dataset,
    data = map2(
      meta_deseq, counts_deseq,
      function(m, counts) {
        m <- m %>%
          rownames_to_column("sample_id")
        m %>%
          filter(drug_id != "control") %>%
          group_nest(condition, control_condition, .key = "treated_meta") %>%
          inner_join(
            m %>%
              filter(drug_id == "control") %>%
              group_nest(control_condition, .key = "control_meta"),
            by = c("control_condition")
          ) %>%
          transmute(
            condition,
            control_condition,
            deseq_meta = map2(treated_meta, control_meta, bind_rows),
            deseq_counts = map2(
              treated_meta, control_meta,
              ~counts[, union(.x[["sample_id"]], .y[["sample_id"]])]
            )
          )
      }
    )
  ) %>%
  unnest(data)

unloadNamespace("synapser")
unloadNamespace("PythonEmbedInR")
deseq_de <- two_rep_deseq %>%
  mutate(
    deseq = map2(
      deseq_meta, deseq_counts,
      ~DESeqDataSetFromMatrix(
        .y,
        .x %>%
          mutate(drug_norm = fct_relevel(drug_norm,  "control")) %>%
          column_to_rownames("sample_id"),
        design = ~drug_norm
      ) %>%
        DESeq()
    )
  )

deseq_res <- deseq_de %>%
  mutate(
    res = map(
      deseq,
      function(x) {
        res_name <- resultsNames(x)[[2]]
        res <- results(x, name = res_name)
        mod <- lfcShrink(x, coef = res_name, res = res, type = "apeglm")
        mod %>%
          as.data.frame() %>%
          rownames_to_column("gene_id") %>%
          full_join(
            res %>%
              as.data.frame() %>%
              rownames_to_column("gene_id") %>%
              select(gene_id, log2FoldChange_MLE = log2FoldChange),
            by = "gene_id"
          )
      }
    )
  )


library(genebabel)
library(clueR)

comp_res_entrez <- deseq_res %>%
  mutate(
    res = map(
      res,
      function(df) {
        df %>%
          dplyr::rename(ensembl_gene_id = gene_id) %>%
          {suppressWarnings(join_hgnc(., "ensembl_gene_id", "ensembl_gene_id", c("symbol", "entrez_id")))} %>%
          dplyr::rename(gene_id = entrez_id)
      }
    )
  )

comp_res_gene_sets <- comp_res_entrez %>%
  filter(
    map_lgl(res, ~nrow(filter(.x, padj < 0.05)) > 20)
  ) %>%
  mutate(
    gene_set = condition %>%
      str_replace_all("\\s", "_") %>%
      str_replace_all("[^\\w]", ""),
    res = map(
      res,
      ~.x %>%
        drop_na(padj) %>%
        arrange(padj) %>%
        filter(padj < 0.05) %>%
        mutate(direction = if_else(log2FoldChange < 0, "down", "up"))
    )
  )

MAX_PER_QUERY = 25
comp_res_clue <- comp_res_gene_sets %>%
  split(ceiling(seq_len(nrow(.)) / MAX_PER_QUERY)) %>%
  map(
    ~bind_rows(set_names(.x[["res"]], .x[["gene_set"]]), .id = "gene_set") %>%
      select(gene_set, gene_id, direction)
  ) %>%
  map(
    clue_gmt_from_df,
    drop_invalid = TRUE
  )

comp_res_jobs <- comp_res_clue %>%
  imap(
    ~clue_query_submit(
      .x[["up"]], .x[["down"]],
      name = paste0("dge_job_", .y),
      use_fast_tool = FALSE
    )
  )

# comp_res_jobs <- c(
#   "5eb305a5d58eba0011f9031c",
#   "5eb305a6ed2e4f0011539a5c",
#   "5eb305a6d58eba0011f9031e",
#   "5eb305a6ed2e4f0011539a5e",
#   "5eb305a6d58eba0011f90320",
#   "5eb305a7ed2e4f0011539a60"
# )


walk(
  comp_res_jobs,
  clue_query_wait
)

comp_res_clue_res <- map_chr(
  comp_res_jobs,
  clue_query_download
)


clue_res <- crossing(
  score_type = "tau",
  result_type = c("pert", "pcl"),
  score_level = c("cell", "summary")
) %>%
  mutate(
    data = pmap(
      .,
      function(...) {
        comp_res_clue_res %>%
          map(clue_parse_result, ...) %>%
          bind_rows()
      }
    )
  )


cmap_signatures <- syn("syn21747571") %>%
  read_rds()

signature_meta <- syn("syn21547101") %>%
  read_csv()

selected_signatures <- signature_meta %>%
  filter(lspci_id %in% sr_drug_repurposing_two_reps)

zscores_selected_signatures <- cmap_signatures %>%
  select(entrez_id, one_of(selected_signatures[["sig_id"]])) %>%
  gather("sig_id", "zscore", -entrez_id) %>%
  inner_join(
    selected_signatures %>%
      select(sig_id, lspci_id, cell_id),
    by = "sig_id"
  )

# Aggregating multiple cells using procedure by 10.1016/j.cell.2017.10.049
zscores_selected_signatures_agg <- zscores_selected_signatures %>%
  as.data.table() %>%
  {
    # First aggregate replicates within same cell line
    .[
      ,
      .(zscore = mean(zscore)),
      keyby = c("entrez_id", "lspci_id", "cell_id")
    ][
      # Then across cell lines
      ,
      .(zscore_agg = quantile(zscore, c(0.67, 0.33), names = FALSE) %>%
          {.[order(abs(.))[2]]}
      ),
      by = c("entrez_id", "lspci_id")
    ]
  } %>%
  as_tibble()

selected_signatures_gene_sets <- tibble(
  cutoff = c(0.7, 0.8, 0.9, 0.95)
) %>%
  mutate(
    data = map(
      cutoff,
      ~zscores_selected_signatures_agg %>%
        filter(abs(zscore_agg) > qnorm(.x)) %>%
        arrange(desc(abs(zscore_agg))) %>%
        transmute(
          lspci_id,
          gene_id = as.character(entrez_id),
          direction = if_else(zscore_agg > 0, "up", "down"),
          gene_set = paste0(lspci_id, "_", .x),
          zscore_agg
        ) %>%
        {suppressWarnings(join_hgnc(., "gene_id", "entrez_id", "symbol"))}
    )
  ) %>%
  unnest(data) %>%
  group_nest(gene_set, lspci_id, cutoff)

two_rep_condition_meta <- two_rep_meta %>%
  pull(meta_deseq) %>%
  bind_rows() %>%
  select(condition, lspci_id, cells, stim, drug_conc) %>%
  distinct() %>%
  as_tibble()

gene_set_comparison <- selected_signatures_gene_sets %>%
  filter(cutoff == 0.7) %>%
  transmute(lspci_id, method = "cmap", data = as.list(data)) %>%
  bind_rows(
    comp_res_gene_sets %>%
      select(condition, data = res) %>%
      left_join(
        two_rep_condition_meta,
        by = "condition"
      ) %>%
      filter(is.na(stim)) %>%
      mutate(method = "dge")
  )

gene_set_comparison_venn <- gene_set_comparison %>%
  group_by(lspci_id) %>%
  summarize(
    venn = set_names(
      map(data, "gene_id"),
      method
    ) %>%
      map(na.omit) %>%
      map(head, n = 300) %>%
      VennDiagram::venn.diagram(filename = NULL) %>%
      list()
  ) %>%
  left_join(
    signature_meta %>%
      group_by(lspci_id) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      select(lspci_id, pert_iname)
  )

pwalk(
  gene_set_comparison_venn,
  function(pert_iname, venn, ...) {
    ggsave(paste0("venn_de_genes_", pert_iname, ".png"), venn)
  }
)



clue_res_l1000 <- syn("syn21907143") %>%
  read_rds()

clue_comparison <-

