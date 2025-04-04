library( tidyverse )
library( seriation )   # For optimal leaf reordering
library( cowplot )
library(synExtra)
library(data.table)
library(qs)
library(powerjoin)

synapser::synLogin()

syn <- synDownloader("~/data", .cache = TRUE)

compound_names <- syn("syn26260344") %>%
  read_csv() %>%
  select(lspci_id, name) %>%
  drop_na() %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()

dge_gene_sets <- syn("syn25303778") %>%
  qread()

GS0 <- syn("syn22105667") %>%
  read_csv()

GS <- GS0 %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "replicates_aggregated",
    is.na(stim),
    str_detect(coll(cell_id, ignore_case = TRUE), "rencell")
  ) %>%
  group_by( DrugID=lspci_id ) %>% summarize( Set = list(entrezgene_id) )

R0 <- syn("syn26468923") %>%
  qread()

## Perform additional wrangling
relevant_gene_sets <- dge_gene_sets %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "replicates_aggregated",
    is.na(stim),
    str_detect(coll(cells, ignore_case = TRUE), "rencell")
    # pert_type == "trt_cp"
  )

R <- R0 %>%
  filter(result_type == "pert", score_level == "summary") %>%
  chuck("data", 1) %>%
  inner_join(
    relevant_gene_sets %>%
      select(idQ = lspci_id, gene_set = gene_set_id),
    by = "gene_set"
  ) %>%
  select( idQ, idT = pert_id, tau )

# For GSEA
expr <- syn("syn25303743") %>%
  fread()

expr_agg <- expr[
  replicate_method == "replicates_aggregated" &
    concentration_method == "concentration_aggregated" &
    stim == "" &
    str_detect(coll(cells, ignore_case = TRUE), "rencell")
]

library(fgsea)
library(msigdbr)

# all_msigdbr <- msigdbr(category = "H") %>%
#   bind_rows(msigdbr(category = "C5", subcategory = "BP"))
all_msigdbr <- msigdbr() %>%
  filter(
    str_starts(gs_subcat, fixed("CP")) |
      gs_cat == "H"
  )

hallmark_gene_sets <- all_msigdbr %>%
  group_nest(gs_name) %>% {
    set_names(
      map(.$data, "human_ensembl_gene"),
      .$gs_name
    )
  }

# expr_for_gsea <- expr_agg %>%
#   filter(is.finite(log2FoldChange)) %>%
#   group_nest(lspci_id)

# fgsea_res <- expr_for_gsea %>%
# head(n = 1) %>%
# unnest(data) %>%
fgsea_res <- expr_agg %>%
  filter(is.finite(log2FoldChange)) %>%
  group_by(lspci_id) %>%
  summarize(
    # gsea_res = map(
    #   gene_set_table,
    #   ~fgseaMultilevel(
    #     hallmark_gene_sets, set_names(.x$log2FoldChange, .x$ensembl_gene_id)
    #   )
    # )
    res = {
      message(lspci_id[1])
      fgseaMultilevel(
        hallmark_gene_sets, set_names(log2FoldChange, ensembl_gene_id),
        nproc = 2
      ) %>%
        list()
    },
    .groups = "drop"
  )

qsave(
  fgsea_res,
  "fgsea_res2.qs"
)
fgsea_res <- qread("fgsea_res2.qs")

library(VIM)

fgsea_res_mat <- reduce2(
  fgsea_res$res, fgsea_res$lspci_id,
  function(agg, res, lspci_id) {
    message(lspci_id)
    lspci_id_s <- sym(as.character(lspci_id))
    full_join(
      agg,
      # select(res, pathway, {{lspci_id_s}} := padj),
      select(res, pathway, {{lspci_id_s}} := NES),
      by = "pathway"
    )
  },
  .init = tibble(pathway = character())
) %>%
  column_to_rownames("pathway") %>%
  as.matrix() %>%
  kNN( imp_var = FALSE, trace = TRUE)
# kNN(numFun = laeken::weightedMean, weightDist = TRUE, imp_var = FALSE, trace = TRUE)

library(pheatmap)
withr::with_pdf(
  "fgsea_res_padj.pdf",
  pheatmap(
    fgsea_res_mat %>%
      as.matrix() %>% {
        -log10(.)
      }
  )
)
withr::with_pdf(
  "fgsea_res_nes.pdf",
  pheatmap(
    fgsea_res_mat %>%
      as.matrix()
  )
)

fgsea_cor <- cor(fgsea_res_mat, method = "pearson", use = "all.obs")
# {
#   na_mat <- is.na(.)
#   .[-which(rowSums(na_mat) > 1), -which(colSums(na_mat) > 1)]
# }

fgsea_cor_df <- fgsea_cor %>%
  as_tibble(rownames = "DrugID1") %>%
  pivot_longer(-DrugID1, names_to = "DrugID2", values_to = "Correlation") %>%
  mutate(across(starts_with("DrugID"), as.numeric))

## Compose tau profiles for each drug
## Join against gene sets
V <- R %>% group_by( DrugID=idQ ) %>% arrange(idT) %>%
  summarize( Vals = list(set_names(tau, idT)) ) %>%
  inner_join(GS, by="DrugID")

## Ensure that all tau values are in the same order
stopifnot( map(V$Vals, names) %>% map_lgl(identical, .[[1]]) %>% all )

## Compute pair-wise similarity for all query drugs
tausim <- function( v1, v2 ) cor(v1,v2,method="pearson", use="complete.obs")
jcrdsim <- function(gs1, gs2)
  length(intersect(gs1, gs2)) / length(union(gs1, gs2))
SM <- crossing(rename_all(V, str_c, "1"),
               rename_all(V, str_c, "2")) %>%
  mutate(TauSim  = map2_dbl(Vals1, Vals2, tausim),
         JcrdSim = map2_dbl(Set1, Set2, jcrdsim))

simorder <- function(dfs, .sim) {
  # Common IDs
  ids <- map(dfs, "DrugID1") %>%
    reduce(intersect)

  dfs_subset <- map(
    dfs,
    \(x) filter(
      x,
      DrugID1 %in% ids,
      DrugID2 %in% ids
    )
  )

  # cluster by df 1
  ## Compose the distance matrix
  DM <- dfs_subset[[1]] %>% select( DrugID1, DrugID2, {{.sim}} ) %>%
    spread( DrugID2, {{.sim}} ) %>% as.data.frame() %>%
    column_to_rownames("DrugID1") %>% dist()

  ## Perform hierarchical clustering with optimal leaf reordering
  lvl <- hclust(DM) %>% reorder(DM) %>% {
    .$labels[.$order]
  }

  order_df <- function(x) {
    ## Fix order of rows and column based on clustering results
    x %>%
      power_left_join(
        compound_names %>%
          select(DrugID1 = lspci_id, DrugName1 = name),
        by = "DrugID1",
        check = check_specs(
          unmatched_keys_left = "warn",
          duplicate_keys_right = "warn"
        )
      ) %>%
      power_left_join(
        compound_names %>%
          select(DrugID2 = lspci_id, DrugName2 = name),
        by = "DrugID2",
        check = check_specs(
          unmatched_keys_left = "warn",
          duplicate_keys_right = "warn"
        )
      ) %>%
      transmute(
        DrugID1 = factor(DrugID1, lvl),
        DrugID2 = factor(DrugID2, rev(lvl)),
        Similarity = {{.sim}},
        DrugName1, DrugName2
      ) %>%
      arrange(DrugID1, DrugID2) %>%
      mutate(
        across(
          starts_with("DrugName"),
          fct_inorder
        )
      )
  }
  map(
    dfs_subset,
    order_df
  )
}

## Plotting elements
pal  <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
##palj <- gray.colors(7, start=1, end=0)
palj <- pal[4:7]
ebl <- element_blank
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Plots a similarity matrix
simplot <- function( .df ) {
  ggplot( .df, aes(DrugID1, DrugID2, fill=Similarity) ) +
    theme_minimal() +
    geom_tile(color="gray") +
    ##        geom_rect(aes(xmin=8.5, xmax=13.5, ymin=ndg-8.5+1, ymax=ndg-13.5+1),
    ##                  fill=NA, color="black", size=1) +
    theme(
      axis.text = ebl(),
      axis.title=ebl(),
      legend.text = etxt(12), legend.title=etxt(14),
      plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
}

## Plots a zoom panel
zoomplot <- function( .df ) {
  ggplot( .df, aes(DrugID1, DrugID2, fill=Similarity) ) +
    theme_minimal() + geom_tile(color="gray") +
    scale_y_discrete( labels = function(x) dnm[x] ) +
    theme(axis.title=ebl(), axis.text.x=ebl(), axis.ticks.x=ebl(),
          axis.text.y=etxt(12), plot.background=element_rect(color="black", size=2),
          panel.grid.minor=ebl(), panel.grid.major=ebl())
}

library(data.table)
tas <- syn("syn26260405") %>%
  fread()


tas_weighted_jaccard <- function(data_tas, query_id, min_n = 6) {
  query_tas <- data_tas[lspci_id == query_id, .(gene_id, tas = 11 - tas)]
  data_tas[
    ,
    .(lspci_id, gene_id, tas = 11 - tas)
  ][
    query_tas,
    on = "gene_id",
    nomatch = NULL
  ][
    ,
    mask := tas > 1 | i.tas > 1
  ][
    ,
    if (sum(mask) >= min_n) .(
      "tas_similarity" = sum(pmin(tas[mask], i.tas[mask])) / sum(pmax(tas[mask], i.tas[mask])),
      "n" = sum(mask),
      "n_prior" = .N
    ) else .(
      "tas_similarity" = numeric(),
      "n" = integer(),
      "n_prior" = integer()
    ),
    by = "lspci_id"
  ]
}

tas_used <- tas %>%
  filter(
    lspci_id %in% SM$DrugID1
  )

all_tas_similarity <- tibble(lspci_id_1 = unique(SM$DrugID1)) %>%
  mutate(
    data = map(
      lspci_id_1,
      ~tas_weighted_jaccard(tas_used, .x) %>%
        rename(lspci_id_2 = lspci_id)
    )
  ) %>%
  unnest(data)

sim_dfs <- list(
  Tau = SM %>%
    select(
      DrugID1, DrugID2, Similarity = TauSim
    ),
  Jcrd_genes = SM %>%
    select(
      DrugID1, DrugID2, Similarity = JcrdSim
    ),
  Cor = fgsea_cor_df %>%
    select(
      DrugID1, DrugID2, Similarity = Correlation
    ),
  Jcrd_targets = all_tas_similarity %>%
    select(
      DrugID1 = lspci_id_1, DrugID2 = lspci_id_2, Similarity = tas_similarity
    )
)

sim_dfs_clustered <- simorder(sim_dfs, Similarity)

write_csv(
  tibble(
    DrugID = levels(sim_dfs_clustered$Tau$DrugID1)
  ),
  "fig6_drug_ids.csv"
)

## Plot similarity matrices
ggjcrd_genes <- simplot(sim_dfs_clustered$Jcrd_genes) + scale_fill_gradientn( colors=palj, limits=c(0,1), name="Gene set\nJaccard\nSimilarity" )
ggtau  <- simplot(sim_dfs_clustered$Tau)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Tau\nPearson\nCorrelation" )
ggcor  <- simplot(sim_dfs_clustered$Cor)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="GSEA\nPearson\nCorrelation" )
ggjcrd_targets  <- sim_dfs_clustered$Jcrd_targets %>%
  complete(DrugID1, DrugID2, fill = list(Similarity = NA)) %>%
  simplot()  +
  scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Target\nJaccard\nSimilarity", na.value = "grey80" )


gg <- egg::ggarrange( plots=list(ggtau, ggplot(mtcars) + theme_nothing(), ggjcrd_genes, ggjcrd_targets), ncol=2,
                      labels=c(" A"," B", " C", " D"),
                      label.args = list(gp = grid::gpar(fontfamily = "Helvetica", cex = 2, fontface = 1))) #%>%
##    ggdraw() %>%
##    + draw_plot( gztau, .68, .7, .17, .25 ) %>%
##    + draw_plot( gzjcrd, .18, .7, .17, .25 )
ggsave( "fig6_new.pdf", gg, width=17, height=13 )


remove_na_iter <- function(df, col1, col2, m) {
  df <- complete(df, {{col1}}, {{col2}})
  while(any(is.na(pull(df, {{m}})))) {
    most_missing <- df %>%
      group_by({{col1}}) %>%
      summarize(n_missing = sum(is.na({{m}})), .groups = "drop") %>%
      arrange(desc(n_missing)) %>%
      pull({{col1}}) %>%
      head(n = 1)
    df <- filter(df, !({{col1}} %in% most_missing), !({{col2}} %in% most_missing))
  }
  df
}

all_tas_similarity_no_missing <- all_tas_similarity %>%
  select(
    DrugID1 = lspci_id_1, DrugID2 = lspci_id_2, Similarity = tas_similarity
  ) %>%
  complete(
    DrugID1, DrugID2,
    fill = list(Similarity = NA)
  ) %>%
  remove_na_iter(
    DrugID1, DrugID2, Similarity
  )


sim_dfs2 <- list(
  Tau = SM %>%
    select(
      DrugID1, DrugID2, Similarity = TauSim
    ),
  Jcrd_genes = SM %>%
    select(
      DrugID1, DrugID2, Similarity = JcrdSim
    ),
  Cor = fgsea_cor_df %>%
    select(
      DrugID1, DrugID2, Similarity = Correlation
    ),
  Jcrd_targets = all_tas_similarity_no_missing
)

sim_dfs_clustered2 <- simorder(sim_dfs2, Similarity)

write_csv(
  sim_dfs_clustered2$Tau %>%
    arrange(DrugID1) %>%
    distinct(DrugID = DrugID1, DrugName = DrugName1),
  "fig6_drug_ids2.csv"
)

## Plot similarity matrices
ggjcrd_genes <- simplot(sim_dfs_clustered2$Jcrd_genes) + scale_fill_gradientn( colors=palj, limits=c(0,1), name="Gene set\nJaccard\nSimilarity" )
ggtau  <- simplot(sim_dfs_clustered2$Tau)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Tau\nPearson\nCorrelation" )
ggcor  <- simplot(sim_dfs_clustered2$Cor)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="GSEA\nPearson\nCorrelation" )
ggjcrd_targets  <- sim_dfs_clustered2$Jcrd_targets %>%
  complete(DrugID1, DrugID2, fill = list(Similarity = NA)) %>%
  simplot()  +
  scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Target\nJaccard\nSimilarity", na.value = "grey80" )


gg <- egg::ggarrange( plots=list(ggtau, ggplot(mtcars) + theme_nothing(), ggjcrd_genes, ggjcrd_targets), ncol=2,
                      labels=c(" A"," B", " C", " D"),
                      label.args = list(gp = grid::gpar(fontfamily = "Helvetica", cex = 2, fontface = 1))) #%>%
##    ggdraw() %>%
##    + draw_plot( gztau, .68, .7, .17, .25 ) %>%
##    + draw_plot( gzjcrd, .18, .7, .17, .25 )
ggsave( "fig6_new2.pdf", gg, width=17, height=13 )


