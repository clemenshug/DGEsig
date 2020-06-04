library( tidyverse )
library( seriation )   # For optimal leaf reordering

pathData <- "~/data/DGEsig"

## Isolate the relevant data chunk
R0 <- file.path(pathData, "clue_results_dge.rds") %>% read_rds() %>%
    pluck( "data", 4 )

## Perform additional wrangling
source("dmap.R")
R <- R0 %>% filter(dataset_query == "sr_repurposing",
                   !is.na(lspci_id_query),
                   !is.na(lspci_id_target),
                   pert_type == "trt_cp" ) %>%
    mutate( drugQ = map_chr(str_split(gene_set, "_"), pluck, 3) ) %>%
    select( drugQ, drugT = pert_iname, tau,
           idQ = lspci_id_query, idT = pert_id ) %>%
    mutate_at( "drugQ", recode, !!!dmap )

## Compose tau profiles for each drug
V <- R %>% group_by( Drug=drugQ ) %>%
    arrange(idT) %>% summarize( Vals = list(set_names(tau, idT)) )

## Ensure that all tau values are in the same order
stopifnot( map( V$Vals, names ) %>% map_lgl(identical, .[[1]]) %>% all )

## Compute pair-wise similarity for all query drugs
DST <- crossing(rename_all(V, str_c, "1"),
                rename_all(V, str_c, "2")) %>%
    mutate( Similarity = map2_dbl(Vals1, Vals2, cor,
                            method="pearson", use="complete.obs") )

## Perform hierarchical clustering with optimal leaf reordering
DM <- DST %>% select( Drug1, Drug2, Similarity ) %>%
    spread( Drug2, Similarity ) %>% as.data.frame %>%
    column_to_rownames("Drug1") %>% dist
lvl <- hclust(DM) %>% reorder(DM) %>% dendextend::order.hclust() %>% labels(DM)[.]

## Fix order of rows and column based on clustering results
DST <- DST %>% mutate(Drug1 = factor(Drug1, lvl),
                      Drug2 = factor(Drug2, rev(lvl)))

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
ebl <- element_blank
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Plot the similarity
ggplot( DST, aes(Drug1, Drug2, fill=Similarity) ) +
    theme_minimal() +
    geom_tile(color="gray") +
    scale_fill_gradientn( colors=pal, limits=c(-1,1) ) +
    theme(axis.text.x = etxt(10, angle=90, hjust=1, vjust=0.5),
          axis.text.y = etxt(10), axis.title=ebl(),
          legend.text = etxt(12), legend.title=etxt(14)) +
    ggsave( "fig4.pdf", width=8.5, height=7.75 ) +
    ggsave( "fig4.png", width=8.5, height=7.75 )
