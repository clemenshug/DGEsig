library( tidyverse )
library( seriation )   # For optimal leaf reordering

pathData <- "~/data/DGEsig"

source("dmap.R")    # Map for standardizing drug names
                    # (not used in the current code version)

## Load gene sets associated with Steve's drugs
ct <- cols(lspci_id = col_integer(), query_type = col_character(),
           cutoff   = col_double(),  cell_id    = col_character())
GS0 <- file.path(pathData, "all_gene_sets.csv.gz") %>% read_csv( col_types=ct )
GS <- GS0 %>% filter( grepl("rencell", gene_set), !is.na(symbol) ) %>%
    group_by( DrugID=lspci_id ) %>% summarize( Set = list(symbol) )

## Isolate the relevant data chunk
R0 <- file.path(pathData, "clue_results_dge.rds") %>% read_rds() %>%
    pluck( "data", 4 )

## Perform additional wrangling
R <- R0 %>% filter(grepl("rencell", gene_set),
                   !is.na(lspci_id_query),
                   !is.na(lspci_id_target),    
                   pert_type == "trt_cp" ) %>%
    select( idQ = lspci_id_query, idT = pert_id, tau )

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

## Perform hierarchical clustering with optimal leaf reordering
## Fixes row order based on the clustering results
## .df - data frame, .sim - similarity column
simorder <- function( .df, .sim )
{
    ## Compose the distance matrix
    DM <- .df %>% select( DrugID1, DrugID2, {{.sim}} ) %>%
        spread( DrugID2, {{.sim}} ) %>% as.data.frame() %>%
        column_to_rownames("DrugID1") %>% dist()
    
    ## Perform hierarchical clustering with optimal leaf reordering
    lvl <- hclust(DM) %>% reorder(DM) %>%
        dendextend::order.hclust() %>% labels(DM)[.]
    
    ## Fix order of rows and column based on clustering results
    .df %>% transmute(DrugID1 = factor(DrugID1, lvl),
                      DrugID2 = factor(DrugID2, rev(lvl)),
                      Similarity = {{.sim}})
}

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
ebl <- element_blank
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Plot a similarity matrix
simplot <- function( .df )
{
    ggplot( .df, aes(DrugID1, DrugID2, fill=Similarity) ) +
        theme_minimal() +
        geom_tile(color="gray") +
        scale_fill_gradientn( colors=pal, limits=c(-1,1) ) +
        theme(axis.text.x = etxt(10, angle=90, hjust=1, vjust=0.5),
              axis.text.y = etxt(10), axis.title=ebl(),
              legend.text = etxt(12), legend.title=etxt(14))
}

ggtau  <- simorder( SM, TauSim )  %>% simplot() %>%
    + ggsave( "fig6a.pdf", width=8.5, height=7.75 ) %>%
    + ggsave( "fig6a.png", width=8.5, height=7.75 )
ggjcrd <- simorder( SM, JcrdSim ) %>% simplot() %>%
    + ggsave( "fig6b.pdf", width=8.5, height=7.75 ) %>%
    + ggsave( "fig6b.png", width=8.5, height=7.75 )
