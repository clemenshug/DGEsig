library( tidyverse )
library( seriation )   # For optimal leaf reordering
library( cowplot )

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
pal  <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
##palj <- gray.colors(7, start=1, end=0)
palj <- pal[4:7]
ebl <- element_blank
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Plots a similarity matrix
simplot <- function( .df )
{
    ggplot( .df, aes(DrugID1, DrugID2, fill=Similarity) ) +
        theme_minimal() +
        geom_tile(color="gray") +
        geom_rect(aes(xmin=8.5, xmax=13.5, ymin=ndg-8.5+1, ymax=ndg-13.5+1),
                  fill=NA, color="black", size=1) +
        theme(axis.text = ebl(),
              axis.title=ebl(),
              legend.text = etxt(12), legend.title=etxt(14),
              plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
}

## Plots a zoom panel
zoomplot <- function( .df )
{
    ggplot( .df, aes(DrugID1, DrugID2, fill=Similarity) ) +
        theme_minimal() + geom_tile(color="gray") +
        scale_y_discrete( labels = function(x) dnm[x] ) +
        theme(axis.title=ebl(), axis.text.x=ebl(), axis.ticks.x=ebl(),
              axis.text.y=etxt(12), plot.background=element_rect(color="black", size=2),
              panel.grid.minor=ebl(), panel.grid.major=ebl())
}

## Compose similarity matrices
XTau <- simorder( SM, TauSim )
XJcrd <- simorder( SM, JcrdSim ) %>%
    mutate(across(DrugID1, fct_relevel, levels(XTau$DrugID1)),
           across(DrugID2, fct_relevel, levels(XTau$DrugID2)))
ndg <- length(unique(SM$DrugID2))

## Plot zoom facets
dnm <- set_names(c("bortezomib", "mg132", "fedratinib", "staurosp. agl.", "nilotinib"),
                 c("57736", "36292", "97896", "14772", "100531"))
ZTau <- XTau %>% filter( DrugID1 %in% names(dnm), DrugID2 %in% names(dnm) )
ZJcrd <- XJcrd %>% filter( DrugID1 %in% names(dnm), DrugID2 %in% names(dnm) )

## Plot similarity matrices
ggjcrd <- simplot(XJcrd) + scale_fill_gradientn( colors=palj, limits=c(0,1), name="Jaccard\nSimilarity" )
ggtau  <- simplot(XTau)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Pearson\nCorrelation" )

## Plot zoom facets
gzjcrd <- zoomplot(ZJcrd) + scale_fill_gradientn( colors=palj, limits=c(0,1), guide=FALSE )
gztau  <- zoomplot(ZTau)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), guide=FALSE )

## Put everything together
gg <- egg::ggarrange( plots=list(ggjcrd, ggtau), ncol=2,
                     labels=c(" A"," B"), padding=unit(2,"line"),
                     label.args = list(gp = grid::gpar(font = 4, cex = 4)) ) %>%
    ggdraw() %>%
    + draw_plot( gztau, .68, .7, .17, .25 ) %>%
    + draw_plot( gzjcrd, .18, .7, .17, .25 )
ggsave( "fig6.pdf", gg, width=17, height=7 )
ggsave( "fig6.png", gg, width=17, height=7 )
