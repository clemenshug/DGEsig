library( tidyverse )

pathData <- "~/data/DGEsig"

myhmap <- function(R)
{
    ## Compose a matrix for plotting and perform row clustering
    X <- R %>% spread( drugQ, tau ) %>% as.data.frame() %>% column_to_rownames("drugT")
    h <- hclust(dist(X))

    ## Resize the viewport and plot the heatmap
    vp <- grid::viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))
    setHook("grid.newpage", function() grid::pushViewport(vp), action="prepend")
    pheatmap::pheatmap(X, cluster_rows=h, cluster_cols=h, treeheight_col=0,
                       main="Querying with DGE signatures")
    setHook("grid.newpage", NULL, "replace")

    ## Add axis labels
    grid::grid.text("Query", y=-0.07, gp=grid::gpar(fontsize=16))
    grid::grid.text("Target", x=-0.07, rot=90, gp=grid::gpar(fontsize=16))
}

fig1a <- function()
{
    ## Extracts drug name from gene_set annotation
    fdrug <- function(.x)
        str_sub(.x, 16) %>% str_split("_") %>% map_chr(pluck, 1)
        
    ## Load all relevant results and identify drugs in common
    RDGE <- file.path(pathData, "clue_results_dge.rds") %>% read_rds() %>%
        filter( result_type == "pert", score_level == "summary" ) %>%
        pluck( "data", 1 ) %>%
        filter(lspci_id_target %in% lspci_id_query,
               lspci_id_query %in% lspci_id_target,
               !is.na(lspci_id_target), !is.na(lspci_id_query),
               pert_id != "BRD-K19687926" ) %>%
        mutate( drugQ = fdrug(gene_set) ) %>%
        select( -id, -source, -pert_id, -pert_type, -gene_set, -dataset_query ) %>%
        rename( lspciQ = lspci_id_query, lspciT = lspci_id_target, drugT = pert_iname )

    ## Extract the drug name <-> ID map
    dmap <- RDGE %>% filter( !duplicated(drugT) ) %>% with( set_names(drugT, lspciT) )

    ## Aggregate across multiple entries to compute master similarity score
    SDGE <- RDGE %>% group_by( lspciQ, lspciT ) %>% summarize_at( "tau", max ) %>% ungroup %>%
        mutate_at( c("lspciQ", "lspciT"), as.character ) %>%
        transmute( drugQ = dmap[lspciQ], drugT = dmap[lspciT], tau )

    pdf( "fig1-dge.pdf", width=9, height=7.5 )
    myhmap( SDGE )
    dev.off()
}
