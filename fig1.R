library( tidyverse )

pathData <- "~/data/DGEsig"

fig1b <- function()
{
    ## Load all results
    R <- file.path(pathData, "clue_results_combined.rds") %>% read_rds() %>%
        filter( result_type == "pert", score_level == "summary" ) %>%
        pluck( "data", 1 ) %>%
        rename(idQ = lspci_id_query, idT = lspci_id_target, drugT = pert_iname)

    ## Identify the common set of drugs between DGE-query, L1000-query and targets
    qcom <- R %>% group_by( source ) %>% summarize_at( "idQ", list ) %>%
        with( lift(intersect)(idQ) )
    dmap <- R %>% select( idT, drugT ) %>% filter( idT %in% qcom ) %>%
        distinct() %>% with( set_names(drugT,idT) )

    ## Isolate the appropriate slice of data
    ## Aggregate across multiple entries to compute master similarity score
    R2 <- R %>% filter(idT %in% names(dmap), idQ %in% names(dmap)) %>%
        select( idQ, idT, tau, source ) %>% group_by( idQ, idT, source ) %>%
        summarize_at( "tau", ~.x[ which.max(abs(.x)) ] ) %>% ungroup() %>%
        mutate_at( c("idQ", "idT"), as.character ) %>%
        mutate(drugT = factor(dmap[idT]),
               drugQ = factor(dmap[idQ],
                              levels=rev(levels(drugT))),
               self = ifelse(idQ == idT, "self", "other"))

    ## colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

    pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
    
    ## Test plot
    X <- R2 %>% filter(source == "dge")
    ggplot( X, aes(x=drugT, y=drugQ, fill=tau) ) +
        theme_minimal() + geom_tile(color="black") +
        geom_tile(data=filter(X, idQ==idT), color="black", size=1) +
        scale_fill_gradientn( colors=pal, name="Tau" ) +
        xlab( "Target" ) + ylab( "Query" ) +
        theme( axis.text.x = element_text(angle=90, hjust=1, vjust=0.5) )
}
