library( tidyverse )

pathData <- "~/data/DGEsig"

fig1b <- function() {
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
        mutate_at( "source", toupper ) %>%
        mutate(drugT = factor(dmap[idT]),
               drugQ = factor(dmap[idQ],
                              levels=rev(levels(drugT))))

    ## Plotting elements
    pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme_bold <- function() {
        theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
              axis.text.y = etxt(12), axis.title = etxt(14),
              legend.text = etxt(12), legend.title = etxt(14))
    }

    ## Plotting a heatmap of clue hits
    fplot <- function(X) {
        ggplot( X, aes(x=drugT, y=drugQ, fill=tau) ) +
            theme_minimal() + theme_bold() +
            geom_tile(color="black") +
            geom_tile(data=filter(X, idQ==idT), color="black", size=1) +
            scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
            xlab( "Clue Target" )
    }
    
    ## DGE plot
    gg1 <- fplot( filter(R2, source == "DGE") ) +
        scale_x_discrete(position = "top") +
        theme(axis.title.x = element_blank()) +
        ylab( "DGE Query" )

    ## L1000 plot
    gg2 <- fplot( filter(R2, source == "L1000") ) +
        ylab( "L1000 Query" )

    ## Summary plot
    S <- R2 %>% filter(idQ == idT) %>%
        mutate_at("source", factor, levels=c("L1000","DGE")) %>%
        mutate_at("source", fct_recode, `Self (DGE)`="DGE", `Self (L1000)`="L1000")
    ggs <- ggplot( S, aes(x=drugT, y=source, fill=tau) ) +
        theme_minimal() + theme_bold() +
        geom_tile(color="black") + ylab("") +
        scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank())

    ## Create the composite plot
    ggcomp <- egg::ggarrange( gg1, ggs, gg2, heights=c(6.5,1,6.5) )
    ggsave( "fig1.pdf", ggcomp, width=7, height=12 )
}
