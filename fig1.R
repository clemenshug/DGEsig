library( tidyverse )
library( seriation )   # For optimal leaf reordering

pathData <- "~/data/DGEsig"

synapser::synLogin()
syn <- synExtra::synDownloader(pathData, ifcollision="overwrite.local")


## Load all results
R <- file.path(pathData, "clue_results_combined.rds") %>% read_rds() %>%
    filter( result_type == "pert", score_level == "summary" ) %>%
    pluck( "data", 1 ) %>%
    rename(idQ = lspci_id_query, idT = lspci_id_target, drugT = pert_iname)

condition_conc_vars <- c("cells", "drug_id", "lspci_id", "stim", "stim_conc", "time")

M <- syn("syn22000707") %>%
    read_rds() %>%
    unnest(meta) %>%
    mutate(
        gene_set = exec(
            paste,
            !!!as.list(.)[condition_conc_vars],
            sep = "_"
        ) %>%
            str_replace_all("\\s", "_") %>%
            str_replace_all("[^\\w]", "")
    )

M_cell_info <- M %>%
    filter(dataset != "fp_transdiff") %>%
    distinct(gene_set, lspci_id, cells)

## Identify the common set of drugs between DGE-query, L1000-query and targets
qcom <- R %>% group_by( source ) %>% summarize_at( "idQ", list ) %>%
    with( lift(intersect)(idQ) )
dmap <- R %>% select( idT, drugT ) %>% filter( idT %in% qcom ) %>%
    distinct() %>% with( set_names(drugT,idT) )

## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2 <- R %>% filter(idT %in% names(dmap), idQ %in% names(dmap)) %>%
    select( idQ, idT, tau, source, z_score_cutoff, query_type, cell_id_query ) %>%
    group_by( idQ, idT, source, z_score_cutoff, query_type ) %>%
    summarize(
        "tau" = tau[ which.max(abs(tau)) ],
        cell_id_query = list(unique(cell_id_query))
    ) %>%
    ungroup() %>%
    mutate_at( c("idQ", "idT"), as.character ) %>%
    mutate_at( "source", toupper ) %>%
    mutate( drugT = dmap[idT], drugQ = dmap[idQ] )
    # nest_join(
    #     transmute(M_cell_info, lspci_id = as.character(lspci_id), cells) %>%
    #         distinct(),
    #     by = c("idQ" = "lspci_id"),
    #     name = "cells"
    # )

## Perform hierarchical clustering on drugT profiles (columns in the final plot)
## Use the DGE slice because it is cleaner and less saturated
DM <- R2 %>% filter(query_type == "aggregated", z_score_cutoff == 0.7) %>% select( source, drugT, drugQ, tau ) %>%
    spread( drugT, tau ) %>% select(-source, -drugQ) %>% as.matrix() %>% t() %>% dist
lvl <- hclust(DM) %>% reorder(DM) %>%  dendextend::order.hclust() %>% labels(DM)[.]

## Fix the order via factor levels
R2 <- R2 %>% mutate(drugT = factor(drugT, lvl),
                    drugQ = factor(drugQ, rev(lvl)))

# Complete missing observations at z-scores that yielded insufficient
# genes for Clue with NA
R2_completed <- bind_rows(
    R2 %>%
        filter(source == "DGE"),
    R2 %>%
        filter(source == "L1000") %>%
        complete(nesting(idQ, idT, source, drugT, drugQ), z_score_cutoff, query_type)
)

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
    theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
          axis.text.y = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          axis.ticks = element_blank())
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

composite_plot <- function(X) {
    ## DGE plot
    gg1 <- fplot( filter(X, source == "DGE") ) +
        scale_x_discrete(position = "top") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(hjust=0, vjust=0.5)) +
        ylab( "DGE Query" )

    ## L1000 plot
    gg2 <- fplot( filter(X, source == "L1000") ) +
        ylab( "L1000 Query" ) +
        scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
        theme(legend.position = "bottom")

    cell_plot <- X %>%
        distinct(drugQ, cell_id_query) %>%
        # unchop(cell_id_query) %>%
        mutate(
            cell_id_query = map_chr(
                cell_id_query,
                ~if (is.null(.x) || length(.x) == 1) .x %||% "" else "multiple"
            )
        ) %>%
        drop_na() %>%
        ggplot(aes(drugQ, fill = cell_id_query)) +
            geom_bar() +
            coord_flip() +
            theme_minimal() + theme_bold() +
            scale_fill_brewer(palette = "Set2", name = "Query cell line") +
            theme(
                # axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                # axis.title.y = element_blank(), axis.title.x = element_blank(),
                axis.title = element_blank(), axis.ticks = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank()
            )

    ## Summary plot
    S <- X %>% filter(idQ == idT) %>%
        mutate_at("source", factor, levels=c("L1000","DGE")) %>%
        mutate_at("source", fct_recode, `Self (DGE)`="DGE", `Self (L1000)`="L1000")
    ggs <- ggplot( S, aes(x=drugT, y=source, fill=tau) ) +
        theme_minimal() + theme_bold() +
        geom_tile(color="black") + ylab("") +
        scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank())

    ## Create the composite plot
    egg::ggarrange( gg1, cell_plot, ggs, egg::.dummy_ggplot, gg2, egg::.dummy_ggplot, heights=c(6.5,0.5,7.5), widths = c(6.5, 0.3), draw=FALSE )
}

ggcomp <- R2_completed %>%
    filter( source == "L1000" ) %>%
    group_nest( z_score_cutoff, query_type ) %>%
    mutate(
        data = map(
            data,
            bind_rows,
            filter( R2_completed, source == "DGE")
        ) %>%
            map(composite_plot)
    )

# ggcomp <- gridExtra::arrangeGrob(
#     grobs = map2(
#         ggcomp$data, ggcomp$z_score_cutoff,
#         ~gridExtra::arrangeGrob(.x, top = paste0("z-cutoff ", .y * 100, "%"))
#     ),
#     # set_names( ggcomp$data, ggcomp$z_score_cutoff ),
#     nrow = 1
# )

pwalk(
    ggcomp,
    function(z_score_cutoff, query_type, data, ...) {
        walk( paste0("fig1_", z_score_cutoff, "_", query_type, c(".png", ".pdf")), ggsave, data, width=8, height=13 )
    }
)

walk( c("fig1.pdf", "fig1.png"), ggsave, ggcomp2, width=28, height=12 )

## Query by cell-line

R3 <- R %>% filter(idT %in% names(dmap), idQ %in% names(dmap)) %>%
    select( idQ, idT, tau, source, z_score_cutoff, query_type, cell_id_query ) %>%
    group_by( idQ, idT, source, z_score_cutoff, query_type, cell_id_query ) %>%
    summarize_at( "tau", ~.x[ which.max(abs(.x)) ] ) %>% ungroup() %>%
    mutate_at( c("idQ", "idT"), as.character ) %>%
    mutate_at( "source", toupper ) %>%
    mutate( drugT = dmap[idT], drugQ = dmap[idQ] )

R_torin <- R3 %>%
    filter(query_type == "per_cell_line") %>%
    mutate(drugQ = cell_id_query) %>%
    bind_rows(
        R2 %>%
            filter(query_type == "aggregated", source == "L1000", drugQ == "torin-1") %>%
            mutate(drugQ = "aggregated")
    ) %>%
    bind_rows(
        crossing(
            drugQ = "",
            drugT = unique(.[["drugT"]]),
            z_score_cutoff = unique(.[["z_score_cutoff"]])
        ) %>%
            mutate(tau = 0)
    )

plot_torin <- R_torin %>%
    group_nest(z_score_cutoff) %>%
    mutate(
        data = map(
            data,
            ~.x %>%
                mutate(
                    drugQ = fct_relevel(drugQ, "aggregated", ""),
                    border = case_when(
                        drugQ == "" ~ 0,
                        idQ == idT ~ 1,
                        TRUE ~ 0.4
                    ),
                    color = if_else(drugQ == "", NA_character_, "black")
                ) %>%
                ggplot(aes(x=drugT, y=drugQ, fill=tau, size=border) ) +
                scale_size_identity() +
                theme_minimal() + theme_bold() +
                geom_tile(aes(color = color)) +
                scale_color_identity() +
                # geom_tile(data=filter(X, idQ==idT), color="black", size=1) +
                scale_fill_gradientn( colors=pal, limits=c(-100,100) ) +
                xlab( "Clue Target" ) +
                ylab("Query cell line")
        )
    )

pwalk(
    plot_torin,
    function(z_score_cutoff, data, ...) {
        walk(
            paste0("figSx_torin-1_per_cell_line_", z_score_cutoff, c(".pdf", ".png")),
            ggsave,
            data,
            width = 7, height = 6
        )
    }
)

## Target by cell-line

R_all <- file.path(pathData, "clue_results_combined.rds") %>% read_rds()

R4 <- R_all %>%
    filter( result_type == "pert", score_level == "cell" ) %>%
    pluck( "data", 1 ) %>%
    rename(idQ = lspci_id_query, idT = lspci_id_target, drugT = pert_iname)%>%
    filter(idT %in% names(dmap), idQ %in% names(dmap)) %>%
    select( idQ, idT, tau, source, z_score_cutoff, query_type, cell_id_query, cell_id_target ) %>%
    group_by( idQ, idT, source, z_score_cutoff, query_type, cell_id_target ) %>%
    summarize_at( "tau", ~.x[ which.max(abs(.x)) ] ) %>% ungroup() %>%
    mutate_at( c("idQ", "idT"), as.character ) %>%
    mutate_at( "source", toupper ) %>%
    mutate( drugT = dmap[idT], drugQ = dmap[idQ] )

R4_torin <- R4 %>%
    filter(query_type == "aggregated", idQ == 101253) %>%
    mutate(drugQ = str_to_upper(cell_id_target)) %>%
    bind_rows(
        R2 %>%
            filter(query_type == "aggregated", drugQ == "torin-1") %>%
            mutate(drugQ = "aggregated")
    ) %>%
    bind_rows(
        crossing(
            drugQ = "",
            drugT = unique(.[["drugT"]]),
            z_score_cutoff = unique(.[["z_score_cutoff"]]),
            source = unique(.[["source"]])
        ) %>%
            mutate(tau = 0)
    )

plot_torin2 <- R4_torin %>%
    group_nest(source, z_score_cutoff) %>%
    mutate(
        data = map(
            data,
            ~.x %>%
                mutate(
                    drugQ = fct_relevel(drugQ, "aggregated", ""),
                    border = case_when(
                        drugQ == "" ~ 0,
                        idQ == idT ~ 1,
                        TRUE ~ 0.4
                    ),
                    color = if_else(drugQ == "", NA_character_, "black")
                ) %>%
                ggplot(aes(x=drugT, y=drugQ, fill=tau, size=border) ) +
                scale_size_identity() +
                theme_minimal() + theme_bold() +
                geom_tile(aes(color = color)) +
                scale_color_identity() +
                # geom_tile(data=filter(X, idQ==idT), color="black", size=1) +
                scale_fill_gradientn( colors=pal, limits=c(-100,100) ) +
                xlab( "Clue Target" ) +
                ylab("Target cell line")
        )
    )

pwalk(
    plot_torin2,
    function(z_score_cutoff, data, source, ...) {
        walk(
            paste0("figSx_torin-1_target_cell_line_", source, "_", z_score_cutoff, c(".pdf", ".png")),
            ggsave,
            data,
            width = 7, height = 3.5
        )
    }
)
