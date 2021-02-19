library( tidyverse )
library( seriation )   # For optimal leaf reordering
library(synapser)

pathData <- "~/data/DGEsig"

synapser::synLogin()
syn <- synExtra::synDownloader(pathData, ifcollision="overwrite.local")

rdbu_colormap <- RColorBrewer::brewer.pal(5, "RdBu")

cmap_nonlinear_colormap <- list(
    colors = c(
        rdbu_colormap[1:3],
        rdbu_colormap[3:5]
    ) %>%
        rev(),
    values = scales::rescale(c(-100, -90, -80, 80, 90, 100), from = c(-100, 100))
)
# scales::gradient_n_pal(
#     c(
#         rdbu_colormap[1:3],
#         rdbu_colormap[3:5]
#     ),
#     scales::rescale(c(-100, -90, -80, 80, 90, 100), from = c(-100, 100))
# )

## Load all results
R <- synapser::synGet("syn21907166", version = 5, ifcollision = "overwrite.local") %>%
    chuck("path") %>%
    # file.path(pathData, "clue_results_combined.rds") %>%
    read_rds() %>%
    filter( result_type == "pert", score_level == "summary" ) %>%
    pluck( "data", 1 ) %>%
    dplyr::rename(idQ = lspci_id_query, idT = lspci_id_target, drugT = pert_iname)

condition_conc_vars <- c("cells", "drug_id", "lspci_id", "stim", "stim_conc", "time")

M <- synGet("syn22000707", version = 8, ifcollision = "overwrite.local") %>%
    chuck("path") %>%
    # syn("syn22000707") %>%
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
        # scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
        exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, guide = FALSE, limits = c(-100, 100)) +
        xlab( "CMap Target" )
}

composite_plot <- function(X) {
    ## DGE plot
    gg1 <- fplot( filter(X, source == "DGE") ) +
        scale_x_discrete(position = "top") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(hjust=0, vjust=0.5)) +
        ylab( "3' DGE Query" )

    ## L1000 plot
    gg2 <- fplot( filter(X, source == "L1000") ) +
        ylab( "L1000 Query" ) +
        # scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
        exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, name = "Tau", limits = c(-100, 100)) +
        theme(legend.position = "bottom")
#
#     cell_plot <- X %>%
#         distinct(drugQ, cell_id_query) %>%
#         # unchop(cell_id_query) %>%
#         mutate(
#             cell_id_query = map_chr(
#                 cell_id_query,
#                 ~if (is.null(.x) || length(.x) == 1) .x %||% "" else "multiple"
#             )
#         ) %>%
#         drop_na() %>%
#         ggplot(aes(drugQ, fill = cell_id_query)) +
#             geom_bar() +
#             coord_flip() +
#             theme_minimal() + theme_bold() +
#             scale_fill_brewer(palette = "Set2", name = "Query cell line") +
#             theme(
#                 # axis.ticks.y = element_blank(), axis.text.y = element_blank(),
#                 # axis.title.y = element_blank(), axis.title.x = element_blank(),
#                 axis.title = element_blank(), axis.ticks = element_blank(),
#                 axis.text.x = element_blank(), axis.text.y = element_blank()
#             )

    ## Summary plot
    S <- X %>% filter(idQ == idT) %>%
        mutate_at("source", factor, levels=c("L1000","DGE")) %>%
        mutate_at("source", fct_recode, `Self (DGE)`="3' DGE", `Self (L1000)`="L1000")
    ggs <- ggplot( S, aes(x=drugT, y=source, fill=tau) ) +
        theme_minimal() + theme_bold() +
        geom_tile(color="black") + ylab("") +
        # scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
        exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, guide = FALSE, limits = c(-100, 100)) +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank())

    ## Create the composite plot
    egg::ggarrange( gg1, ggs, gg2, heights=c(6.5,0.5,7.5), widths = c(6.5), draw=FALSE )
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
        walk( paste0("fig1_", z_score_cutoff, "_", query_type, "_", "cmap_limits", c(".png", ".pdf")), ggsave, data, width=6.5, height=13 )
    }
)

walk( c("fig1.pdf", "fig1.png"), ggsave, ggcomp2, width=28, height=12 )

# Test if self-similarity is significantly higher than other similarities

wilcox_self_similarities <- R2 %>%
    filter(query_type == "aggregated", is.na(z_score_cutoff) | z_score_cutoff == 0.7) %>%
    group_by(source, idQ, drugQ) %>%
    summarize(
        wilcox = wilcox.test(
            tau[idQ != idT],
            mu = tau[idQ == idT],
            alternative = "less"
        ) %>%
            list(),
        .groups = "drop"
    ) %>%
    mutate(
        p.value = map_dbl(wilcox, "p.value")
    )

library(ggbeeswarm)

self_similarity_beeswarm <- R2 %>%
    filter(query_type == "aggregated", is.na(z_score_cutoff) | z_score_cutoff == 0.7) %>%
    mutate(
        self_similarity = if_else(
            idQ == idT,
            "self_similarity",
            "cross_similarity"
        ) %>%
            fct_relevel("cross_similarity")
    ) %>%
    arrange(self_similarity) %>%
        ggplot(
            aes(tau, drugQ, color = self_similarity)
        ) +
        geom_quasirandom(
            data = ~filter(.x, self_similarity == "cross_similarity"),
            groupOnX = FALSE
        ) +
        geom_point(
            data = ~filter(.x, self_similarity == "self_similarity")
        ) +
        facet_wrap(~source, nrow = 1) +
        scale_color_manual(
            values = c(
                self_similarity = "#FF0000",
                cross_similarity = "#00000088"
            ),
            guide = FALSE
        )

dir.create("self_similarity")
ggsave(
    file.path("self_similarity", "self_similarity_beeswarm.pdf"),
    self_similarity_beeswarm,
    width = 5, height = 10
)

self_similarity_stats <- R2 %>%
    filter(query_type == "aggregated", is.na(z_score_cutoff) | z_score_cutoff == 0.7) %>%
    group_by(source, idQ, drugQ) %>%
    summarize(
        n_greater = sum(
            tau[idQ != idT] > tau[idQ == idT]
        ),
        .groups = "drop"
    )

wilcox_self_similarities_col_plot <- wilcox_self_similarities %>%
    mutate(log_p = -log10(p.value)) %>%
    ggplot(aes(x = log_p, y = drugQ, fill = source)) +
    geom_col() +
    facet_wrap(~source, nrow = 1) +
    guides(fill = FALSE)

ggsave(
    file.path("self_similarity", "wilcox_cols.pdf"),
    wilcox_self_similarities_col_plot,
    width = 4, height = 8
)

wilcox_self_similarities_scatterplot <- wilcox_self_similarities %>%
    mutate(log_p = -log10(p.value)) %>%
    select(idQ, source, p.value) %>%
    spread(source, p.value) %>%
    ggplot(aes(x = DGE, y = L1000)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    coord_equal()



# Plot TAS vector similarity

# Use old version of TAS vectors that matches used lspci_ids
tas <- syn("syn20830942.4") %>%
    read_csv()

library(data.table)

lspci_id_name_map <- R2 %>%
    distinct(lspci_id = idQ, name = drugQ) %>%
    mutate(across(lspci_id, as.double))

tas_used <- tas %>%
    filter(fp_name == "morgan_normal", lspci_id %in% R2$idQ) %>%
    distinct(lspci_id, gene_id = entrez_gene_id, tas) %>%
    setDT()

tas_weighted_jaccard <- function(data_tas, query_id, min_n = 6) {
    query_tas <- data_tas[lspci_id == query_id, .(gene_id, tas)]
    data_tas[
        ,
        .(lspci_id, gene_id, tas)
    ][
        query_tas,
        on = "gene_id",
        nomatch = NULL
    ][
        ,
        mask := tas < 10 | i.tas < 10
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

all_similarity <- tibble(lspci_id_1 = unique(tas_used$lspci_id)) %>%
    mutate(
        data = map(
            lspci_id_1,
            ~tas_weighted_jaccard(tas_used, .x) %>%
                rename(lspci_id_2 = lspci_id)
        )
    ) %>%
    unnest(data) %>%
    left_join(
        lspci_id_name_map %>%
            rename(name_1 = name, lspci_id_1 = lspci_id)
    ) %>%
    left_join(
        lspci_id_name_map %>%
            rename(name_2 = name, lspci_id_2 = lspci_id)
    ) %>%
    mutate(across(starts_with("name"), factor, levels = lvl)) %>%
    mutate(across(name_2, fct_rev))

tas_similarity_plot <- ggplot( all_similarity, aes(x=name_1, y=name_2, fill=tas_similarity) ) +
    theme_minimal() + theme_bold() +
    geom_tile(color="black") +
    geom_tile(data=filter(all_similarity, lspci_id_1==lspci_id_2), color="black", size=1) +
    # scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(0, 1) ) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    # xlab( "CMap Target" ) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(hjust=0, vjust=0.5),
          axis.title.y = element_blank())
    # ylab( "3' DGE Query")

dir.create("tas_similarity")
ggsave(
    file.path("tas_similarity", "tas_similarity_heatmap.pdf"),
    width = 8, height = 6
)

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
