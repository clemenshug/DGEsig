library( tidyverse )
library( ggrepel )
library( ggbeeswarm )


pathData <- "~/data/DGEsig"


synapser::synLogin()
syn <- synExtra::synDownloader(pathData, ifcollision="overwrite.local")

pertubation_meta <- syn("syn21547097") %>%
    read_csv()

## Cell line metadata
CM <- tribble(
    ~Short,    ~Tissue,    ~Subtype,
    "A375",   "Skin",     "melanoma",
    "A549",   "Lung",     "adeno",
    "HA1E",   "Kindey",   "normal",
    "HCC515", "Lung",     "adeno",
    "HepG2",  "Liver",    "hcc",
    "HT29",   "Colon",    "adeno",
    "MCF-7",  "Breast",   "lumA",
    "PC3",    "Prostate", "small cell",
    "VCaP",   "Prostate", "AR-V7" ) %>%
    mutate( Label = glue::glue("{Short} ({Tissue}/{Subtype})") )

## Dictionary for cell name standardization
cndict <- c(a375 = "A375", a549 = "A549", ha1e = "HA1E", hcc515 = "HCC515",
            hepg2 = "HepG2", ht29 = "HT29", mcf7 = "MCF-7", pc3 = "PC3", vcap = "VCaP")

## Isolate the relevant data chunk
R0 <- syn("syn21907139") %>% read_rds()

R <- R0 %>% filter( pert_type == "trt_cp", cell_id_query == "mcf7" ) %>%
    select(cellT = cell_id_target, drugT = pert_iname, tau, gene_set,
           idT = lspci_id_target, idQ = lspci_id_query) %>%
    filter( !is.na(idQ), !is.na(idT) ) %>%
    mutate_at( "cellT", recode, !!!cndict ) %>%
    left_join(
        pertubation_meta %>%
            distinct(lspci_id, drugQ = pert_iname),
        by = c("idQ" = "lspci_id")
    )

## MCF7 @ Alvocidib / Flavopiridol / Palbo
S1 <- R %>% filter(idT == idQ) %>%
    inner_join( CM, by = c("cellT"="Short") ) %>%
    mutate_at( c("drugQ","Label"), factor ) %>%
    rename( Tau = tau )

## Split by tau interval
SS1 <- S1 %>% mutate( region = fct_rev(cut(Tau, breaks=c(-100, 95, 100))) ) %>%
    split(.$region)

## Plotting elements
pal <- set_names( ggthemes::few_pal()(7), unique(S1$Tissue) )
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
    theme(axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14))
}

## Additional parameters for each plot
fplot <- function( .df, isTop ) {
    gg <- ggplot(.df, aes(x=drugQ, y=Tau, color=Tissue, group=drugQ)) +
        theme_bw() + theme_bold() +
        geom_beeswarm(cex=3, beeswarmArgs=list(side=-1)) +
        geom_vline(xintercept=c(1.7,2.7), color="lightgray") +
        scale_color_manual( values=pal, drop=FALSE ) +
        scale_x_discrete( drop=FALSE, expand=expansion(add=c(0.3,0.6)),
                         name="3' DGE Query Signature" ) +
        geom_text_repel(aes(label=cellT), show.legend=FALSE,
                        fontface="bold", nudge_x = 0.3,
                        direction="y", segment.color=NA ) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank())

    if( isTop == TRUE )
        gg <- gg + theme(axis.title.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.text.x = element_blank(),
                         legend.position="top")
    else
        gg <- gg + theme(axis.text.x = element_text(hjust=0.3)) +
            guides( color=FALSE )
    gg
}

ggs <- map2( SS1, c(TRUE,FALSE), fplot )
ggcomp <- egg::ggarrange( plots=ggs, ncol=1, heights=c(0.5,0.5), draw=FALSE )
ggsave("fig2.pdf", ggcomp, width=6, height=9)
ggsave("fig2.png", ggcomp, width=6, height=9)


# Comparing
clue_results_l1000 <- syn("syn21907143") %>%
    read_rds() %>%
    filter(result_type == "pert", score_level == "summary") %>%
    chuck("data", 1) %>%
    left_join(
        pertubation_meta %>%
            distinct(lspci_id, drugQ = pert_iname),
        by = c("lspci_id_query" = "lspci_id")
    )

clue_results_dge_old <- synapser::synGet("syn21907139", version = 1) %>%
    chuck("path") %>%
    read_rds() %>%
    filter(result_type == "pert", score_level == "summary") %>%
    chuck("data", 1)

clue_results_alvo <- clue_results_dge_old %>%
    filter( pert_type == "trt_cp", grepl("MCF7", gene_set) ) %>%
    select(drugT = pert_iname, tau, gene_set,
           idT = lspci_id_target, idQ = lspci_id_query) %>%
    filter( !is.na(idQ), !is.na(idT) ) %>%
    mutate(drugQ = str_split(gene_set, "_") %>%
               map_chr(pluck, 5) %>% str_to_lower,
           gene_set = NULL) %>%
    mutate_at( c("drugQ"), factor ) %>%
    rename( Tau = tau ) %>%
    mutate_at(c("drugT","drugQ"), recode,
              alvocidib    = "alvocidib (DGE R1)",
              flavopiridol = "alvocidib (DGE R2)")

corr <- bind_rows(
    clue_results_l1000 %>%
        filter(drugQ == "alvocidib", query_type == "aggregated", z_score_cutoff == 0.7, pert_type == "trt_cp") %>%
        select(drugT = pert_iname, Tau = tau, drugQ) %>%
        group_by(drugT, drugQ) %>%
        summarize(Tau = mean(Tau)) %>%
        ungroup() %>%
        mutate(drugQ = recode(drugQ, alvocidib = "alvocidib (CMap)")),
    clue_results_alvo %>%
        select(drugT, drugQ, Tau) %>%
        group_by(drugT, drugQ) %>%
        summarize(Tau = mean(Tau)) %>%
        ungroup() %>%
        filter(drugQ %in% c("alvocidib (DGE R1)", "alvocidib (DGE R2)"))
) %>%
    spread(drugT, Tau) %>%
    select_if(~!any(is.na(.x)))

corr_mat <- cor(
    corr %>%
        column_to_rownames("drugQ") %>%
        as.matrix() %>%
        t()
)

corr_mat_plot <- pheatmap::pheatmap(
    corr_mat,
    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                        "RdBu")))(100),
    breaks = seq(from = -1, to = 1, length.out = 101),
    scale = "none"
)

ggplot2::ggsave(
    "fig3_b.png",
    corr_mat_plot,
    width = 4, height = 4
)



