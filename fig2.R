library( tidyverse )
library( ggrepel )

pathData <- "~/data/DGEsig"

## Cell line metadata
CM <- tribble(
    ~short,    ~tissue,    ~subtype,
    "A375",   "skin",     "melanoma",
    "A549",   "lung",     "adeno",
    "HA1E",   "kindey",   "normal",
    "HCC515", "lung",     "adeno",
    "HepG2",  "liver",    "hcc",
    "HT29",   "colon",    "adeno",
    "MCF-7",  "breast",   "lumA",
    "PC3",    "prostate", "small cell",
    "VCaP",   "prostate", "AR-V7" ) %>%
    mutate( Cell = glue::glue("{short} ({tissue}/{subtype})") )

## Dictionary for cell name standardization
cndict <- c(a375 = "A375", a549 = "A549", ha1e = "HA1E", hcc515 = "HCC515",
            hepg2 = "HepG2", ht29 = "HT29", mcf7 = "MCF-7", pc3 = "PC3", vcap = "VCaP")

## Isolate the relevant data chunk
R0 <- file.path(pathData, "clue_results_dge.rds") %>% read_rds() %>%
    pluck( "data", 3 )

R <- R0 %>% filter( pert_type == "trt_cp", grepl("MCF7", gene_set) ) %>%
    select(cellT = cell_id_target, drugT = pert_iname, tau, gene_set,
           idT = lspci_id_target, idQ = lspci_id_query) %>%
    filter( !is.na(idQ), !is.na(idT) ) %>%
    mutate_at( "cellT", recode, !!!cndict )

## MCF7 @ Alvocidib / Flavopiridol / Palbo
S1 <- R %>% filter(idT == idQ) %>%
    mutate(drugQ = str_split(gene_set, "_") %>%
               map_chr(pluck, 5) %>% str_to_lower,
           gene_set = NULL) %>%
    inner_join( CM, by = c("cellT"="short") ) %>%
    mutate_at( c("drugQ","Cell"), factor ) %>%
    rename( Tau = tau )

## Split by tau interval
SS1 <- S1 %>% mutate( region = fct_rev(cut(Tau, breaks=c(-100, -1, 25, 50, 95, 100))) ) %>%
    split(.$region)

## Plotting elements
pal <- set_names( c("black", ggthemes::few_pal()(8)), unique(S1$Cell) )
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
    theme(axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14))
}

## Additional parameters for each plot
pparam <- list(list(ybr = waiver(),    keepx=FALSE, keepy=FALSE, minor=TRUE),
               list(ybr = c(80,82,84), keepx=FALSE, keepy=TRUE,  minor=FALSE),
               list(ybr = c(34,39,44), keepx=FALSE, keepy=FALSE, minor=TRUE),
               list(ybr = 0,           keepx=FALSE, keepy=FALSE, minor=TRUE),
               list(ybr = -80.8,       keepx=TRUE,  keepy=FALSE, minor=TRUE)
               )

fplot <- function( .df, lparam ) {
    gg <- ggplot(.df, aes(x=drugQ, y=Tau, color=Cell)) + 
        geom_point() + theme_bw() + theme_bold() +
        scale_color_manual( values=pal, drop=FALSE ) +
        scale_x_discrete( drop=FALSE, name="Query Signature" ) +
        scale_y_continuous( breaks = lparam$ybr ) +
        geom_text_repel(aes(label=cellT), show.legend=FALSE, fontface="bold")

    if( lparam$keepx == FALSE )
        gg <- gg + theme(axis.title.x = element_blank(),
                         axis.text.x = element_blank())

    if( lparam$keepy == FALSE )
        gg <- gg + theme(axis.title.y = element_blank()) +
            guides( color=FALSE )

    if( lparam$minor == FALSE )
        gg <- gg + theme(panel.grid.minor = element_blank())
        
    gg
}

ggs <- map2( SS1, pparam, fplot )
ggcomp <- egg::ggarrange( plots=ggs, ncol=1, heights=c(0.5,0.2,0.2,0.1,0.1) )
ggsave("fig2.pdf", ggcomp, width=8, height=9)
ggsave("fig2.png", ggcomp, width=8, height=9)

