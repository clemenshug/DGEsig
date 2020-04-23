library( tidyverse )

pathData <- "~/data/DGEsig"

R0 <- file.path(pathData, "clue_results_combined.rds") %>% read_rds() %>%
    pluck( "data", 3 )

R <- R0 %>% filter( pert_type == "trt_cp", source == "l1000" ) %>%
    select(cellT = cell_id_target, drugT = pert_iname, drugQ = name_query,
           tau, idT = lspci_id_target, idQ = lspci_id_query) %>%
    filter( !is.na(idT) )

fztau <- function(.df)
    .df %>% summarize( fz = sum(tau==0)/n() ) %>% arrange( desc(fz) )

R %>% group_by( idQ ) %>% fztau()
R %>% group_by( cellT ) %>% fztau()

R %>% group_by( idT ) %>% fztau()
## # A tibble: 2,243 x 2
##       idT    fz
##     <dbl> <dbl>
##  1    469 0.611
##  2   9409 0.590
##  3  82152 0.556
##  4 414384 0.545
##  5  87721 0.542
##  6  66592 0.531
##  7 101244 0.528
##  8  21311 0.524
##  9 149375 0.523
## 10 113909 0.514

SM <- file.path(pathData, "cmap_signature_meta.csv.gz") %>% read_csv() %>%
    select( dataset, idQ=lspci_id ) %>% distinct() %>% na.omit()

Z <- R %>% filter( idT == "469" ) %>% inner_join(SM)

M1 <- read_tsv( "ext/GSE92742_Broad_LINCS_sig_metrics.txt" )
M2 <- read_tsv( "ext/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz" )

MZ <- M %>% filter( pert_iname == "PIK-75" )
