library( tidyverse )

pathData <- "~/data/DGEsig"

fig1 <- function()
{
    ## Load all relevant results
    RDGE <- file.path(pathData, "clue_results_dge.rds") %>% read_rds()

    ## Summary-level slice for perturbation queries
    R1 <- RDGE %>% filter( result_type == "pert", score_level == "summary" ) %>%
        pluck( "data", 1 ) %>% select(-id, -source)

    ## Extracts drug name from gene name
    fdrug <- function(.x)
        str_sub(.x, 16) %>% str_split("_") %>% map_chr(pluck, 1)
        
    ## Identify drugs in common
    R2 <- R1 %>% filter(lspci_id_target %in% lspci_id_query,
                        lspci_id_query %in% lspci_id_target,
                        !is.na(lspci_id_target), !is.na(lspci_id_query)) %>%
        filter(pert_id != "BRD-K19687926") %>% 
        mutate( drugQ = fdrug(gene_set) ) %>%
        select( -pert_id, -pert_type, -gene_set, -dataset_query ) %>%
        rename( lspciQ = lspci_id_query, lspciT = lspci_id_target, drugT = pert_iname )

    ## Extract the drug name <-> ID map
    dmap <- R2 %>% filter( !duplicated(drugT) ) %>% with( set_names(drugT, lspciT) )

    ## Aggregate across multiple entries
    R3 <- R2 %>% group_by( lspciQ, lspciT ) %>% summarize_at( "tau", max ) %>% ungroup
}
