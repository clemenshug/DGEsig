library( tidyverse )

pathData <- "~/data/DGEsig"

## Isolate the relevant data chunk
R0 <- file.path(pathData, "clue_results_dge.rds") %>% read_rds() %>%
    pluck( "data", 4 )

## Perform additional wrangling
source("dmap.R")
R <- R0 %>% filter(dataset_query == "sr_repurposing",
                   !is.na(lspci_id_query),
                   !is.na(lspci_id_target),
                   pert_type == "trt_cp" ) %>%
    mutate( drugQ = map_chr(str_split(gene_set, "_"), pluck, 3) ) %>%
    select( drugQ, drugT = pert_iname, tau,
           idQ = lspci_id_query, idT = lspci_id_target ) %>%
    mutate_at( "drugQ", recode, !!!dmap )

