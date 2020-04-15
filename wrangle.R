library( tidyverse )

synapser::synLogin()
dloc <- "~/data/DGEsig"
syn <- synExtra::synDownloader(dloc, ifcollision="overwrite.local")

compound_name_map <- synapser::synGet("syn21586544", downloadLocation=dloc,
                                      ifcollision="overwrite.local", version = 1)

cmap_gene_meta    <- syn("syn21547102")
clue_res_dge      <- syn("syn21907139")
clue_res_l1000    <- syn("syn21907143")
clue_res_combined <- syn("syn21907166")
diff_exp_by_conc  <- syn("syn21559856")
pertubation_meta  <- syn("syn21547097")
signature_meta    <- syn("syn21547101")
