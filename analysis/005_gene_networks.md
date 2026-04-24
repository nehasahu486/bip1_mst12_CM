Magnaporthe oryzae all_DEGs from Bip1 and Mst12 comparisons at 0h, 4h
and 24h - Gene Function Network Analysis
================
Neha Sahu
24 April, 2026

- [Description](#description)
  - [GoHREP](#gohrep)
  - [Load Required Packages](#load-required-packages)
  - [Network layout](#network-layout)
  - [Base Network](#base-network)
  - [Mapping differentially expressed genes on the Base
    Network](#mapping-differentially-expressed-genes-on-the-base-network)
    - [Downreg in Guy11 at 0h (Upregulated in Bip1/Mst12 at
      0h)](#downreg-in-guy11-at-0h-upregulated-in-bip1mst12-at-0h)
    - [Downreg in Guy11 at 4h (Upregulated in Bip1/Mst12 at
      4h)](#downreg-in-guy11-at-4h-upregulated-in-bip1mst12-at-4h)
    - [Downreg in Guy11 at 24h (Upregulated in Bip1/Mst12 at
      24h)](#downreg-in-guy11-at-24h-upregulated-in-bip1mst12-at-24h)
    - [Upreg in Guy11 at 0h (Downregulated in Bip1/Mst12 at
      0h)](#upreg-in-guy11-at-0h-downregulated-in-bip1mst12-at-0h)
    - [Upreg in Guy11 at 4h (Downregulated in Bip1/Mst12 at
      4h)](#upreg-in-guy11-at-4h-downregulated-in-bip1mst12-at-4h)
    - [Upreg in Guy11 at 24h (Downregulated in Bip1/Mst12 at
      24h)](#upreg-in-guy11-at-24h-downregulated-in-bip1mst12-at-24h)

# Description

Camilla wants to have a visual summary of which gene functions are
coming up most among the DEGs and if their expression is co-regulated by
bip1 and mst12 OR is specific.

For this we discussed and planned on creating networks . I will first
create the base network file in GEXF format for import into Gephi, where
it will be visualized with various layouts.

## GoHREP

**Goal:**

To prepare gene function networks from DEGs from this project and see
their temoral profiles during in planta infection, and to identify key
hubs that may coordinate virulence functions over the time course.

**Hypothesis:**  
We think that based on the protein domains, their will be hubs of
different functions of genes expressed during specific stages of plant
infection that work together to suppress host immunity and/or facilitate
pathogen colonization. We also think there will be with a hub of
effectors playing central roles in coordinating virulence.

**Rationale:**  
Network analysis of differentially expressed genes when annotated
broadly into functional categories based on InterPro domains and known
characterized gene info can reveal groups of gene that may act together
over the time course to coordinate pathogenicity or virulence-related
functions.

**Experimental plan:**

1.  Discuss and create a broad fucntional annotations file with Cami -
    using the protein domain information from InterPRo annotations, and
    other known pathways in Magnaporthe. This file will be saved at
    `data/processed/network_base_manual.tsv`

2.  Generate a network file of DEGs based on RNA-seq data across bip1
    and mst12 samples.

3.  Import the gexf file to be further used in Gephi and plot various
    networks based on time points

## Load Required Packages

``` r
library(tidyverse)
```

## Network layout

Base network colored based on functional categories

``` r
dat <- read_tsv(here("data/processed/network_base_manual.tsv")) %>%
  filter(!Function_broad_final %in% c("other_unknown","uncharacterized","metabolic_enzymes")) # remove categories which have more than 300 genes, to avoid skewing the netwrok
unique(dat$Function_broad_final)
```

    ##  [1] "transport"                     "transcription_factor"         
    ##  [3] "binding_scaffold_other"        "kinases"                      
    ##  [5] "translation_ribosome"          "cell_wall_enzymes"            
    ##  [7] "effector"                      "oxidoreductases"              
    ##  [9] "ubiquitin_proteasome"          "cytoskeleton"                 
    ## [11] "dna_rna_protein_binding"       "dna_rna_processing_and_repair"
    ## [13] "signal_transduction"           "phosphatases"                 
    ## [15] "proteases_peptidases"          "nlrs"                         
    ## [17] "vesicular_transport"           "secreted_extracellular"       
    ## [19] "mitochondria_related"          "protein_folding_chaperones"   
    ## [21] "lipid_metabolism"              "cell_division"                
    ## [23] "sln1_pathway"                  "septin_oligomerization"       
    ## [25] "pmk1_pathway"

``` r
dat <- dat %>%
  filter(!is.na(Id),trimws(Id) != "",
         !is.na(Function_broad_final),  trimws(Function_broad_final)  != "") %>%
  mutate(Function_broad_final = trimws(Function_broad_final),
         Id = trimws(Id))

# Assign colors for the various functional categories - Used ChatGPT to generate color-blind friendly palette

palette <- c(
  "binding_scaffold_other"         = "#6BA3A8",
  "transport"                      = "#E89B6F",
  "signal_transduction"            = "#C77D9E",
  "transcription_factor"           = "#7B9ACA",
  "effector"                       = "#D97B7B",
  "dna_rna_protein_binding"        = "#95B57D",
  "cell_wall_enzymes"              = "#C9A76B", 
  "translation_ribosome"           = "#9B85BA",
  "oxidoreductases"                = "#6BADA3", 
  "proteases_peptidases"           = "#D88FA3",
 
  "mitochondria_related"           = "#5A7A8C", 
  "ubiquitin_proteasome"           = "#9C6B87",
  "kinases"                        = "#6D9178",
  "protein_folding_chaperones"     = "#B87D50", 
  "cytoskeleton"                   = "#8C7D5E", 
  "vesicular_transport"            = "#6B8BA8", 
  "phosphatases"                   = "#A66D6D", 
  "cell_division"                  = "#5E7E7A",
  "lipid_metabolism"               = "#9A7B97", 
  "secreted_extracellular"         = "#7A9A6E", 

  "sln1_pathway"                   = "#B8D4A8", 
  "nlrs"                           = "#E8C5A0", 
  "pmk1_pathway"                   = "#D4B8D8", 
  "septin_oligomerization"         = "#C8D8E0"
)

# hex buidler
hex2rgb <- function(hex) {
  h <- gsub("#", "", hex)
  as.integer(c(strtoi(substr(h, 1, 2), 16L),
               strtoi(substr(h, 3, 4), 16L),
               strtoi(substr(h, 5, 6), 16L)))
}

# Build edge table (shared Function_broad_final) 
# Two proteins are joined by an edge when they belong to the
# same functional group.  Groups with > 25 members would
# produce a lot edges each, so we random-sample
# MAX_EDGES_PER_GROUP pairs instead (for this I used 800 to make the hub edges a bit more denser - number of nodes will not be affected by this).  If this number is raised the network will show denser intra-group connections. If number is reduced, layout will be speedy with less connections.

MAX_EDGES_PER_GROUP <- 800
set.seed(148)

edge_list <- lapply(
  split(dat$Id, dat$Function_broad_final),
  function(pids) {
    n <- length(pids)
    if (n < 2) return(NULL)

    if (n * (n - 1) / 2 <= MAX_EDGES_PER_GROUP) {
      # small group  →  full clique
      idx <- combn(n, 2)                        # 2 × C(n,2)
      data.frame(source = pids[idx[1, ]],
                 target = pids[idx[2, ]],
                 stringsAsFactors = FALSE)
    } else {
      # large group  →  random sample of pairs
      pairs <- t(replicate(MAX_EDGES_PER_GROUP, sample(pids, 2)))
      data.frame(source = pairs[, 1],
                 target = pairs[, 2],
                 stringsAsFactors = FALSE)
    }
  }
)

edges <- do.call(rbind, edge_list)

# deduplicate undirected edges (A→B==B→A)
edges <- edges %>%
  mutate(s = pmin(source, target),
         t = pmax(source, target)) %>%
  distinct(s, t) %>%
  dplyr::rename(source = s, target = t)

# this will give a warning if any Function_broad_final is missing from palette
missing_cats <- setdiff(unique(dat$Function_broad_final), names(palette))
if (length(missing_cats) > 0) {
  warning("These categories have no colour assigned and will appear grey:\n  ",
          paste(missing_cats, collapse = "\n  "))
  # fill missing with a neutral grey so the code and network doesnt crash
  missing_hex <- rep("#D3D3D3", length(missing_cats))
  names(missing_hex) <- missing_cats
  palette <- c(palette, missing_hex)
}

# Create GEXF XML
# Node elements (one per protein)
node_fmt <- paste0(
  '    <node id="%s" label="%s">\n',
  '      <attvalues>\n',
  '        <attvalue for="0" value="%s"/>\n',
  '      </attvalues>\n',
  '      <viz:color r="%d" g="%d" b="%d"/>\n',
  '    </node>'
)

node_blocks <- mapply(function(pid, func) {
  rgb <- hex2rgb(palette[func])
  sprintf(node_fmt, pid, pid, func, rgb[1], rgb[2], rgb[3])
}, dat$Id, dat$Function_broad_final, SIMPLIFY = TRUE)


# Edge format: added label, type, and weight attributes
edge_fmt <- '      <edge id="%s" label="%s" source="%s" target="%s" type="undirected" weight="1.0" />'
eids        <- as.character(seq_len(nrow(edges)))  # IDs as strings, starting from 1
edge_blocks <- mapply(function(eid, src, tgt) sprintf(edge_fmt, eid, eid, src, tgt),
                      eids, edges$source, edges$target, SIMPLIFY = TRUE)
# write and save the file
gexf_xml <- paste0(
  '<?xml version="1.0" encoding="UTF-8"?>\n',
  '<gexf xmlns="http://www.gexf.net/1.3" ',
        'xmlns:viz="http://www.gexf.net/1.3/viz" version="1.3">\n', 
  '  <meta lastmodifieddate="', Sys.Date(), '">\n',
  '    <creator>R</creator>\n',
  '  </meta>\n',
  '  <graph defaultedgetype="undirected">\n',
  '    <attributes class="node" mode="static">\n',
  '      <attribute id="0" title="Function_broad_final" type="string"/>\n',
  '    </attributes>\n',
  '    <nodes>\n',
       paste(node_blocks, collapse = "\n"), "\n",
  '    </nodes>\n',
  '    <edges>\n',
       paste(edge_blocks, collapse = "\n"), "\n",
  '    </edges>\n',
  '  </graph>\n',
  '</gexf>\n'
)

writeLines(gexf_xml, here("data/results/005_gene_networks/","gephi_network.gexf"))

cat("✓  gephi_network.gexf written\n",
    "   Nodes :", nrow(dat),  "\n",
    "   Edges :", nrow(edges), "\n")
```

    ## ✓  gephi_network.gexf written
    ##     Nodes : 2742 
    ##     Edges : 13067

This is the raw netowrk, and will look like a blob when opened in
Gephi - as shown below:

![](../data/results/005_gene_networks/00_figures/00_raw_network.png)

## Base Network

This network was then force directed based on the ForceAtlas2 layout
parameters run, the set layout was saved at
`data/results/005_gene_networks/01_network_with_genenames.gexf.gephi`.
Since Gephi doesn’t allow adding text/annotations outside the nodes, I
added the functional categories for the network hubs in powerpoint.

![](../data/results/005_gene_networks/00_figures/01_network_with_genenames.png)

## Mapping differentially expressed genes on the Base Network

The above network was the base network layout that was then used to
visualize different categories of upreg and downreg genes (all performed
on Gephi using the node sizing option). The information required for
sizing the nodes is saved at
`data/processed/network_base_manual_with_degs.tsv`. The set layout with
node sizing info for DEGs (in data laboratory) was saved at
`data/results/005_gene_networks/02_network_fin.gexf.gephi`.

For showing the regulation in Bip1, Mst12 or both, I used the colors
consistent with previous analyses:

``` r
Differentially regulated in Dbip1="#fdc086"
Differentially regulated in Dmst12="#beaed4"
Differentially regulated in Both Dbip1 and Dmst12 : "#ffff99"
```

### Downreg in Guy11 at 0h (Upregulated in Bip1/Mst12 at 0h)

![](../data/results/005_gene_networks/00_figures/downreg_in_guy11_0h.png)

### Downreg in Guy11 at 4h (Upregulated in Bip1/Mst12 at 4h)

![](../data/results/005_gene_networks/00_figures/downreg_in_guy11_4h.png)

### Downreg in Guy11 at 24h (Upregulated in Bip1/Mst12 at 24h)

![](../data/results/005_gene_networks/00_figures/downreg_in_guy11_24h.png)

### Upreg in Guy11 at 0h (Downregulated in Bip1/Mst12 at 0h)

![](../data/results/005_gene_networks/00_figures/upreg_in_guy11_0h.png)

### Upreg in Guy11 at 4h (Downregulated in Bip1/Mst12 at 4h)

![](../data/results/005_gene_networks/00_figures/upreg_in_guy11_4h.png)

### Upreg in Guy11 at 24h (Downregulated in Bip1/Mst12 at 24h)

![](../data/results/005_gene_networks/00_figures/upreg_in_guy11_24h.png)

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS 15.7.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.4 forcats_1.0.1   stringr_1.6.0   dplyr_1.1.4    
    ##  [5] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.3.0   
    ##  [9] ggplot2_3.5.2   tidyverse_2.0.0 here_1.0.2     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit_4.6.0          archive_1.1.12     gtable_0.3.6       crayon_1.5.3      
    ##  [5] compiler_4.3.1     tidyselect_1.2.1   parallel_4.3.1     dichromat_2.0-0.1 
    ##  [9] scales_1.4.0       yaml_2.3.10        fastmap_1.2.0      R6_2.6.1          
    ## [13] generics_0.1.4     knitr_1.50         rprojroot_2.1.1    pillar_1.11.1     
    ## [17] RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.6        stringi_1.8.7     
    ## [21] xfun_0.52          bit64_4.6.0-1      timechange_0.3.0   cli_3.6.5         
    ## [25] withr_3.0.2        magrittr_2.0.3     digest_0.6.37      grid_4.3.1        
    ## [29] vroom_1.6.5        rstudioapi_0.17.1  hms_1.1.4          lifecycle_1.0.4   
    ## [33] vctrs_0.6.5        evaluate_1.0.5     glue_1.8.0         farver_2.1.2      
    ## [37] rmarkdown_2.29     tools_4.3.1        pkgconfig_2.0.3    htmltools_0.5.8.1
