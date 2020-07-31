miReact example
================
Rasmus Rydbirk
31-07-2020

# Downsampled Tabula Muris example

Setup

``` r
library(miReact)
library(magrittr)
library(pagoda2)
```

    ## Loading required package: Matrix

    ## Loading required package: igraph

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(ggplot2)
```

Here, we assume that motif p-values and motif counts have been prepared
according to the tutorial, or the data have been downloaded.

First, we load the downsampled Tabula Muris count matrix. Since
scRNA-seq data are sparse by nature, we make it sparse to save time and
memory. Then, we prepare our data object with counts and run parameters:

``` r
cm <- readRDS("~/miReact/data/mm.exp1000downsample.rds") %>% Matrix::drop0()

runparameters <- list(motifs = "7",
                      alpha = 1e-10)

sco <- list(exp = cm,
            runparameters = runparameters)
```

Next, we run the prepareData() wrapper function to add motif sequences,
patterns, p-values and counts:

``` r
sco <- prepareData(sco = sco,
                    seqs.path = "~/miReact/seqs/mm.utr3.seqs.rds", 
                    seqlist.path = "~/miReact/seqs/mm.utr3.seqlist.rds", 
                    patterns.path = "~/miReact/motif.models/patterns.7mer.Rdata", 
                    pval.path = "~/miReact/motif.models/mm.seqXmot.utr3_mrs_7mer.rds", 
                    counts.path = "~/miReact/motif.models/mm.seqXmot.counts.utr3_mrs_7mer.rds", 
                    verbose = T)
```

    ## 2020-07-31 03:45:29 Loading data ... done!
    ##  2020-07-31 03:46:01 Calculating medians ... done!
    ##  2020-07-31 03:46:01 Ordering ... 
    ## 2020-07-31 03:46:03 All done!

Then, we calculate the miRNA activity:

``` r
res <- miReactC(sco, n.cores = 50, verbose = T)
```

    ## 2020-07-31 03:46:03 Analyzing with 50 cores in 20 chunks ... 
    ## 
    ##  2020-07-31 03:58:45 All done!

Next, we try to cluster the cells based on the estimated miRNA activity.
For the example dataset with 1,000 cells we use Pagoda2, but for larger
datasets it is possible to subdivide cells back into multiple samples
for more sophisticated analyses, e.g., using Conos or Seurat.

Since negative activity is arbitrary, and since single-cell pipelines
are not designed to handle negative numbers, we set negative numbers to
0 before making data sparse:

``` r
res0 <- res
res0[res0 < 0] <- 0
res0 %<>% Matrix::drop0()
```

Then, we cluster our data:

``` r
p2 <- Pagoda2$new(res0, 
                  log.scale=T, 
                  n.cores=50)
```

    ## 1000 cells, 16384 genes; normalizing ... using plain model log scale ... done.

``` r
p2$adjustVariance(plot=F, 
                  gam.k=10)
```

    ## calculating variance fit ... using gam 1340 overdispersed genes ... 1340persisting ... done.

``` r
p2$calculatePcaReduction(nPcs=50, 
                         n.odgenes=3e3)
```

    ## running PCA using 3000 OD genes .... done

``` r
p2$makeKnnGraph(k=40,type='PCA', 
                center=T, 
                distance='cosine', 
                n.cores=50)
p2$getKnnClusters(method=infomap.community, 
                  type='PCA', 
                  n.cores=50)
p2$getEmbedding(type='PCA', 
                embeddingType='tSNE', 
                perplexity=50, 
                verbose=F, 
                n.cores=50)
```

    ## calculating distance ... pearson ...running tSNE using 50 cores:

We load the annotation included in miReact:

``` r
annotation <- readRDS("~/miReact/data/mm.annotations1000downsample.rds")
annotation.tissue <- annotation$tissue %>% setNames(annotation$cell)
```

We plot the annotation:

``` r
p2$plotEmbedding(type='PCA',
                 embeddingType='tSNE',
                 mark.clusters=T,
                 min.group.size=1,
                 mark.cluster.cex=1,
                 alpha=0.5,
                 main='clusters (tSNE)', 
                 groups = annotation.tissue)
```

    ## using provided groups as a factor

![](figure-gfm/unnamed-chunk-8-1.png)<!-- -->

We see that clustering based on miRNA activity provides a fair
distinction between most of the cell types. Let’s plot an a priori known
liver-specific miRNA in our embedding:

``` r
p2$plotEmbedding(type='PCA', 
                 embeddingType='tSNE', 
                 colors=p2$counts[,"ACACTCC"], 
                 main="miR-122-5p activity")
```

    ## treating colors as a gradient with zlim: 0.0112585 0.5542113

![](/tmp/Rtmp1uc4rL/preview-d12713e4e23d.dir/miReactC_vignette_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

We can also show this with a dot plot:

``` r
plot.df <- data.frame(Activity = p2$counts[,"ACACTCC"], anno = annotation.tissue)

ggplot(plot.df, aes(anno, Activity, col=anno)) + 
  geom_jitter() +
  labs(title="miR-122-5p activity", x="", y="Normalized activity") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=90))
```

![](/tmp/Rtmp1uc4rL/preview-d12713e4e23d.dir/miReactC_vignette_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Further, we can try to identify cell type-specific miRNAs:

``` r
p2$getDifferentialGenes(groups = annotation.tissue, 
                        upregulated.only = T, 
                        verbose = T, 
                        append.specificity.metrics = T, 
                        append.auc = T)
```

    ## running differential expression with  20  clusters ... adjusting p-values ... done.

Let’s look at the best markers for liver cells:

``` r
p2$diffgenes$counts$customClustering$Liver %>% 
  dplyr::filter(Specificity > 0.7, ExpressionFraction > 0.5) %>% 
  dplyr::arrange(desc(Precision)) %>% 
  head(10)
```

    ##                 Z        M highest        fe    Gene Specificity Precision
    ## ACGACTC  7.582282 2.879764    TRUE 0.5140187 ACGACTC   0.9693122 0.6547619
    ## CGCCTCT  9.774723 2.795655    TRUE 0.6542056 CGCCTCT   0.9548387 0.6250000
    ## ATATATA  9.629946 2.972329    TRUE 0.6355140 ATATATA   0.9527897 0.6071429
    ## CACCCCA  8.410033 3.031692    TRUE 0.5607477 CACCCCA   0.9585106 0.6060606
    ## ACAAACA  8.620456 2.876478    TRUE 0.5794393 ACAAACA   0.9530917 0.5849057
    ## GCCCGCC  9.836744 2.924600    TRUE 0.6542056 GCCCGCC   0.9462366 0.5833333
    ## GGCGCGG  7.906069 3.020321    TRUE 0.5327103 GGCGCGG   0.9565217 0.5816327
    ## TCCGTTT  7.529671 2.719547    TRUE 0.5233645 TCCGTTT   0.9555085 0.5714286
    ## ACTCCCG  8.007992 2.949962    TRUE 0.5420561 ACTCCCG   0.9532909 0.5686275
    ## CCACCCC 10.852452 2.714102    TRUE 0.7289720 CCACCCC   0.9305857 0.5492958
    ##         ExpressionFraction
    ## ACGACTC          0.5140187
    ## CGCCTCT          0.6542056
    ## ATATATA          0.6355140
    ## CACCCCA          0.5607477
    ## ACAAACA          0.5794393
    ## GCCCGCC          0.6542056
    ## GGCGCGG          0.5327103
    ## TCCGTTT          0.5233645
    ## ACTCCCG          0.5420561
    ## CCACCCC          0.7289720

We plot the top hit in the embedding:

``` r
p2$plotEmbedding(type='PCA', 
                 embeddingType='tSNE', 
                 colors=p2$counts[,"ACGACTC"], 
                 main="ACGACTC activity")
```

    ## treating colors as a gradient with zlim: 0 0.02275152

![](/tmp/Rtmp1uc4rL/preview-d12713e4e23d.dir/miReactC_vignette_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

And in a dot plot:

``` r
plot.df <- data.frame(Activity = p2$counts[,"ACGACTC"], anno = annotation.tissue)

ggplot(plot.df, aes(anno, Activity, col=anno)) + 
  geom_jitter() +
  labs(title="ACGACTC activity", x="", y="Normalized activity") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=90))
```

![](/tmp/Rtmp1uc4rL/preview-d12713e4e23d.dir/miReactC_vignette_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Lastly, let’s plot a heatmap of the top markers for all cell types to
investigate how well they distinguish our cells:

``` r
genes <- sapply(p2$diffgenes$counts$customClustering, function(x) rownames(x)[1:5]) %>% 
  c() %>% 
  .[!is.na(.)]

p2$plotGeneHeatmap(genes=genes, 
                   groups=annotation.tissue, 
                   gradient.range.quantile = 0.9)
```

![](/tmp/Rtmp1uc4rL/preview-d12713e4e23d.dir/miReactC_vignette_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_3.3.2 pagoda2_0.1.1 igraph_1.2.5  Matrix_1.2-18 magrittr_1.5 
    ## [6] miReact_0.2  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0        xfun_0.15               Rook_1.1-1             
    ##  [4] purrr_0.3.4             pbapply_1.4-2           splines_4.0.2          
    ##  [7] lattice_0.20-41         colorspace_1.4-1        vctrs_0.3.2            
    ## [10] generics_0.0.2          expm_0.999-4            htmltools_0.5.0        
    ## [13] sccore_0.1              yaml_2.2.1              mgcv_1.8-31            
    ## [16] base64enc_0.1-3         rlang_0.4.7             pillar_1.4.6           
    ## [19] glue_1.4.1              withr_2.2.0             matrixStats_0.56.0     
    ## [22] lifecycle_0.2.0         stringr_1.4.0           MatrixGenerics_1.0.2   
    ## [25] munsell_0.5.0           gtable_0.3.0            evaluate_0.14          
    ## [28] labeling_0.3            knitr_1.29              irlba_2.3.3            
    ## [31] parallel_4.0.2          Regmex_1.0              urltools_1.7.3         
    ## [34] triebeard_0.3.0         Rcpp_1.0.5              scales_1.1.1           
    ## [37] RcppParallel_5.0.2-9000 farver_2.0.3            brew_1.0-6             
    ## [40] rjson_0.2.20            digest_0.6.25           Rtsne_0.15             
    ## [43] stringi_1.4.6           dplyr_1.0.0             grid_4.0.2             
    ## [46] tools_4.0.2             tibble_3.0.3            crayon_1.3.4           
    ## [49] pkgconfig_2.0.3         MASS_7.3-51.6           ellipsis_0.3.1         
    ## [52] sparseMatrixStats_1.0.5 rmarkdown_2.3           R6_2.4.1               
    ## [55] nlme_3.1-148            compiler_4.0.2
