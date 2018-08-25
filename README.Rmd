# scRNA-seq Iterative Pooling with DESeq2 

Single cell transcriptomics (scRNA-seq) has in recent years become a popular method to assay heterogeneity of cells in a cell population and uncover hidden subpopulations. In contrast to bulk RNA-seq, scRNA-seq is thought to suffer from large heterogeneity between samples and "dropout events" (false zeros). Moreover, over-dispersion and technical variation are more evident as gene expression is inherently stochasatic in nature. 

Many packages exist allowing to analyze and cluster scRNA-seq data produced using different methods. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a commonly used R package for identifying differentially expressed genes based on a model using the negative binomial distribution. scRNA-seq reads are known to follow a Poisson distribution, and therefore utilising DESeq2 and similar packages destined for bulk RNA-seq analysis is not recommended. To determine the effect of using such methods, an iterative method similar to a "bootstrap" was implemented using DESeq2. 

The goal was to determine whether individual cells expressing the transcription factor FoxD3 in the zebrafish neural crest at 5-6ss somite stage carried out in the TSS lab and single cells expressing FoxD3 at 50% epiboly co-expressed similar sets of genes and to what extent these could be uncovered using DESeq2. Zebrafish scRNA-seq 50% epiboly dataset was kindly provided by R.Satija. The original study describing the 50% epiboly dataset was published in :
  
*Satija, R., Farrell, J.A., Gennert, D., Schier, A.F. and Regev, A., 2015. Spatial reconstruction of single-cell gene expression data. Nature biotechnology, 33(5), p.495.[doi:10.1038/nbt.3192](https://www.nature.com/articles/nbt.3192)*


A concatenated fasta file genome containing both the danRer10 sequence and ERCCs and a corresponding gtf file was produced (data/Zv10_ERCCs.gtf). A STAR (v.2.4.2) index was made using default parameters. 

Paired end reads were mapped to the zebrafish genome (danRer10) using STAR using the following parameters:

```{r eval = FALSE} 
STAR --genomeDir $GDIR --readFilesIn $A1 $A2 --runThreadN 8 --outFileNamePrefix $AN --outSAMstrandField intronMotif
```

Read counts were obtained using subread featureCounts (v.1.4.6) package using the following parameters: 


```{r eval = FALSE} 
featureCounts -T 16 -p -t exon -g gene_id -a Zv10_ERCCs.gtf -o $NAME
```

This repository contains the count tables obtained using featureCounts for both single cell datasets in data/

The raw reads from the neural crest dataset has been uploaded to NCBI's GEO under accession number GSE106676 that will be rendered public upon pear-reviewed publication of manuscript. 
  
All of the cells in the neural crest single cell dataset were FACs-sorted based on FoxD3-Citrine reporter expression. The 50% epiboly dataset however contained all cells from this stage, many of which did not express FoxD3. Cells expressing FoxD3 transcript XXX => 1 FPKM were considered as FoxD3-expressing cells. 201 cells out of the total 743 expressed FoxD3. One cell expressed FoxD3 < 1 FPKM. 

Raw reads were converted to FPKM values using  R (inspired from [this blog](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/)):

Low read counts and lack of reporoducibilty of replicate results in DESeq2's failure to estimation of the variance in the null distribution. Using the [fdrtool R library](https://cran.r-project.org/web/packages/fdrtool/index.html) this can be correct. For much better explanation of this link [here](http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html)

To run : 

```{r eval = FALSE} 
single_cell_iterative_pooling2.py -f STAR_mapped_zv10_single_cell_satija_counts -gi ENSDARG00000021032 -o FoxD3 -c 50  -l 1  -i 100 
```

This will use the featureCounts count table 'our_single_cell_with_ERCCs_zv10.counts.txt' and subselect all cells expressing FoxD3 (Ensembl geneid=ENSDARG00000021032). The limit to consider a cell as expressing FoxD3 is defined by -l of 1 FPKM. Using the -c 50 and -i 100, it will pool 50 cells together arbitrarily and repeat this 100 times. 

#### Reference manuscript

This was carried out as part of the FoxD3 project carried out in the laboratory of T. Sauka-Spengler at the Weatherall Institute of Molecular Medicine at the University of Oxford, in the UK. 

Preliminary manuscript for this project is available on bioRxiv: https://www.biorxiv.org/content/biorxiv/early/2017/11/22/213611.full.pdf

For additional information on project please see [tsslab/foxd3](https://github.com/tsslab/foxd3) github repository.  

