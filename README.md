Cancer prognosis with shallow tumor RNA sequencing
==================================================

Here, we describe a simplified version of the analysis that is at the
core of our paper
(<a href="doi:10.1038/s41591-019-0729-3" class="uri">https://www.nature.com/articles/s41591-019-0729-3</a>).
With the help of one example, we show how a dramatic reduction in
RNA sequencing depth has little to no impact on the performance of machine
learning-based linear Cox models that predict disease outcome based on
tumor gene expression.

Since this analysis is peformed in R, if you have not installed it yet,
you can follow the intructions in
<a href="https://cran.r-project.org/" class="uri">https://cran.r-project.org/</a>.

In case R is installed, it needs to be version 3.6.1 or higher for this
example to work. The following code can help determine if R needs to be
updated.

    if(sessionInfo()$R.version$version.string < '3.6.1'){
      stop(paste0("This will not run for R versions older than 3.6.1. ",
                  "Your version is: ",
                  sessionInfo()$R.version$version.string,
                  ". Please, update R and try again."))
    } else {
      paste0("Your version is: ",
             sessionInfo()$R.version$version.string,
             ". This R version should be able to handle this example")
    }

In this example, we will use adrenocortical carcinoma (ACC) to
demonstrate how a drastic reduction in RNA-seq depth still gives enough
information to predict the relative risk of adverse outcome of disease. You can change the cancer
type by changing “ACC” here to any of the standard cancer type name
abbreviations of TCGA:
<a href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations" class="uri">https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations</a>. Please keep in mind we did not perform the analysis for DLBC, KICH, PCPG, CHOL, SKCM and SARC as described in the methods our paper.

    type <- "ACC"

Importantly, in our paper, two outcomes of disease were used, overall
survival (OS) or progression-free interval (PFI), depending on cancer
type. TCGA’s reccommendations as in
<a href="doi:10.1016/j.cell.2018.02.052" class="uri">doi:10.1016/j.cell.2018.02.052</a>
were followed.

With the following code, PFI can be used as outcome of disease for the
appropirate cancer types.

    #define cancer types where progression-free interval should be used instead of overall survival
    PFI <- c("BRCA", "LGG",  "PRAD", "READ", "TGCT", "THCA", "THYM")

Next, a few packages need to be installed.

Depending on the internet connection and machine configuration, this can
take up to several minutes.

    tryCatch(library("caret"), 
             error = function(e){
               install.packages(pkgs =  "caret", 
                                repos = 'http://cran.us.r-project.org')
               library("caret")
             })

    ## Loading required package: lattice

    ## Loading required package: ggplot2

    tryCatch(library("openxlsx"), 
             error = function(e){
               install.packages(pkgs =  "openxlsx", 
                                repos = 'http://cran.us.r-project.org')
               library("openxlsx")
             })
    tryCatch(library("doParallel"), 
             error = function(e){
               install.packages(pkgs =  "doParallel", 
                                repos = 'http://cran.us.r-project.org')
               library("doParallel")
             })

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

    tryCatch(library("rms"), 
             error = function(e){
               install.packages(pkgs =  "rms", 
                                repos = 'http://cran.us.r-project.org')
               library("rms")
             })

    ## Loading required package: Hmisc

    ## Loading required package: survival

    ## 
    ## Attaching package: 'survival'

    ## The following object is masked from 'package:caret':
    ## 
    ##     cluster

    ## Loading required package: Formula

    ## 
    ## Attaching package: 'Hmisc'

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

    ## Loading required package: SparseM

    ## 
    ## Attaching package: 'SparseM'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

    tryCatch(library("dplyr"), 
             error = function(e){
               install.packages(pkgs =  "dplyr", 
                                repos = 'http://cran.us.r-project.org')
               library("dplyr")
             })

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:Hmisc':
    ## 
    ##     src, summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    tryCatch(library("survival"), 
             error = function(e){
               install.packages(pkgs =  "survival", 
                                repos = 'http://cran.us.r-project.org')
               library("survival")
             })
    tryCatch(library("glmnet"), 
             error = function(e){
               install.packages(pkgs =  "glmnet", 
                                repos = 'http://cran.us.r-project.org')
               library("glmnet")
             })

    ## Loading required package: Matrix

    ## Loaded glmnet 3.0-2

    tryCatch(library("SummarizedExperiment"), 
             error = function(e){
               if (!requireNamespace("BiocManager",
                                     quietly = TRUE))
                 install.packages("BiocManager", 
                                  repos = 'http://cran.us.r-project.org')
               BiocManager::install("SummarizedExperiment",
                                    update = FALSE,
                                    ask = FALSE)
               library("SummarizedExperiment")
             })

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     which

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     expand

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:Hmisc':
    ## 
    ##     contents

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

    tryCatch(library("TCGAbiolinks"), 
             error = function(e){
               if (!requireNamespace("BiocManager",
                                     quietly = TRUE))
                 install.packages("BiocManager", 
                                  repos = 'http://cran.us.r-project.org')
               BiocManager::install("TCGAbiolinks",
                                    update = FALSE,
                                    ask = FALSE)
               library("TCGAbiolinks")
             })
    tryCatch(library("biomaRt"), 
             error = function(e){
               if (!requireNamespace("BiocManager",
                                     quietly = TRUE))
                 install.packages("BiocManager", 
                                  repos = 'http://cran.us.r-project.org')
               BiocManager::install("biomaRt",
                                    update = FALSE,
                                    ask = FALSE)
               library("biomaRt")
             })
    tryCatch(library("subSeq"), 
             error = function(e){
               if (!requireNamespace("BiocManager",
                                     quietly = TRUE))
                 install.packages("BiocManager", 
                                  repos = 'http://cran.us.r-project.org')
               BiocManager::install("subSeq",
                                    update = FALSE,
                                    ask = FALSE)
               library("subSeq")
             })
    tryCatch(library("edgeR"), 
             error = function(e){
               if (!requireNamespace("BiocManager",
                                     quietly = TRUE))
                 install.packages("BiocManager", 
                                  repos = 'http://cran.us.r-project.org')
               BiocManager::install("edgeR",
                                    update = FALSE,
                                    ask = FALSE)
               library("edgeR")
             })

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    tryCatch(library("limma"), 
             error = function(e){
               if (!requireNamespace("BiocManager",
                                     quietly = TRUE))
                 install.packages("BiocManager", 
                                  repos = 'http://cran.us.r-project.org')
               BiocManager::install("limma",
                                    update = FALSE,
                                    ask = FALSE)
               library("limma")
             })

    ### It's really essential that your version of TCGAbiolinks is 2.12.3 or newer

    if(packageVersion("TCGAbiolinks") < '2.12.3'){
      stop(paste0("This will not run with versions of TCGAbiolinks older than 2.12.3.",
                  " Your version is: ",
                  packageVersion("TCGAbiolinks"),
                  ". Update TCGAbiolinks and try again.",
                  " Importantly, TCGAbiolinks 2.12.3 and higher only run on R 3.6.1."))
    } else {
      paste0("Your version is: ",
             packageVersion("TCGAbiolinks"),
             ". This version should be able to handle this example")
    }

    ## [1] "Your version is: 2.14.0. This version should be able to handle this example"

Now, we can get and pre-process the gene expression data.

    #get gene expression data
    query <- GDCquery(project = paste0("TCGA-",
                                       as.character(type)),
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - Counts",
                      legacy = FALSE)

    ## --------------------------------------

    ## o GDCquery: Searching in GDC database

    ## --------------------------------------

    ## Genome of reference: hg38

    ## --------------------------------------------

    ## oo Accessing GDC. This might take a while...

    ## --------------------------------------------

    ## ooo Project: TCGA-ACC

    ## --------------------

    ## oo Filtering results

    ## --------------------

    ## ooo By data.type

    ## ooo By workflow.type

    ## ----------------

    ## oo Checking data

    ## ----------------

    ## ooo Check if there are duplicated cases

    ## ooo Check if there results for the query

    ## -------------------

    ## o Preparing output

    ## -------------------

    GDCdownload(query, 
                method = "api",
                files.per.chunk = 10, 
                directory = "GDCdata")

    ## Downloading data for project TCGA-ACC

    ## Of the 79 files for download 79 already exist.

    ## All samples have been already downloaded

    data <- GDCprepare(query,
                       save = TRUE,
                       save.filename = paste0("RangSummExp.", 
                                              as.character(type),
                                              ".Rdata"))

    ## |                                                    |  0%                      |                                                    |1.265823% ~58 s remaining |=                                                   |2.531646% ~34 s remaining |=                                                   |3.797468% ~24 s remaining |==                                                  |5.063291% ~19 s remaining |===                                                 |6.329114% ~16 s remaining |===                                                 |7.594937% ~14 s remaining |====                                                |8.860759% ~12 s remaining |=====                                               |10.12658% ~11 s remaining |=====                                               |11.39241% ~15 s remaining |======                                              |12.65823% ~14 s remaining |=======                                             |13.92405% ~13 s remaining |=======                                             |15.18987% ~12 s remaining |========                                            |16.4557% ~11 s remaining  |=========                                           |17.72152% ~10 s remaining |=========                                           |18.98734% ~10 s remaining |==========                                          |20.25316% ~9 s remaining  |===========                                         |21.51899% ~9 s remaining  |===========                                         |22.78481% ~9 s remaining  |============                                        |24.05063% ~8 s remaining  |=============                                       |25.31646% ~8 s remaining  |=============                                       |26.58228% ~8 s remaining  |==============                                      |27.8481% ~7 s remaining   |===============                                     |29.11392% ~7 s remaining  |===============                                     |30.37975% ~7 s remaining  |================                                    |31.64557% ~7 s remaining  |=================                                   |32.91139% ~6 s remaining  |=================                                   |34.17722% ~6 s remaining  |==================                                  |35.44304% ~6 s remaining  |===================                                 |36.70886% ~6 s remaining  |===================                                 |37.97468% ~6 s remaining  |====================                                |39.24051% ~5 s remaining  |=====================                               |40.50633% ~5 s remaining  |=====================                               |41.77215% ~5 s remaining  |======================                              |43.03797% ~5 s remaining  |=======================                             |44.3038% ~5 s remaining   |=======================                             |45.56962% ~5 s remaining  |========================                            |46.83544% ~4 s remaining  |=========================                           |48.10127% ~4 s remaining  |=========================                           |49.36709% ~5 s remaining  |==========================                          |50.63291% ~5 s remaining  |==========================                          |51.89873% ~4 s remaining  |===========================                         |53.16456% ~4 s remaining  |============================                        |54.43038% ~4 s remaining  |============================                        |55.6962% ~4 s remaining   |=============================                       |56.96203% ~4 s remaining  |==============================                      |58.22785% ~4 s remaining  |==============================                      |59.49367% ~3 s remaining  |===============================                     |60.75949% ~3 s remaining  |================================                    |62.02532% ~3 s remaining  |================================                    |63.29114% ~3 s remaining  |=================================                   |64.55696% ~3 s remaining  |==================================                  |65.82278% ~3 s remaining  |==================================                  |67.08861% ~3 s remaining  |===================================                 |68.35443% ~3 s remaining  |====================================                |69.62025% ~2 s remaining  |====================================                |70.88608% ~2 s remaining  |=====================================               |72.1519% ~2 s remaining   |======================================              |73.41772% ~2 s remaining  |======================================              |74.68354% ~2 s remaining  |=======================================             |75.94937% ~2 s remaining  |========================================            |77.21519% ~2 s remaining  |========================================            |78.48101% ~2 s remaining  |=========================================           |79.74684% ~2 s remaining  |==========================================          |81.01266% ~1 s remaining  |==========================================          |82.27848% ~1 s remaining  |===========================================         |83.5443% ~1 s remaining   |============================================        |84.81013% ~1 s remaining  |============================================        |86.07595% ~1 s remaining  |=============================================       |87.34177% ~1 s remaining  |==============================================      |88.60759% ~1 s remaining  |==============================================      |89.87342% ~1 s remaining  |===============================================     |91.13924% ~1 s remaining  |================================================    |92.40506% ~1 s remaining  |================================================    |93.67089% ~1 s remaining  |=================================================   |94.93671% ~0 s remaining  |==================================================  |96.20253% ~0 s remaining  |==================================================  |97.46835% ~0 s remaining  |=================================================== |98.73418% ~0 s remaining  |====================================================|100% ~0 s remaining       |====================================================|100%                      Completed after 8 s

    ## Starting to add information to samples

    ##  => Add clinical information to samples

    ## Add FFPE information. More information at: 
    ## => https://cancergenome.nih.gov/cancersselected/biospeccriteria 
    ## => http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html

    ##  => Adding subtype information to samples

    ## acc subtype information from:doi:10.1016/j.ccell.2016.04.002

    ## Accessing www.ensembl.org to get gene information

    ## Downloading genome information (try:0) Using: Human genes (GRCh38.p13)

    ## Cache found

    ## From the 60483 genes we couldn't map 3984

    ## Saving file:RangSummExp.ACC.Rdata

    ## File saved

    count <- assay(data)[,colData(data)$shortLetterCode == "TP" |
                           colData(data)$shortLetterCode == "TB" |
                           colData(data)$shortLetterCode == "TBM"]
    map_ens_sym <- rowData(data)

    count <- count[!duplicated(rownames(count)),]

    #function to keep genes detected in at least 0.1% of samples and to get log2-counts per million
    log.cpm <- function(valid.count){
      vc.dge <- DGEList(counts = valid.count)
      vc.dge.isexpr <- rowSums(cpm(vc.dge) > 1) >= round(dim(vc.dge)[2]*0.001)
      vc.dge <- vc.dge[vc.dge.isexpr,]
      vc.dge <- calcNormFactors(vc.dge)
      vc.voom <- voom(vc.dge)
      vlc <- t(vc.voom$E)
      vlc <- vlc[complete.cases(vlc),]
      return(vlc)
    }

    logCPM <- log.cpm(count)

    #get rid of samples sequenced more than once
    duplicate.samples <-
      sort(rownames(logCPM)[
        duplicated(substr(rownames(logCPM), 
                          1, 
                          12)) |
          duplicated(substr(rownames(logCPM), 
                            1, 
                            12),
                     fromLast = TRUE)])

    logCPM <- logCPM[!rownames(logCPM) %in%
                       duplicate.samples[duplicated(substr(duplicate.samples, 
                                                           1, 
                                                           12))],]

In our paper, the updated disease outcomes data published by TCGA in
<a href="doi:10.1016/j.cell.2018.02.052" class="uri">doi:10.1016/j.cell.2018.02.052</a>
were used. Let’s acccess, load and pre-process these data here.

    #download outcome data from Liu et al. Cell 2018
    #freely accessible on PubMed Central: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/
    upd.Surv <- read.xlsx("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/bin/NIHMS978596-supplement-1.xlsx",
                          sheet = "TCGA-CDR") 
    #clean up data
    upd.Surv <- upd.Surv[,-1]
    upd.Surv$OS <- as.character(upd.Surv$OS) %>%
      as.numeric()
    upd.Surv$OS.time <- as.character(upd.Surv$OS.time) %>%
      as.numeric()
    upd.Surv$PFI <- as.character(upd.Surv$PFI) %>%
      as.numeric()
    upd.Surv$PFI.time <- as.character(upd.Surv$PFI.time) %>%
      as.numeric()

    #keep only data for cancer type analyzed here
    clin <- upd.Surv[upd.Surv$type == type
                     ,c("bcr_patient_barcode",
                        "OS", "OS.time",
                        "PFI", "PFI.time")] %>%
      droplevels(.)
    rownames(clin) <- clin$bcr_patient_barcode
    rm(upd.Surv)

    #clean up data
    if(type %in% PFI) {
      clin.cov <- colnames(clin)
      clin.cov[clin.cov == "PFI"] <- "status"
      clin.cov[clin.cov == "PFI.time"] <- "time"
      colnames(clin) <- clin.cov
    } else {
      clin.cov <- colnames(clin)
      clin.cov[clin.cov == "OS"] <- "status"
      clin.cov[clin.cov == "OS.time"] <- "time"
      colnames(clin) <- clin.cov
    }
    clin <- clin[!is.na(clin$time),]
    clin <- clin[clin$time > 0,]
    clin <- clin[substr(clin$bcr_patient_barcode, 
                        1, 
                        12) %in% 
                   substr(rownames(logCPM), 
                          1, 
                          12),]
    logCPM <- logCPM[substr(rownames(logCPM), 
                            1, 
                            12) %in%
                       substr(rownames(clin), 
                              1, 
                              12),]
    clin <- clin[
      match(substr(rownames(logCPM), 
                   1,
                   12),
            substr(clin$bcr_patient_barcode, 
                   1,
                   12)),
      c("bcr_patient_barcode",
        "time",
        "status")]

Now that we have loaded and preprocessed the raw data, we can start
training and testing our machine learning models. Remember, the aim is
to predict outcome of disease based on tumor gene expression data
generated by RNA-seq. We will do that using Cox proportional hazards
regression with an elastic net penalty.

Let’s create the indices of the samples which will be either in the
training set or in the test set.

    #create data split (50/50 split)
    testindex <- foreach(repetitions = 1:100) %do%{
      set.seed(repetitions + 2020)
      createFolds(clin[,"status"], k = 2)
    }

With the code aboove, as in our paper, we can create 100 different data
splits into training and testing samples. However, for computational
reasons, we will only perform the analysis for one of these 100
repetitions here. If desired, you can change the number below to chose a
different data split for training and testing. Here we picked repetition
number 42, but you can pick any from 1-100.

    #pick a data split (change to any number from 1 to 100 to run on a different data split)
    repetition <- 42

Now let’s use the indices created above to actually split our datasets.

    #select actual data split to be used here
    testindex <- lapply(testindex,
                        function(repetitions)
                          repetitions[[1]])
    trainindex <- seq(dim(clin)[1])[
      !seq(dim(clin)[1]) %in% 
        testindex[[repetition]]]

    #outcome of disease for test samples
    test.clin <- clin[
      testindex[[repetition]],]
    test.clin <- droplevels(test.clin)

    #gene expression for test samples
    testset <- logCPM[
      testindex[[repetition]],]

    #outcome of disease for training samples
    train.surv <- Surv(clin[trainindex, "time"], 
                       clin[trainindex, "status"])

    #gene expression for training samples
    trainset <- logCPM[
      trainindex,]

Before we can train our models on our training data, and test on our
test samples, we need to scale the gene expression data. As in our
paper, we first scaled the training data and then used the center and
scale of each gene in the training set to scale the test set. By doing
this we ensure that training and testing data are on the same scale.

    #scale gene expression of training samples
    trainset <- scale(trainset, 
                      center = TRUE, 
                      scale = TRUE)

    #scale gene expression of test samples using center and scale of train samples
    testset <- scale(testset, 
                     center = attr(trainset, "scaled:center"), 
                     scale = attr(trainset, "scaled:scale"))

Now it’s finally time to train our model! We will make a function
(build.model) that creates 5 cross-validation folds and feeds our
training algorithm with a range of alpha values. In case you’re not
familiar with these terms used here, this is a good start to find more
information:
<a href="https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html" class="uri">https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html</a>.

Our “build.model” function keeps the same cross-validation folds across
different alpha values to ensure that the performances of each alpha are
compared based on the same data. If you are using a unix system (but not
on windows), the function will run in parallel, but it will still take
several minutes to hours to train depending on the number of samples
used (adrenocortical carcinoma \[ACC\] has relatively few samples and
should run much faster than breast cancer \[BRCA\], for example). Also,
training takes a lot of memory! If you are running out of memory, make
sure to change the number of cores used in in “detectCores()-1” (go from
-1 to -2 or -3 to reduce the number of cores) and try again.

    #function to train elastic net Cox model
    build.model <- function(scaled.log.cpm, surv) {
      set.seed(2020)
      fold.id <-
        createFolds(surv[,2], k = 5, list = FALSE)
      alpha <- c(0, 10^seq(-5, -1, 1), seq(0.2, 0.9, 0.1), c(0.95, 0.99, 1))
      model <- mclapply(alpha,
                        function(a)
                          cv.glmnet(x = scaled.log.cpm,
                                            y = surv,
                                            family = "cox",
                                            type.measure = "deviance",
                                            alpha = a,
                                            foldid = fold.id,
                                            parallel = FALSE,
                                            standardize = FALSE),
                        mc.cores = ifelse(Sys.info()[['sysname']] == "Windows",
                                          yes = 1,
                                          no = (detectCores()-1)))
      names(model) <- alpha
      best.alpha <- lapply(model, 
                           function(x)
                             min(x$cvm)) %>%
        unlist(.) %>%
        which.min(.) %>%
        names(.)
      best.model <- model[[best.alpha]]
      best.model$"best.alpha" <- best.alpha
      return(best.model)
    }

    #actual training of the elastic net Cox model
    model <- build.model(scaled.log.cpm = trainset,
                         surv = train.surv)

Now we can predict the relative risk of death (or relative risk of
recurrence for cancer types that use PFI as measure of outcome) for
samples in the test set, and see how the prediction compares to actual
survival. We will do that with Cox regression after testing for the
proportional hazards assumption.

    #predict relative risk of event (RRE) using enet model
    test.clin$pred.resp <- 
      predict(model,
              newx = testset,
              s = "lambda.min",
              type = "response") %>%
      log(.) %>%
      .[,1]

    #build validation model using RRE
    cox.model <- coxph(Surv(time, status) ~
                         pred.resp,
                       data = test.clin)

    #test proportional hazards assumption
    cox.zph(cox.model)

    ##            chisq df    p
    ## pred.resp 0.0581  1 0.81
    ## GLOBAL    0.0581  1 0.81

    #since alpha < 0.05, check validation model
    summary(cox.model)

    ## Call:
    ## coxph(formula = Surv(time, status) ~ pred.resp, data = test.clin)
    ## 
    ##   n= 39, number of events= 15 
    ## 
    ##             coef exp(coef) se(coef)     z Pr(>|z|)    
    ## pred.resp 0.8889    2.4325   0.2510 3.541 0.000398 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##           exp(coef) exp(-coef) lower .95 upper .95
    ## pred.resp     2.432     0.4111     1.487     3.978
    ## 
    ## Concordance= 0.791  (se = 0.058 )
    ## Likelihood ratio test= 14.76  on 1 df,   p=1e-04
    ## Wald test            = 12.54  on 1 df,   p=4e-04
    ## Score (logrank) test = 14.85  on 1 df,   p=1e-04

In this example, gene expression RNA-seq data can be used to predict
outcome of disease in adrenocortical carcinoma with a concordance index
of 0.79 and a p-value (likelihood ratio test) of 1e-04.

Next we can subsample our count matrix to simulate a 100-fold reduction
in sequencing depth and see how that impacts predictive performance. The
fold reduction can be controlled here by changing “proportion” in
“generateSubsampledMatrix” (for example, proportion = 0.001 for a
1000-fold reduction, or proportion = 0.1 for 10-fold reduction).

    #subsample gene expression data
    #here, 100-fold reduction was used
    sub.count <- generateSubsampledMatrix(counts = count,
                                          proportion = 0.01,
                                          seed = 2020)
    sub.logCPM <- log.cpm(sub.count)
    duplicate.samples <-
      sort(rownames(sub.logCPM)[
        duplicated(substr(rownames(sub.logCPM), 
                          1, 
                          12)) |
          duplicated(substr(rownames(sub.logCPM), 
                            1, 
                            12),
                     fromLast = TRUE)])
    sub.logCPM <- sub.logCPM[!rownames(sub.logCPM) %in%
                               duplicate.samples[duplicated(substr(duplicate.samples, 
                                                                   1, 
                                                                   12))],]
    sub.logCPM <- sub.logCPM[substr(rownames(sub.logCPM), 
                                    1, 
                                    12) %in%
                               substr(rownames(clin), 
                                      1, 
                                      12),]

After splitting and scaling the subsampled data, we can build a model
that takes the subsampled gene expression data as input and aims at
predicting the relative risk of death (or recurrence, depending on
outcome measure used).

    sub.testset <- sub.logCPM[
      testindex[[repetition]],]

    sub.trainset <- sub.logCPM[
      trainindex,]

    sub.trainset <- scale(sub.trainset, 
                          center = TRUE, 
                          scale = TRUE)

    sub.testset <- scale(sub.testset, 
                         center = attr(sub.trainset, 
                                       "scaled:center"), 
                         scale = attr(sub.trainset, 
                                      "scaled:scale"))

    sub.model <- build.model(scaled.log.cpm = sub.trainset,
                             surv = train.surv)

With the model built on subsampled gene expression data we can try to
predict the relative risk of death (or recurrence) of patients in the
subsampled test set. Just as for full-coverage models, we can test
whether the predicted relative risk correlates with actual outcome using
Cox regression if the ph-assumption is not violated.

    test.clin$sub.pred.resp <- 
      log(predict(sub.model,
                  newx = sub.testset,
                  s = "lambda.min",
                  type = "response"))

    sub.cox.model <- coxph(Surv(time, status) ~
                             sub.pred.resp,
                           data = test.clin)
    cox.zph(sub.cox.model)

    ##               chisq df    p
    ## sub.pred.resp 0.115  1 0.73
    ## GLOBAL        0.115  1 0.73

    #since alpha < 0.05, check validation model
    summary(sub.cox.model)

    ## Call:
    ## coxph(formula = Surv(time, status) ~ sub.pred.resp, data = test.clin)
    ## 
    ##   n= 39, number of events= 15 
    ## 
    ##                 coef exp(coef) se(coef)     z Pr(>|z|)    
    ## sub.pred.resp 0.9382    2.5555   0.2529 3.709 0.000208 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##               exp(coef) exp(-coef) lower .95 upper .95
    ## sub.pred.resp     2.555     0.3913     1.557     4.195
    ## 
    ## Concordance= 0.818  (se = 0.054 )
    ## Likelihood ratio test= 16.37  on 1 df,   p=5e-05
    ## Wald test            = 13.76  on 1 df,   p=2e-04
    ## Score (logrank) test = 16.68  on 1 df,   p=4e-05

Please note that gene expression RNA-seq data at much shallower depths
could also be used to predict outcome of disease in adrenocortical
carcinoma. The concordance index was 0.82 and the p-value (likelihood
ratio test) was 5e-05.

These results demonstrate, with the help of an single example, that a
strong reduction on sequencing depth does not strongly impact the
performance of models that take RNA-seq data and try to predict disease
outcome.

    #performance stats
    ifelse(summary(sub.cox.model)$coef < 1e-3,
           yes = formatC(summary(sub.cox.model)$coef,
                         format = "e",
                         digits = 1),
           no = round(summary(sub.cox.model)$coef,
                      digits = 2) %>%
             as.numeric())

    ##               coef   exp(coef) se(coef) z      Pr(>|z|) 
    ## sub.pred.resp "0.94" "2.56"    "0.25"   "3.71" "2.1e-04"

    ifelse(summary(cox.model)$coef < 1e-3,
           yes = formatC(summary(cox.model)$coef,
                         format = "e",
                         digits = 1),
           no = round(summary(cox.model)$coef,
                      digits = 2) %>%
             as.numeric())

    ##           coef   exp(coef) se(coef) z      Pr(>|z|) 
    ## pred.resp "0.89" "2.43"    "0.25"   "3.54" "4.0e-04"

    #median library size of original dataset
    count %>% t(.) %>% rowSums() %>% median(.) %>% formatC(., format = "e", digits = 1) %>% paste0("median library size of original dataset: ", ., " counts")

    ## [1] "median library size of original dataset: 4.1e+07 counts"

    #median library size of subsampled dataset
    sub.count %>% t(.) %>% rowSums() %>% median(.) %>% formatC(., format = "e", digits = 1) %>% paste0("median library size of subsampled dataset: ", ., " counts")

    ## [1] "median library size of subsampled dataset: 4.1e+05 counts"

In case, this analysis can be run on any cancer type available in TCGA
(see second code chunk in the beginning; make sure to pick the correct
disease outcome \[OS vs PFI\]). In addition, any repetition from 1 to
100 can be run to see whether the results are biased due to unequal
sample distributions.
