if(sessionInfo()$R.version$version.string < '3.6.1'){
  stop(paste0("\n##### This will not run for R versions older than 3.6.1 #####\n",
              "##### your version is: ",
              sessionInfo()$R.version$version.string,
              " #####\n",
              "##### update R and rerun"))
}

#load (and install, if necessary) packages
tryCatch(library("caret"), 
         error = function(e){
           install.packages(pkgs =  "caret", 
                            repos = 'http://cran.us.r-project.org')
           library("caret")
         })
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
tryCatch(library("rms"), 
         error = function(e){
           install.packages(pkgs =  "rms", 
                            repos = 'http://cran.us.r-project.org')
           library("rms")
         })
tryCatch(library("dplyr"), 
         error = function(e){
           install.packages(pkgs =  "dplyr", 
                            repos = 'http://cran.us.r-project.org')
           library("dplyr")
         })
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

if(packageVersion("TCGAbiolinks") < '2.12.3'){
  stop(paste0("\n##### This will not run with versions of TCGAbiolinks older than 2.12.3 #####\n",
              "##### your version is: ",
              packageVersion("TCGAbiolinks"),
              " #####\n",
              "##### update TCGAbiolinks and rerun #####\n",
              "##### TCGAbiolinks 2.12.3 runs only on R 3.6.1 #####\n"))
}

#In this example, we will run the core analysis with ACC
type <- "ACC"

#define cancer types where progression-free interval should be used instead of overall survival
PFI <- c("BRCA", "LGG",  "PRAD", "READ", "TGCT", "THCA", "THYM")

#get gene expression data
query <- GDCquery(project = paste0("TCGA-",
                                   as.character(type)),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type  = "HTSeq - Counts",
                  legacy = FALSE)
GDCdownload(query, 
            method = "api",
            files.per.chunk = 10, 
            directory = "GDCdata")

data <- GDCprepare(query,
                   save = TRUE,
                   save.filename = paste0("RangSummExp.", 
                                          as.character(type),
                                          ".Rdata"))

count <- assay(data)[,colData(data)$shortLetterCode == "TP" |
                       colData(data)$shortLetterCode == "TB" |
                       colData(data)$shortLetterCode == "TBM"]
map_ens_sym <- rowData(data)

count <- count[!duplicated(rownames(count)),]

#function to keep genes detected in 0.1% of samples and get log2-counts per million
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

#create data split
testindex <- foreach(repetitions = 1:100) %do%{
  set.seed(repetitions + 2017)
  createFolds(clin[,"status"], k = 2)
}

#change 1 to any number from 1 to 100 to run a different data split
repetition <- 42

#select actual data split to be used here
testindex <- lapply(testindex,
                    function(repetitions)
                      repetitions[[1]])
trainindex <- seq(dim(clin)[1])[
  !seq(dim(clin)[1]) %in% 
    testindex[[repetition]]]

#outcome of test samples
test.clin <- clin[
  testindex[[repetition]],]
test.clin <- droplevels(test.clin)
#gene expression of test samples
testset <- logCPM[
  testindex[[repetition]],]

#outcome of training samples
train.surv <- Surv(clin[trainindex, "time"], 
                   clin[trainindex, "status"])

#gene expression of training samples
trainset <- logCPM[
  trainindex,]

#scale gene expression of training samples
trainset <- scale(trainset, 
                  center = T, 
                  scale = T)

#scale gene expression of test samples using center and scale of train samples
testset <- scale(testset, 
                 center = attr(trainset, "scaled:center"), 
                 scale = attr(trainset, "scaled:scale"))

#function to train elastic net Cox model
build.model <- function(scaled.log.cpm, surv) {
  set.seed(2018)
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

#since alpha < 0.05, check validation model
summary(cox.model)

#subsample gene expression data
#here, 100-fold reduction was used
sub.count <- generateSubsampledMatrix(counts = count,
                                      proportion = 0.01,
                                      seed = 2019)
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
sub.testset <- sub.logCPM[
  testindex[[repetition]],]

sub.trainset <- sub.logCPM[
  trainindex,]

sub.trainset <- scale(sub.trainset, 
                      center = T, 
                      scale = T)

sub.testset <- scale(sub.testset, 
                     center = attr(sub.trainset, 
                                   "scaled:center"), 
                     scale = attr(sub.trainset, 
                                  "scaled:scale"))

sub.model <- build.model(scaled.log.cpm = sub.trainset,
                         surv = train.surv)

test.clin$sub.pred.resp <- 
  log(predict(sub.model,
              newx = sub.testset,
              s = "lambda.min",
              type = "response"))

sub.cox.model <- coxph(Surv(time, status) ~
                         sub.pred.resp,
                       data = test.clin)
cox.zph(sub.cox.model)

#since alpha < 0.05, check validation model
summary(sub.cox.model)$coef
#pretty close to the resut of validation of model trained on full set:
summary(cox.model)$coef

#even though the library size is reduced by 100-fold:
#median library size of original dataset
count %>% t(.) %>% rowSums() %>% median(.)
#median library size of subsampled dataset
sub.count %>% t(.) %>% rowSums() %>% median(.)
