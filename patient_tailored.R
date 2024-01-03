# libraries and default functions
library(devtools)
install_github('oganm/geneSynonym')
lapply(c("parallel","geneSynonym"), library, character.only = T)
species_ = list(species_ = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Macaca mulatta"),
                code_ = c(9606, 10090, 10116, 7955, 9544))
selectedSpecies_ = "Homo sapiens"



# Get compound info (demo data)
Dr_Info_FIMM = openxlsx::read.xlsx("data/compound.xlsx")

### GET COMPOUND IDs, no need to uncomment in case of demo data
###############################################################

library(stringr)
library(httr)
library(purrr)
get_cid <- function(name) {
  url <- paste0(
    "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    str_replace_all(name, " ", ""),
    "/cids/JSON"
  )
  res <- GET(url)
  if (status_code(res) == 200) {
    cid <- content(res)[[1]][[1]][[1]]
    return(cid)
  } else {
    return("")
  }
}

get_cid_safe <- possibly(get_cid, otherwise = "")
Dr_Info_FIMM$CID <- map_chr(Dr_Info_FIMM$compounds_, get_cid_safe)
# Dr_Info_FIMM <- data.frame(compounds_ = as.character(unique(Dr_Info_FIMM$Compound)))
# Dr_Info_FIMM$CID = mclapply(Dr_Info_FIMM$compounds_, function(i){
#   tryCatch({
#     jsonlite::fromJSON(paste0("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",gsub(" ","",i),"/cids/JSON"))[[1]][[1]][[1]]
#   }, error = function(e) {
#     ""
#   })
# }, mc.cores = 1)

### GET TARGETS from DGIdb, no need to uncomment in case of demo data
##########################

Dr_Info_FIMM$TARGETS = mclapply(1:nrow(Dr_Info_FIMM), function(i){

  db = tryCatch({
    tmp_ = jsonlite::fromJSON(paste0("http://www.dgidb.org/api/v2/interactions.json?drugs=",gsub(" ","",tolower(as.character(Dr_Info_FIMM$compounds_[i][[1]])))))$matchedTerms$interactions[[1]]
    paste0(unique(tmp_[tmp_$score>2,]$geneName),collapse=",")
    }, error = function(e) {
    ""
  })
  if(length(db) == 0 || is.null(db) || db == ""){
    db = tryCatch({
      synonyms_ = tryCatch({
        tolower(jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",gsub(" ","",Dr_Info_FIMM$CID[i][[1]]),"/synonyms/JSON"))$InformationList$Information$Synonym[[1]])
      }, error = function(e) {
        ""
      })
      synonyms_ = synonyms_[as.numeric(sapply(synonyms_, nchar)) < 50]
      tmp_ = jsonlite::fromJSON(readLines(curl::curl(URLencode(paste0("http://www.dgidb.org/api/v2/interactions.json?drugs=",gsub(" ","",tolower(paste0(synonyms_,collapse = ","))))))))$matchedTerms$interactions[[1]]
      paste0(unique(tmp_[tmp_$score>2,]$geneName),collapse=",")
    }, error = function(e) {
      ""
    })
  }
  db
}, mc.cores = 1)


### GET TARGETS FROM DTC, no need to uncomment in case of demo data Drug target combinations
###################################################################

dataDTC = readr::read_delim("~/Downloads/DtcDrugTargetInteractions.csv", delim = ",")
dataDTC = dataDTC[toupper(dataDTC$standard_type) %in% c("KI","KD","IC50","EC50"),]
dataDTC = dataDTC[!is.na(dataDTC$gene_names) & !is.na(dataDTC$compound_name) & !is.na(dataDTC$standard_units),c(3,7,11,13,14)]
dataDTC = dataDTC[toupper(dataDTC$standard_units) == "NM", ]; dataDTC$standard_units = NULL; dataDTC = unique(dataDTC)

Dr_Info_FIMM$TARGETSDTC = mclapply(1:nrow(Dr_Info_FIMM), function(i){

  synonyms_ = toupper(tryCatch({
    c(Dr_Info_FIMM$compounds_[i][[1]], tolower(jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",gsub(" ","",Dr_Info_FIMM$CID[i][[1]]),"/synonyms/JSON"))$InformationList$Information$Synonym[[1]]))
    }, error = function(e) {
      Dr_Info_FIMM$compounds_[i][[1]]
  }))
  dataDTCdrug = dataDTC[toupper(dataDTC$compound_name) %in% synonyms_, ];  genes_DTCdrug = "";

  if(nrow(dataDTCdrug) > 0){
    dataDTCdrug = aggregate(standard_value~., dataDTCdrug, median)

    # keep only the most common bioactivity values
    dataDTCdrug = dataDTCdrug[dataDTCdrug$standard_type == names(sort(table(dataDTCdrug$standard_type),decreasing = !0))[[1]], ]
    genes_DTCdrug = paste0(unique(unlist(dataDTCdrug$gene_names[dataDTCdrug$standard_value < 1000])), collapse=",")

    if(length(genes_DTCdrug)==0) genes_DTCdrug = "";
  }

  paste0(genes_DTCdrug,collapse=",")
}, mc.cores = 1)


# SCTransformed scRNA-seq data from SEURAT pipeline here is the use of sc Rna Seq data
aml_data_Diagnosis = readRDS("aml_data_Diagnosis5238.RDS")

#
library(HGNChelper)
g_symb = checkGeneSymbols(rownames(aml_data_Diagnosis))
aml_data_Diagnosis = aml_data_Diagnosis[!is.na(g_symb$Suggested.Symbol), ]; g_symb = g_symb[!is.na(g_symb$Suggested.Symbol), ]
rownames(aml_data_Diagnosis) = g_symb$Suggested.Symbol; rownames_aml_data_Diagnosis = rownames(aml_data_Diagnosis)

#############
dss_aml1 = read.csv2("./DSRT_3_patient/DSRT_5238.tsv", sep="\t", stringsAsFactors = F)
dss_aml1 = dss_aml1[grepl("Live",dss_aml1$readout),]
dss_aml1$DSS = as.numeric(as.character(dss_aml1$DSS)) 

#dssamlcp = do.call("rbind", lapply(1:nrow(dss_aml1), function(i){
# data.frame(DRUG_NAME = dss_aml1$drug[i], CONCENTRATION = c(0.1,1,10,100,1000), SCREEN_NAME = "5750", PERCENT_INHIBITION = unlist(strsplit(dss_aml1$dose.response[i], "; ")))
#}))
#openxlsx::write.xlsx(dssamlcp,"forDSS_RF6333.xlsx")
dssamlcp2 = openxlsx::read.xlsx("./5238_DSRT_analysis_table_Rpipeline.xlsx")
dssamlcp2$EC50 = round(log10(dssamlcp2$EC50)+2)
dssamlcp2 = dssamlcp2[,c("DRUG_NAME","EC50")]
dss_aml1 = merge(dss_aml1, dssamlcp2, by.x = "drug", by.y = "DRUG_NAME")
dss_aml1$DSScp =dss_aml1$DSS

for(i in 1:nrow(dss_aml1)){
  dss_aml1$DSS[i] = 100 - as.numeric(as.character(unlist(strsplit(dss_aml1$dose.response[i], "; "))))[dss_aml1$EC50[i]]
}
print(paste0("Disease stage: ", unique(dss_aml1$disease.stage)))

# The result of the above code is a data frame named dss_aml1 that contains information about the drug sensitivity and dose response of different drugs for a patient with acute myeloid leukemia (AML). The data frame has 12 columns and 120 rows. The columns are:
#   
#   drug: The name of the drug tested on the patient’s cells.
# dose.response: A vector of percentage inhibition values for different concentrations of the drug, separated by semicolons.
# DSS: The drug sensitivity score, which is a measure of how effective the drug is at inhibiting the cell growth. It is calculated as 100 minus the percentage inhibition at the concentration that corresponds to the EC50 value.
# readout: The type of cell viability assay used to measure the cell growth. It can be either Live or Dead.
# disease.stage: The stage of the disease for the patient. It can be either Relapse or Remission.
# EC50: The half maximal effective concentration, which is the concentration of the drug that produces 50% of the maximum possible inhibition of the cell growth. It is transformed by applying the logarithm base 10 function and adding 2, then rounding to the nearest integer.
# DRUG_NAME: The same as the drug column, but with different capitalization. It is used to merge the data frames dss_aml1 and dssamlcp2.
# DSScp: A copy of the original DSS column before updating it with the new values.
# SCREEN_NAME: The name of the screen that the data belongs to. It is always 5750 in this case.
# CONCENTRATION: The concentration of the drug in nanomolar units. It can be one of the following values: 0.1, 1, 10, 100, or 1000.
# PERCENT_INHIBITION: The percentage of cell growth inhibition at the given concentration of the drug. It is a numeric value between 0 and 100.
# DSS2: The same as the DSS column, but with a different name. It is used to avoid duplicate column names when merging the data frames.
# The code also prints the disease stage of the data in dss_aml1, which is Relapse.

# # prepare for predictions
dss_aml1$drug = toupper(as.character(dss_aml1$drug)); Dr_Info_FIMM$compounds_ = toupper(as.character(Dr_Info_FIMM$compounds_));
dss_aml1 = merge(dss_aml1, Dr_Info_FIMM, by.x = "drug", by.y = "compounds_", all = F)
dss_aml1 = dss_aml1[dss_aml1$dubious == 0,]
dss_aml1 = dss_aml1[dss_aml1$DSS>=0,]

# gene sets
dss_aml1 = dss_aml1[dss_aml1[["TARGETS_JOINED2"]]!="",];
gs_ = mclapply(dss_aml1$drug, function(d_){ 
  genesInSelectedSets = as.character(dss_aml1[dss_aml1$drug == d_, "TARGETS_JOINED2"])
  
  # Find gene synonyms
  geneSynonymsMix = unique(na.omit(checkGeneSymbols(unique(unlist(strsplit(genesInSelectedSets,","))))$Suggested.Symbol))
  
  # Genes to keep
  GeneIndToKeep = rownames_aml_data_Diagnosis %in% geneSynonymsMix; 
  ids_ = rownames_aml_data_Diagnosis[GeneIndToKeep]; ids_
}, mc.cores = 20)
gc(T)

# exclude gene sets with genes not found in data
dss_aml1 = dss_aml1[!sapply(gs_, function(i) length(i)==0),]; gs_ = gs_[!sapply(gs_, function(i) length(i)==0)]
#waterfall::waterfallchart( drug ~ DSS, data = dss_aml1)
names(gs_) = as.character(dss_aml1$drug)


#FeaturePlot("")
# 
# ################################################## GENE sets for combis
# AllCombinations = expand.grid(names(gs_),names(gs_),stringsAsFactors = F)
# AllCombinations = AllCombinations[AllCombinations$Var1 != AllCombinations$Var2, ]
# 
# AllCombinations$merged = sapply(1:nrow(AllCombinations), function(i) paste0(sort(AllCombinations[i,]), collapse = ","))
# AllCombinations$Var1 = AllCombinations$Var2 = NULL
# AllCombinations_unique = unique(AllCombinations$merged)
# 
# # gs_combis = sapply(1:length(AllCombinations_unique), function(i_){
# #   drugs_ = unlist(strsplit(AllCombinations_unique[i_], ","))
# #   unique(unname(unlist(c(gs_[names(gs_) == drugs_[1]], gs_[names(gs_) == drugs_[2]]))))
# # })
# # 
# # names(gs_combis) = AllCombinations_unique

########################################################################
# GSVA
#gs_combis_all = c(gs_combis, gs_)

# library(SingleR)
# esrnaseqCombi = calculateSignatures(aml_data_Diagnosis, species = "Human", signatures = gs_combis_all,
#                     n.break = 1000, numCores = 6)

library(GSVA)
esrnaseqCombiSing <- gsva(aml_data_Diagnosis, gs_, parallel.sz=22, method = "gsva", kcdf = "Gaussian")
#   rownames(esrnaseqCombi) = names(gs_combis_all)
rownames(esrnaseqCombiSing) = names(gs_)
#esrnaseqCombiSingle = esrnaseqCombi[(length(gs_combis)+1):length(gs_combis_all), ]
###############################################################################################esrnaseqCombiSing <- readRDS("esrnaseqCombiSing_5238.RDS")
# 
# esrnaseq <- gsva(data.matrix(pbmc@data), gs_, parallel.sz=26)
# rownames(esrnaseq) = dss_aml1$drug

#cor_mat = cor(t(esrnaseqCombiSing))

# plot a heatmap
# esrnaseqcp= esrnaseq
# colnames(esrnaseqcp) = 1:ncol(esrnaseqcp)

# pbmcVAR2 = data.matrix(pbmc@data)[rownames(data.matrix(pbmc@data)) %in% pbmcVAR@var.genes,]
#esrnaseqCombiSing = esrnaseqCombi#[(length(gs_combis)+1):(length(gs_combis_all)),]
#esrnaseqCombiComb = esrnaseqCombi[1:length(gs_combis),]
#pdf("combi5249.pdf", width = 20, height = 50); pheatmap(esrnaseqCombiSing); dev.off()


# ###############################################
# .libPaths("/fs/projects/breeze/code/aleks/RNAseqDATAall/171020/libs")
# setwd("/fs/projects/breeze/code/aleks/scRNAseqGSVA")
# lapply(c("GSVA", "xgboost","parallel"), library, character.only = T)
# load("FIMMready_DGidb_DTCpotent_aml1_GSVA.RDATA")
# #############################################


processed_data = as.data.frame(esrnaseqCombiSing); processed_data$labeloutput = dss_aml1$DSS
#processed_dataTest = data.matrix(as.data.frame(esrnaseqCombiComb));
# dupl = which(duplicated(processed_data))
# 
# sstmp  = apply(processed_data, 1, function(r) all(r == processed_data[dupl[2],]))
# aaaaa = processed_data[sstmp, ]
# dss_aml1[dss_aml1$drug %in% rownames(aaaaa),]

library(ParamHelpers); library(xgboost)
### parameter set 
gdes <- function() generateDesign(n = 50, par.set = makeParamSet(
  makeNumericVectorParam("logNtree", len = 1, lower = 6, upper = 11), makeNumericVectorParam("lambda", len = 1, lower = 0, upper = 1),
  makeNumericVectorParam("alpha", len = 1, lower = 0, upper = 1), makeIntegerVectorParam("maxdepth", len = 1, lower = 3, upper = 10),
  makeNumericVectorParam("colsample_bytree", len = 1, lower = .75, upper = 1), makeNumericVectorParam("subsample", len = 1, lower = .75, upper = 1), 
  makeNumericVectorParam("eta", len = 1, lower = .01, upper = .5)
), fun = lhs::randomLHS)
des <- gdes()
CORvalgl <<- list()

### generate an objective function and start training
obj.fun = function(x) {
  
  logNtree = x[1]; lambda = x[2]; alpha = x[3]; maxdepth = x[4]; colsample_bytree = x[5]; subsample = x[6]; eta = x[7]; MAE_ <- 0; RMSE_ <- 0; 
  COR_ <- 0; COR2_ <- 0; CORval = as.data.frame(cbind(rep(NA, nrow(processed_data)),rep(NA, nrow(processed_data)),rep(NA, nrow(processed_data))))
  
  # repeated CV
  for(repCv in 1:6){ #saveRDS(Sys.time(),as.character(Sys.time()))
    
    flds1 <- caret::createFolds(1:nrow(processed_data), k = 5, list = T, returnTrain = F);
    
    for(k in 1:length(flds1)){
      testData <- processed_data[flds1[[k]],]; trainData <- processed_data[-flds1[[k]],];
      #trainData = rbind(trainData, rep(1,ncol(trainData))); trainData$labeloutput[nrow(trainData)] = 100
      
      
      fit = xgboost(data=data.matrix(trainData[, -which(names(trainData) == "labeloutput")]),label = trainData$labeloutput, verbose = F, nrounds=round(2**logNtree[[1]]), nthread = (parallel::detectCores()-4),
                    params=list(objective = "reg:linear", max.depth=maxdepth[[1]], eta=eta[[1]], lambda = lambda[[1]], 
                                alpha = alpha[[1]], colsample_bytree = colsample_bytree[[1]], subsample = subsample[[1]]))
      
      ypred = predict(fit, data.matrix(testData[, -which(names(testData) == "labeloutput")]));
      
      #F1_ <- F1_ + f1Score(actual = testData$labeloutput, predicted = ypred, cutoff = .5); AUC_ <- AUC_ + auc(actual = testData$labeloutput, predicted = ypred); 
      #MAE_ <- MAE_ + mae(actual = testData$labeloutput, predicted = ypred); RMSE_ <- RMSE_ + rmse(actual = testData$labeloutput, predicted = ypred);
      MAE_ <- MAE_ + ModelMetrics::mae(actual = testData$labeloutput, predicted = ypred);
      RMSE_ <- RMSE_ + ModelMetrics::rmse(actual = testData$labeloutput, predicted = ypred);
      COR_ <- COR_ + cor(ypred,testData$labeloutput)
      COR2_ <- COR2_ + cor(ypred,testData$labeloutput, method = "spearman")
      CORval[flds1[[k]],repCv] = ypred
    }
  }
  
  #paste0(c(F1_, AUC_, MAE_, RMSE_, logLOSS_), collapse = ",")
  CORvalgl <<- append(CORvalgl,list(CORval))
  paste0(c(MAE_, RMSE_, COR_, COR2_), collapse = ",")
}
gc(T); des$y = apply(des, 1, obj.fun)



##########################################################des = readRDS("des5238new.RDS")
##########################################################CORvalgl = readRDS("CORvalgl3258.RDS")
#saveRDS(CORvalgl, "CORvalgl5249.RDS")
#saveRDS(des,"des5249.RDS")

####################Add. func
# CORvalgl = c(readRDS("CORvalgl60.RDS"), readRDS("CORvalgl2.RDS"), readRDS("CORvalgl15.RDS"), readRDS("CORvalgl40.RDS"))
# des = rbind(readRDS("gdes60.RDS"), readRDS("gdes2.RDS"), readRDS("gdes15.RDS"), readRDS("gdes40.RDS"))
#CORvalgl = readRDS("CORvalgl5249.RDS")
#des = readRDS("des5249.RDS")

# # 
# # CORvalgl = c(readRDS("CORvalgl10fold.RDS"), readRDS("CORvalgl10fold_101.RDS"), readRDS("CORvalgl10fold_102.RDS"))
# # des =  rbind(readRDS("des10fold.RDS"), readRDS("des10fold_101.RDS"), readRDS("des10fold_102.RDS"))
# 
# 
cols_des = c()
for(i in 1:length(des$y)){ cols_des = rbind(cols_des, as.numeric(as.character(unlist(strsplit(des$y[i],","))))) }
cols_des = as.data.frame(cols_des); colnames(cols_des) = c("MAE","RMSE","PCC","SCC")
des = cbind(des, cols_des)


##########


CORvalglcp = CORvalgl[order(des$MAE)][1:3]
models = des[order(des$MAE),][1:3,]



# #
pred_CV = do.call("cbind",lapply(CORvalglcp, function(i){
  Biobase::rowMedians(data.matrix(i))
}))
# #
# plot(rowMeans(pred_CV), as.numeric(as.character(processed_data$labeloutput)))
#
# aa = as.data.frame((cbind(processed_data$labeloutput, pred_CV)))
# aa$c <- predict(prcomp(~V1+V2, aa,center=TRUE,scale.=TRUE))[,1]
# #
# # rot1.scaled <- scale(pca.object$rotation[,1])
# #
# # predict(prcomp(t(aa),center=TRUE,scale.=TRUE))[,1]
#
# #xx = data.frame(x = rowMeans(pred_CV), y = as.numeric(as.character(processed_data$labeloutput)), z = aaaa$scale )
# #ggplot(aa, aes(V1,V2,color = c))+geom_point()
#
# plot(rowMeans(pred_CV), as.numeric(as.character(processed_data$labeloutput)),col=aaaa$center)
#
# summary(pca.object)
# par(mfrow=c(1,2))
# plot(pca.object)
# biplot(pca.object)
#
# #
# # # conformal
err_ = abs(scales::rescale(Biobase::rowMedians(pred_CV),c(0,1)) - scales::rescale(as.numeric(as.character(processed_data$labeloutput)),c(0,1)))
err_2 = abs(Biobase::rowMedians(pred_CV) - as.numeric(as.character(processed_data$labeloutput)))#**6
processed_data_err = processed_data; processed_data_err$labeloutput = err_;
processed_data_err2 = processed_data; processed_data_err2$labeloutput = err_2;
#
fit_err = xgboost(data=data.matrix(processed_data_err[, -which(names(processed_data_err) == "labeloutput")]), nrounds = 500,
                  label = processed_data_err$labeloutput, verbose = F, nthread = (parallel::detectCores()-1))
fit_err2 = xgboost(data=data.matrix(processed_data_err2[, -which(names(processed_data_err2) == "labeloutput")]), nrounds = 500,
                   label = processed_data_err2$labeloutput, verbose = F, nthread = (parallel::detectCores()-1))

errs_pred = predict(fit_err, data.matrix(processed_data_err[, -which(names(processed_data_err) == "labeloutput")]))
errs_pred2 = predict(fit_err2, data.matrix(processed_data_err2[, -which(names(processed_data_err2) == "labeloutput")]))

aaa = data.frame(real = Biobase::rowMedians(pred_CV), pred = as.numeric(as.character(processed_data$labeloutput)),
                 errs_pred = errs_pred, d_ = rownames(processed_data), errs_pred2 = errs_pred2)
#
library(ggplot2)
ggplot(aaa, aes(pred, real, color = errs_pred)) + geom_point(size = 3) + theme_minimal()
#
# cc  = c();
#
# for(i in seq(0.9,0.999,.001)){
# aaa2 = aaa[aaa$errs_pred < quantile(aaa$errs_pred, i),]
# cc = c(cc, cor(aaa2$real,aaa2$pred))
# }
# #plot(seq(0.9,0.999,.001), cc)

# aaa = data.frame(pred1 = rowMeans(CORvalgl[[1]]), pred2 = rowMeans(CORvalgl[[2]]), pred3 = rowMeans(CORvalgl[[3]]), 
#                  real = as.numeric(as.character(processed_data$labeloutput)),  d_ = rownames(processed_data), 
#                  bel1 = (matrixStats::rowVars(data.matrix(CORvalgl[[1]])) < quantile(round(matrixStats::rowVars(data.matrix(CORvalgl[[1]]))), .75)),
#                  bel2 = (matrixStats::rowVars(data.matrix(CORvalgl[[2]])) < quantile(round(matrixStats::rowVars(data.matrix(CORvalgl[[2]]))), .75)),
#                  bel3 = (matrixStats::rowVars(data.matrix(CORvalgl[[3]])) < quantile(round(matrixStats::rowVars(data.matrix(CORvalgl[[3]]))), .75)))#,
#                  # bel4 = (matrixStats::rowVars(data.matrix(CORvalgl[[4]])) < quantile(round(matrixStats::rowVars(data.matrix(CORvalgl[[4]]))), .75)),
#                  # bel5 = (matrixStats::rowVars(data.matrix(CORvalgl[[5]])) < quantile(round(matrixStats::rowVars(data.matrix(CORvalgl[[5]]))), .75)),
#                  # bel6 = (matrixStats::rowVars(data.matrix(CORvalgl[[6]])) < quantile(round(matrixStats::rowVars(data.matrix(CORvalgl[[6]]))), .75)),
#                  # bel7 = (matrixStats::rowVars(data.matrix(CORvalgl[[7]])) < quantile(round(matrixStats::rowVars(data.matrix(CORvalgl[[7]]))), .75)))
# aaa$dif1 = aaa$pred1 - aaa$real; aaa$dif2 = aaa$pred2 - aaa$real; aaa$dif3 = aaa$pred3 - aaa$real;
#   


aaa$excl = 0;
#aaa$excl[aaa$errs_pred >= quantile(aaa$errs_pred,seq(0.9,0.999,.001)[which.max(cc)])] = 1
#aaa$excl[aaa$errs_pred2 >= quantile(aaa$errs_pred2,.7) ]=1#& aaa$errs_pred2 >= quantile(aaa$errs_pred2,.5)] = 1 #| aaa$errs_pred2 >= quantile(aaa$errs_pred2,.5) & aaa$real <5 | aaa$errs_pred2 >= quantile(aaa$errs_pred2,.5) & aaa$pred <5] = 1

aaa$excl[abs((aaa$real - aaa$pred))> des$MAE[order(des$MAE)[1]]/30] = 1 #5 fold 6 time repeated CV - 5*6 = 30
# cor(aaa$real,aaa$pred)
# cor(aaa2$real,aaa2$pred)
g_ = ggplot(aaa, aes(pred, real, color = as.factor(excl))) + geom_point(size = 3) + theme_minimal()
library(plotly)
ggplotly(g_)
# 
# 
cor(aaa$real,aaa$pred)
cor(aaa$real[aaa$excl==0],aaa$pred[aaa$excl==0])
bbb = aaa[aaa$excl==0,]


# 
# # cors_models = data.frame(real = processed_data$labeloutput,
# #   predicted = rowMeans(cbind(CORvalgl[[14]],CORvalgl[[282]],CORvalgl[[177]],CORvalgl[[282]],CORvalgl[[10]],CORvalgl[[172]],CORvalgl[[190]])))
# # 
# # ggplot(cors_models, aes(real, predicted, color = real)) +
# #   geom_point(shape = 16, size = 3, show.legend = F) +
# #   theme_minimal()
# # 
# # 
# # zz = c();
# # for(ll in 1:303)
# # zz = c(zz, cor(rowMeans(CORvalgl[[ll]]), processed_data$labeloutput))
# # # zz = c(zz, (cor(rowMeans(CORvalgl[[ll]])[processed_data$labeloutput>10 | rowMeans(CORvalgl[[ll]])>10], processed_data$labeloutput[processed_data$labeloutput>10 | rowMeans(CORvalgl[[ll]])>10])))
# # 
# # aa = data.frame(aa = rowMeans(data.matrix(cbind(CORvalgl[[5]],CORvalgl[[4]],CORvalgl[[3]],CORvalgl[[2]],CORvalgl[[1]]))), bb = processed_data$labeloutput)
# # cor(aa$aa,aa$bb)
# 


#dss_aml2 = dss_aml1[dss_aml1$DSS<15,]
dss_aml2 = merge(dss_aml1,bbb,by.x="drug",by.y="d_")

################################################## GENE sets for combis
AllCombinations = expand.grid(dss_aml2$drug,dss_aml2$drug,stringsAsFactors = F)
AllCombinations = AllCombinations[AllCombinations$Var1 != AllCombinations$Var2, ]

AllCombinations$merged = sapply(1:nrow(AllCombinations), function(i) paste0(sort(AllCombinations[i,]), collapse = ","))
AllCombinations$Var1 = AllCombinations$Var2 = NULL
AllCombinations_unique = unique(AllCombinations$merged)

gs_combis = sapply(1:length(AllCombinations_unique), function(i_){
  drugs_ = unlist(strsplit(AllCombinations_unique[i_], ","))
  unique(unname(unlist(c(gs_[names(gs_) == drugs_[1]], gs_[names(gs_) == drugs_[2]]))))
})

names(gs_combis) = AllCombinations_unique


#processed_data_wo_l = processed_data[, -which(names(processed_data) == "labeloutput")]
#processed_data_s = scale(processed_data)
#processed_data_s = scales::rescale(data.matrix(processed_data_wo_l), c(0,1))

processed_datat <- do.call("rbind", mclapply(names(gs_combis), function(i){
  apply(processed_data[unlist(strsplit(i, ",")),], 2, max)  
}, mc.cores = 16))



rownames(processed_datat) = names(gs_combis)

models_ <- lapply(1:3, function(i){
  fit_ = xgboost(data=data.matrix(processed_data[, -which(names(processed_data) == "labeloutput")]),label = processed_data$labeloutput, verbose = F, nrounds=round(2**models[i,"logNtree"]), nthread = (parallel::detectCores()-1),
                 params=list(objective = "reg:linear", max.depth=models[i,"maxdepth"], eta=models[i,"eta"], 
                             lambda = models[i,"lambda"], alpha = models[i,"alpha"], colsample_bytree = models[i,"colsample_bytree"]))
  fit_
})
gc(T)



importance_matrix <- xgb.importance(colnames(processed_data), model = models_[[1]])

xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")
#  
# importance <- xgb.importance(feature_names = names(processed_data)[names(processed_data) != "labeloutput"], model = fit)
# print(xgb.plot.importance(importance_matrix = importance, top_n = 15))
#   
processed_datat = as.data.frame(processed_datat); processed_datat$labeloutput = NULL; processed_datat = data.matrix(processed_datat)
ypredComb = predict(models_, processed_datat);
gc(T)
# fit_errComb = predict(fit_err, processed_dataTest)
# Combi5249 = data.frame(pred = rowMeans(ypredComb), err = fit_errComb, combis = rownames(processed_dataTest), stringsAsFactors = F)
# 
# Combi5249$excl = 0;
# #aaa$excl[aaa$errs_pred >= quantile(aaa$errs_pred,seq(0.9,0.999,.001)[which.max(cc)])] = 1
# Combi5249$excl[Combi5249$err >= quantile(Combi5249$err,.95)] = 1
# Combi5249Clean = Combi5249[Combi5249$excl == 0, ]
# 

Combi5249 = data.frame(pred = rowMeans(do.call("cbind",ypredComb)), combis = rownames(processed_datat), stringsAsFactors = F)
Combi5249Clean = Combi5249


aaa = openxlsx::read.xlsx("tsne_DF_5238_healthystatus.xlsx")
#aaa$Status = NA; aaa$Status[aaa$cluster==8]="H" #################################################################################################################
H_cells = colnames(processed_datat) %in% aaa$cellnames[!is.na(aaa$Status)]
processed_datatH = processed_datat[, H_cells]

Combi5249Clean$Helath = Biobase::rowMedians(processed_datatH)
Combi5249Clean$DSS_r = scales::rescale(rank(Combi5249Clean$pred), c(0,1))
Combi5249Clean$Health_r = scales::rescale(rank(Combi5249Clean$Helath), c(1,0))
#Combi5249Clean$rank = Combi5249Clean$DSS_r + Combi5249Clean$Health_r


Combi5249Clean$HSA = sapply(strsplit(Combi5249Clean$combis, "\\,"), function(i){
  max(dss_aml1[dss_aml1$drug == i[[1]], "DSS"], dss_aml1[dss_aml1$drug == i[[2]], "DSS"])
})
Combi5249Clean$DSS1 = sapply(strsplit(Combi5249Clean$combis, "\\,"), function(i){
  dss_aml1[dss_aml1$drug == i[[1]], "DSS"]
})
Combi5249Clean$DSS2 = sapply(strsplit(Combi5249Clean$combis, "\\,"), function(i){
  dss_aml1[dss_aml1$drug == i[[2]], "DSS"]
})

C_cells = !(colnames(processed_datat) %in% aaa$cellnames[!is.na(aaa$Health)])
processed_datatC = processed_datat[, C_cells]
Combi5249Clean$Can = rowSums(processed_datatC>0)
Combi5249Clean$Can_r = scales::rescale(rank(Combi5249Clean$Can), c(0,1))

Combi5249Clean$Synergy = Combi5249Clean$pred - Combi5249Clean$HSA
Combi5249Clean$Synergy_r = scales::rescale(rank(Combi5249Clean$Synergy), c(0,1))


#Combi5249Clean = Combi5249Clean[Combi5249Clean$DSS1 > 2 | Combi5249Clean$DSS2 > 2, ]

Combi5249Clean$rank = Combi5249Clean$DSS_r + (Combi5249Clean$Health_r + Combi5249Clean$Can_r)/2 + Combi5249Clean$Synergy_r 

Combi5249Clean = Combi5249Clean[Combi5249Clean$Synergy > quantile(Combi5249Clean$Synergy, .8),]#quantile(Combi5249Clean$Synergy, .75), 
Combi5249Clean = Combi5249Clean[Combi5249Clean$pred > quantile(Combi5249Clean$pred, .8), ]
Combi5249Clean = Combi5249Clean[Combi5249Clean$Helath < quantile(Combi5249Clean$Helath,.2),]

#Combi5249Clean = Combi5249Clean[Combi5249Clean$Can > 0 ,]



# 
# dss_aml1$Tox = as.numeric(as.character(dss_aml1$DSS)) - as.numeric(as.character(dss_aml1$sDSS))
# 
# Combi5249Clean$Toxicity = sapply(strsplit(Combi5249Clean$combis, "\\,"), function(i){
#   max(dss_aml1[dss_aml1$drug == i[[1]], "Tox"], dss_aml1[dss_aml1$drug == i[[2]], "Tox"])
# })
# 
# Combi5249Clean$rank = (scales::rescale(rank(Combi5249Clean$Synergy)) + scales::rescale(rank(Combi5249Clean$pred)) + (1-scales::rescale(rank(Combi5249Clean$Toxicity))))/3
# Combi5249Clean = Combi5249Clean[order(Combi5249Clean$rank, decreasing = T), ]
# 
# Combi5249Clean$ApprovedBoth = sapply(strsplit(Combi5249Clean$combis, "\\,"), function(i){
#   sum(c(dss_aml1[dss_aml1$drug == i[[1]], "development.phase"], dss_aml1[dss_aml1$drug == i[[2]], "development.phase"]) == "Approved")==2
# })
# 
# Combi5249CleanApp=Combi5249Clean[Combi5249Clean$ApprovedBoth == !0, ]

# Combi5249Clean$ApprovedAny = sapply(strsplit(Combi5249Clean$combis, "\\,"), function(i){
#   sum(c(dss_aml1[dss_aml1$drug == i[[1]], "development.phase"], dss_aml1[dss_aml1$drug == i[[2]], "development.phase"]) == "Approved")>0
# })
# Combi5249Clean = Combi5249Clean[Combi5249Clean$ApprovedAny == 1,]

#Combi5249Sort = Combi5249Clean[order(Combi5249Clean$Health_r, decreasing = T), ][1:100,]
#head(Combi5249Sort)
#processed_data$labeloutput = NULL

# Combi5249Clean = Combi5249Clean[Combi5249Clean$Can_r > quantile(Combi5249Clean$Can_r,.5),]



processed_data_s = data.matrix(processed_data[,colnames(processed_data) != "labeloutput"])

Combi5249Clean$diff <- do.call("rbind", mclapply(Combi5249Clean$combis, function(i){
  dt_ = unlist(strsplit(i, ","))
  sum((processed_data_s[dt_[1],] - processed_data_s[dt_[2],])>0.75)
}, mc.cores = 16))

Combi5249Clean$diff2 <- do.call("rbind", mclapply(Combi5249Clean$combis, function(i){
  dt_ = unlist(strsplit(i, ","))
  sum((processed_data_s[dt_[2],] - processed_data_s[dt_[1],])>0.75)
}, mc.cores = 16))

Combi5249Clean$diff1r = scales::rescale(rank(Combi5249Clean$diff), c(0,1))
Combi5249Clean$diff2r = scales::rescale(rank(Combi5249Clean$diff2), c(0,1))
Combi5249Clean$diffr = Combi5249Clean$diff1r + Combi5249Clean$diff2r



# Combi5249Sort = Combi5249Clean
# Combi5249Sort$DSS_r = scales::rescale(rank(Combi5249Sort$pred), c(0,1))
# Combi5249Sort$Synergy_r = scales::rescale(rank(Combi5249Sort$Synergy), c(0,1))
# Combi5249Sort$Health_r = scales::rescale(rank(Combi5249Sort$Helath), c(1,0))
# Combi5249Sort$rank = Combi5249Sort$DSS_r + Combi5249Sort$Synergy_r + Combi5249Sort$Health_r
Combi5249Sort = Combi5249Clean[order(Combi5249Clean$diffr, decreasing = !0), ]



Combi5249Sort$ApprovedBoth = sapply(strsplit(Combi5249Sort$combis, "\\,"), function(i){
  sum(c(dss_aml1[dss_aml1$drug == i[[1]], "development.phase"], dss_aml1[dss_aml1$drug == i[[2]], "development.phase"]) == "Approved")==2
})

Combi5249CleanApp=Combi5249Sort[Combi5249Sort$ApprovedBoth == !0, ]
# head(Combi5249Sort)
# 
# 
# Combi5249Sort = tidyr::separate(Combi5249Sort, combis, c("D1","D2"), sep = ",", remove = !1)
# 
# 

# 
# dt_ = data.frame(d_ = unique(c(Combi5249Sort$D1[1:15], Combi5249Sort$D2[1:15])), occ = 0)
# 
# # combis_ = c()
# # # for(i in 1:nrow(dt_)){
# # #   combis_ = c( combis_, head(which(Combi5249Sort$D1 == as.character(dt_$d_[i]) |  Combi5249Sort$D2 == as.character(dt_$d_[i])),1))
# # # }
# # 
# # combis_ = unique(combis_)
# # Combi5249Sort2 = Combi5249Sort[1,]; 
# # dt_$occ[dt_$d_==Combi5249Sort[1,"D1"]]=dt_$occ[dt_$d_==Combi5249Sort[1,"D1"]]+1
# # dt_$occ[dt_$d_==Combi5249Sort[1,"D2"]]=dt_$occ[dt_$d_==Combi5249Sort[1,"D2"]]+1
# # # 
#  for(i in 1:nrow(Combi5249Sort)){
#    if(dt_$occ[dt_$d_==Combi5249Sort[i,"D2"]] == 0 & dt_$occ[dt_$d_==Combi5249Sort[i,"D2"]] = 0){
#      dt_$occ[dt_$d_==Combi5249Sort[1,"D1"]]=dt_$occ[dt_$d_==Combi5249Sort[1,"D1"]]+1
#      dt_$occ[dt_$d_==Combi5249Sort[1,"D2"]]=dt_$occ[dt_$d_==Combi5249Sort[1,"D2"]]+1
#      Combi5249Sort2 = rbind(Combi5249Sort2, Combi5249Sort[i,])
#    }
#  }
# # 

# library(pheatmap)
# pred_1 = processed_data_s[rownames(processed_data_s) %in% unlist(strsplit(Combi5249Sort[1,"combis"], "\\,")), ]
# c_ = data.frame(H = as.character(H_cells)); rownames(c_) = colnames(pred_1)
# print(pheatmap(pred_1, annotation_col = c_))


# medians_ = sapply(1:1000, function(i){
#   pred_1 = processed_data_scaled[rownames(processed_data_scaled) %in% unlist(strsplit(Combi5249Sort[i,"combis"], "\\,")), ]
#   #pred_1 = rbind(scale(processed_data_scaled[1,]), scale(processed_data_scaled[2,]))
#   median((apply(pred_1, 2, max) - apply(pred_1, 2, min)) / apply(pred_1, 2, min))
# })
# hist(medians_,100)
# 
# medians_ = sapply(1:2000, function(i){
#   pred_1 = processed_data[rownames(processed_data) %in% unlist(strsplit(Combi5249Clean[i,"combis"], "\\,")), ]
#   sum((apply(pred_1, 2, max)/apply(pred_1, 2, min)) > 3)
# })
# hist(medians_,100)


#  medians_ = sapply(1:500, function(i){
#    pred_1 = processed_data_scaled[rownames(processed_data_scaled) %in% unlist(strsplit(Combi5249Sort[i,"combis"], "\\,")), ]
#    median((apply(pred_1, 2, max) - apply(pred_1, 2, min)) / apply(pred_1, 2, min))
#  })
#  hist(medians_,100)
# 
# # order(medians_, decreasing = T)[1:20]
#  #
#  Combi5249Sort[order(medians_, decreasing = T)[1:20],"combis"]
#  

# 

# 
# #processed_datatCp = processed_datatC; processed_datatCp[processed_datatCp>-0.5] = 1; processed_datatCp[processed_datatCp<0] = 0;
#ind_ =order(Biobase::rowMedians(processed_datatC), decreasing = T)#[1:10]
# 
# #$rank = Combi5249Sort$DSS_r + Combi5249Sort$Health_r
#Combi5249Sort = Combi5249Sort[order(Combi5249Sort$Canc_r, decreasing = T),]


for(i in 1:nrow(Combi5249Sort)){
  png(paste0("./5750/",i,".png"))
  pred_1 = processed_data_s[rownames(processed_data_s) %in% unlist(strsplit(Combi5249Sort[i,"combis"], "\\,")), ]
  library(pheatmap)
  
  c_ = data.frame(H = as.character(H_cells)); rownames(c_) = colnames(pred_1)
  print(pheatmap(pred_1, annotation_col = c_, cluster_rows = F, clustering_distance_cols = "minkowski",color= colorRampPalette(c("lightgrey","#ff0058"))(100) ))
  xdev.off()
}

(abs(abs(min(pred_1[2,]))-abs(min(pred_1[1,])))/(max(pred_1[2,]) - min(pred_1[2,])))*100
print(pheatmap(pred_1, annotation_col = c_, cluster_rows = F, clustering_distance_cols = "minkowski",
               color= c(rep("lightgrey",84),colorRampPalette(c("lightgrey","#520dff"))(1000-83)) ))



pdf("Vis5238_A_A_1.pdf", width=7, height = 6)
print(pheatmap(pred_1, annotation_col = c_, cluster_rows = F, clustering_distance_cols = "minkowski",color= colorRampPalette(c("lightgrey","#ff0058"))(100) ))
dev.off()

pdf("Vis5238_A_A_2.pdf", width=7, height = 6)



#### NEW FOR PLOTTING COMPLEMENT
library(RColorBrewer)
aaa = openxlsx::read.xlsx("tsne_DF_5238_healthystatus.xlsx")
H_cells = colnames(esrnaseqCombiSingCPCP) %in% aaa$cellnames[aaa$type=="CD8+ T cells"]
c_ = data.frame(H = as.character(H_cells)); rownames(c_) = colnames(esrnaseqCombiSingCPCP)
pdf("VisVen5238.pdf", width=7, height = 6)
print(pheatmap(esrnaseqCombiSingCPCP, annotation_col = c_,color= colorRampPalette(c("lightgrey","#ff0058")) ))
dev.off()

# 
#  pred_1 = processed_data[rownames(processed_data) %in% unlist(strsplit(Combi5249[1 ,"combis"], "\\,")), ]
#  #pred_1 = pred_1[,apply(pred_1, 2, max)<10]
#  library(pheatmap)
#  pheatmap(pred_1)
#  








# ypredComb$ApprovedBoth = sapply(strsplit(ypredComb$pair_names, "\\+"), function(i){
#   sum(c(dss_aml1[dss_aml1$drug == i[[1]], "development.phase"], dss_aml1[dss_aml1$drug == i[[2]], "development.phase"]) == "Approved")==2
# }) 
# 

#  
# 
# 
# ################################################## GENE sets for combis
# AllCombinations = expand.grid(names(gs_),names(gs_),stringsAsFactors = F)
# AllCombinations = AllCombinations[AllCombinations$Var1 != AllCombinations$Var2, ]
# AllCombinations$name = paste0(AllCombinations$Var1, "+", AllCombinations$Var2)
# 
# gs_combis = sapply(1:nrow(AllCombinations), function(row_){
#   unique(unname(unlist(c(gs_[names(gs_) == AllCombinations[row_,1]], gs_[names(gs_) == AllCombinations[row_,2]]))))
# })
# 
# names(gs_combis) = AllCombinations$name
# 
# esrnaseqCombi <- gsva(data.matrix(pbmc@data), gs_combis, parallel.sz=4)
# rownames(esrnaseqCombi) = names(gs_combis)
# 
# 
# ################################################## 
# 
# # genes_CD_list = list(genes_CD); names(genes_CD_list) = "CDs";
# # CD_data_GSVA <- gsva(data.matrix(pbmc@data), genes_CD_list, parallel.sz=4)
# CD_data = data.matrix(pbmc@data[rownames(pbmc@data) %in% unlist(genes_CD_list), ])
# CD_data = CD_data[, colSums(CD_data != 0) > 0]
# colnames(CD_data)
# 
# 
# ################################################## FIND best combis
# 
# esrnaseqCombi = readRDS("esrnaseqCombi.RDS")
# processed_data_test = as.data.frame(esrnaseqCombi[1:nrow(AllCombinations),]);
# processed_data = as.data.frame(esrnaseqCombi[(nrow(AllCombinations)+1):nrow(esrnaseqCombi),])
# rm(esrnaseqCombi); gc(T)
# processed_data_test = data.matrix(processed_data_test);  gc(T)
# 
# ypredComb = sapply(1:length(models_), function(i){
#   predict(models_[[i]], processed_data_test)
# })
# 
# ypredComb = data.frame(ypredComb = apply(ypredComb,1,median), stringsAsFactors = F); ypredComb$pair_names = unlist(rownames(processed_data_test))
# #ypredComb$pair_names = unlist(rownames(processed_data_test))
# #ypredComb2 = data.frame(ypredComb = ypredComb, pair_names = unlist(rownames(processed_data_test))); ypredComb$pair_names = as.character(ypredComb$pair_names)
# ypredComb$HSA = sapply(strsplit(ypredComb$pair_names, "\\+"), function(i){
#   max(dss_aml1[dss_aml1$drug == i[[1]], "DSS"], dss_aml1[dss_aml1$drug == i[[2]], "DSS"])
# })  
# ypredComb$DSS1 = sapply(strsplit(ypredComb$pair_names, "\\+"), function(i){
#   dss_aml1[dss_aml1$drug == i[[1]], "DSS"]
# }) 
# ypredComb$DSS2 = sapply(strsplit(ypredComb$pair_names, "\\+"), function(i){
#   dss_aml1[dss_aml1$drug == i[[2]], "DSS"]
# }) 
# 
# ypredComb$Synergy = ypredComb$ypredComb - ypredComb$HSA
# 
# ypredComb$ApprovedBoth = sapply(strsplit(ypredComb$pair_names, "\\+"), function(i){
#   sum(c(dss_aml1[dss_aml1$drug == i[[1]], "development.phase"], dss_aml1[dss_aml1$drug == i[[2]], "development.phase"]) == "Approved")==2
# }) 
# 
# ypredCombApproved = ypredComb[ypredComb$ApprovedBoth == !0,]
# 
# ypredComb = ypredComb[order(ypredComb$Synergy, ypredComb$ypredComb, decreasing = T),]
# ypredCombApproved = ypredCombApproved[order(ypredCombApproved$rank, decreasing = T),]
# 
# ypredComb = ypredComb[order(ypredComb$rank, decreasing = T),]
# 
# ypredComb = ypredComb[c(!0, !1),]
# ypredComb$rank_ef = rank(ypredComb$ypredComb) / nrow(ypredComb); ypredComb$rank_syn = rank(ypredComb$ypredComb) / nrow(ypredComb);
# #ypredComb$D2 = sapply(strsplit(ypredComb$pair_names, "\\+"), function(i) i[[2]]);
# ypredComb$rank = (ypredComb$rank_ef + ypredComb$rank_syn) / 2
# 
# #########################################################
# # CD markers


geneSynonymList = geneSynonym(unlist(strsplit(c("CD3","CD4","CD8","CD19","CD20"),",")), tax = species_$code_[species_$species_ == selectedSpecies_])
geneSynonymsMix = toupper(as.character(unlist(geneSynonymList)))
genes_CD = rownames(pbmc@data)[(toupper(rownames(pbmc@data)) %in% geneSynonymsMix)]

#pbmc@data[rownames(pbmc@data) %in% rownames(pbmc@data),]


################################################## SC3
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(pbmc),
    logcounts = log2(as.matrix(pbmcRealScales) + 1)
  )
)
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

num_clust = metadata(sc3_estimate_k(sce))$sc3$k_estimation
# acll = sc3(sce, ks=num_clust, n_cores=4)






########################
library(ggplot2)
des_ggplot = data.frame(x = rep(des$RMSE,6), y = c(des$logNtree,des$lambda,des$alpha,des$maxdepth,des$colsample_bytree,des$eta), 
                        z = c(rep("log N tree", nrow(des)), rep("alpha", nrow(des)), rep("lambda", nrow(des)), rep("maxdepth", nrow(des)), 
                              rep("colsample_bytree", nrow(des)), rep("learning rate", nrow(des))))
ggplot(des, aes(eta, RMSE)) + geom_point() + geom_smooth() + theme_minimal() 

library(plotly)
p <- plot_ly(x = des$logNtree, y = des$maxdepth, z = des$RMSE) %>% add_surface()




# 
# esrnaseqCombiSing =  esrnaseqCombiSing[, colnames(esrnaseqCombiSing) %in% gsub("_2","", colnames(pbmc5238@assays[["RNA"]]))]
## colnames(esrnaseqCombiSing) = paste0(colnames(esrnaseqCombiSing),"_2")
#pbmc5238@assays[["RNA"]] = pbmc5238@assays[["RNA"]][,colnames(esrnaseqCombiSing)]
#esrnaseqCombiSing = esrnaseqCombiSing[,colnames(pbmc5238@assays[["RNA"]])]
gg1= esrnaseqCombiSing[rownames(esrnaseqCombiSing)=="VENETOCLAX",]; gg1[gg1<0]=0;
pbmc6333_2 = AddMetaData(object = pbmc6333_2, metadata = gg1, col.name = "VENETOCLAX")
gg2= esrnaseqCombiSing[rownames(esrnaseqCombiSing)=="LOSMAPIMOD",]; gg2[gg2<0]=0;
pbmc6333_2 = AddMetaData(object = pbmc6333_2, metadata = gg2, col.name = "LOSMAPIMOD")
FeaturePlot(pbmc6333_2, features = c("VENETOCLAX", "LOSMAPIMOD"), blend = T, cols = c("blue", "red"))


gg1= esrnaseqCombiSing[rownames(esrnaseqCombiSing)=="ABEMACICLIB",]; gg1[gg1<0]=0;
pbmc6333_2 = AddMetaData(object = pbmc6333_2, metadata = gg1, col.name = "ABEMACICLIB")
gg2= esrnaseqCombiSing[rownames(esrnaseqCombiSing)=="AZD-5438",]; gg2[gg2<0]=0;
pbmc6333_2 = AddMetaData(object = pbmc6333_2, metadata = gg2, col.name = "AZD_5438")
FeaturePlot(pbmc6333_2, features = c("AZD_5438", "ABEMACICLIB"), blend = T, cols = c("blue", "red"))


gg1= esrnaseqCombiSing[rownames(esrnaseqCombiSing)=="NAVITOCLAX",]; gg1[gg1<0]=0;
pbmc6333_2 = AddMetaData(object = pbmc6333_2, metadata = gg1, col.name = "NAVITOCLAX")
gg2= esrnaseqCombiSing[rownames(esrnaseqCombiSing)=="TENIPOSIDE",]; gg2[gg2<0]=0;
pbmc6333_2 = AddMetaData(object = pbmc6333_2, metadata = gg2, col.name = "TENIPOSIDE")
FeaturePlot(pbmc6333_2, features = c("NAVITOCLAX", "TENIPOSIDE"), blend = T, cols = c("blue", "red"))







gg1= esrnaseqCombiSing5750[rownames(esrnaseqCombiSing5750)=="ABEMACICLIB",]; #gg1[gg1<0]=0;
pbmc5750 = AddMetaData(object = pbmc5750, metadata = gg1, col.name = "ABEMACICLIB")
gg2= esrnaseqCombiSing5750[rownames(esrnaseqCombiSing5750)=="AZD-5438",]; #gg2[gg2<0]=0;
pbmc5750 = AddMetaData(object = pbmc5750, metadata = gg2, col.name = "AZD-5438")
FeaturePlot(pbmc5750, features = c("ABEMACICLIB", "AZD-5438"), blend = T, cols = c("blue", "red"))


gg1= esrnaseqCombiSing5238[rownames(esrnaseqCombiSing5238)=="ABEMACICLIB",]; gg1[gg1<0]=0;
pbmc5238 = AddMetaData(object = pbmc5238, metadata = gg1, col.name = "ABEMACICLIB")
gg2= esrnaseqCombiSing5238[rownames(esrnaseqCombiSing5238)=="AZD-5438",]; gg2[gg2<0]=0;
pbmc5238 = AddMetaData(object = pbmc5238, metadata = gg2, col.name = "AZD_5438")
FeaturePlot(pbmc5238, features = c("ABEMACICLIB", "AZD_5438"), blend = T, cols = c("blue", "red"))



#pdf("ABEMACICLIB.pdf", width=4, height=3)
# FeaturePlot(pbmc5750, features = "markers_percent", cols = c("lightgrey", "blue"))
# dev.off()
# 
# pdf("AZD-5438.pdf", width=4, height=3)
# FeaturePlot(pbmc5750, features = "markers_percent2", cols = c("lightgrey", "red"))
# dev.off()

# 
# pdf("Background.pdf", width=4, height=3)
# FeaturePlot(pbmc5238, features = "markers_percent2", cols = c("lightgrey", "lightgrey"))
# dev.off()
# 