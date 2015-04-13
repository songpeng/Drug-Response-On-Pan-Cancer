# Drug Response on CNV.
# Songpeng Zu
# 2015-03-21

#-- NOTE
# We ignore the stages of the patients in our analysis.
# We assume that the genetic informatin in the patients keep stable on
# average during the drug treatment.

#-- Load the library
library(data.table)
library(stringr)
library(ggplot2)
#-- Load the patient information.
patient_infor <-
    read.table('cancer.patient.drug.response.gistic2.20150330.txt',quote="",
               colClasses=c("factor","character","factor","factor","factor"),
               sep="\t")
colnames(patient_infor) <- c("cancerType","patientID","drugName","response",
                              "dataExist")
#rownames(patient_infor) <- patient_infor$patientID
# Change the levels of patients' responses'
levels(patient_infor$response) <- c("Stable","Effective","Effective","Stable")
delrownan <- is.nan(patient_infor$response)
patient_infor[-delrownan,]
# Get the drug treatment data of Cisplatin on Copy Number Variation (CNV)
CisplatinCNV <- fread("Cisplatin.gistic2.gdac_20141206.txt")

#-- Analyse CNV per Cancer Type

#-- Analyse CNV in all the Cancer Type
# Match name from drug treatment file.
# Note: for an example, HNSA.TCGA-2G-AAFM-01A, TCGA-2G-AAFM is the patient ID in
# patien information file. 01A means "Primary Tumor", and "HNSA" means the tumor
# Type.

patient_cisplatin <- patient_infor[which(patient_infor$drugName == "Cisplatin"),]
rownames(patient_cisplatin) <- patient_cisplatin$patientID

# Statistics of responses in each tumor type.
# In this step, we aim to choose several tumors for single task learning.
responseEachTumorCis <- function(cancer_type){
    table(patient_cisplatin$response[which(patient_cisplatin$cancerType==cancer_type)])
}


# Transform the Samples' names in CisplatinCNV to the ones in patient_infor
nametrans <- function(name_drugtreat){
   m = str_match_all(name_drugtreat, "\\.(.+)-01A")
   return(m[[1]][2]) # Note: NA will be returned if no matching.
}

sampleCistreat <- colnames(CisplatinCNV)[2:ncol(CisplatinCNV)]
patientIDinCis <- unlist(lapply(sampleCistreat,function(x) nametrans(x)))
# Delete NA
tmpisna <- is.na(patientIDinCis)
patientIDinCis <- patientIDinCis[!tmpisna]

patient_response_cis <- patient_cisplatin[patientIDinCis,"response"]
responseCisReorder <- as.matrix(CisplatinCNV[1:nrow(CisplatinCNV),
                                             2:ncol(CisplatinCNV),with=FALSE])
responseCisReorder <- responseCisReorder[,-which(tmpisna==TRUE)]
row.names(responseCisReorder) <- unlist(CisplatinCNV[1:nrow(CisplatinCNV),1,
                                                     with=FALSE])

zeronumcutoff <- 0.05
samplenm <- length(patientIDinCis)
pairdetectSenseGene <- function(genenm){
    # Filter genes with too many zeros.
    tmparray <- responseCisReorder[genenm,]
    numequalzero <- sum(tmparray==0)
    if (numequalzero/samplenm <= zeronumcutoff){return(1)}
    else{
        # Get the p-value for one gene.
        pvalue <- pairwise.t.test(tmparray,patient_response_cis)$p.value
        return(pvalue)
    }
}

pvaluearrayCisCNV <- unlist(lapply(rownames(responseCisReorder),function(x)
                            pairdetectSenseGene(x)))
summaryCisCNV <- matrix(pvaluearrayCisCNV,
                           dimnames=list(
                               rownames(responseCisReorder),c("pvalue")))
# False Discovery Control
# pvalueajustCisCNV <- p.adjust(summaryCisCNV,"fdr")
# The mim(pvalueajustCisCNV) is 0.26.
pvaluesrtCisCNV <- sort(summaryCisCNV,index.return = TRUE)
pvaluesumCisCNV <- as.matrix(summaryCisCNV[pvaluesrtCisCNV$ix,])
row.names(pvaluesumCisCNV) <- row.names(summaryCisCNV)[pvaluesrtCisCNV$ix]

# Generate the boxplot
cnvboxplotcis <- function(genenm){
    dat <- data.frame(CopyNumVar=responseCisReorder[genenm,],
                      Treatment=patient_response_cis)
    p <- ggplot(dat,aes(Treatment,CopyNumVar))
    p + geom_boxplot()
}
