# Author: Vitalii Stebliankin (vsteb002@fiu.edu)

# Part of the code was taken from MetaboDiff tutorial:
#   https://rawgit.com/andreasmock/MetaboDiff/master/vignettes/MetaboDiff_tutorial.html

# Installation: 
# library('updateR')
# updateR()
# install.packages("WGCNA")
# library("devtools")
# install_github("andreasmock/MetaboDiff")

#library(MetaboDiff)

# Read input metabolomics matrix:
#   assay - a matrix containing the relative metabolic measurements
#   rowData - a dataframe containing the available metabolite annotation
#   colData - a dataframe containing sample metadata
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(ggplot2)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
source("RIO.R")


create_mae = function(assay,rowData,colData){
    rownames(colData) = colnames(assay)
    se = SummarizedExperiment(assays=as.matrix(assay),
                              rowData=rowData)
    experiment_list = list(raw=se)
    sampleMap = data.frame(primary=rownames(colData),
                           colname=colnames(se))
    sampleMap_list = listToMap(list(raw=sampleMap))
    met = MultiAssayExperiment(experiments = experiment_list,
                               colData = colData,
                               sampleMap = sampleMap_list)
    met
}
get_SMPDBanno <- function(met,
                          column_kegg_id,
                          column_hmdb_id,
                          column_chebi_id) {
    rowData = rowData(met[["raw"]])
    #db = read.csv(system.file("extdata", "metabolites.csv", package = "MetaboDiff"))  
    db = read.csv(paste(prefix(),"metabolites.csv",sep="/"))
    res = matrix(NA,nrow=nrow(rowData),ncol=ncol(db))
    colnames(res) = colnames(db)
    colnames(res) = paste0("SMPDB|",colnames(res))

    if(!is.na(column_kegg_id)){
        temp1 = as.matrix(db[match(rowData[,column_kegg_id],db$KEGG.ID),])
    }
    if(!is.na(column_hmdb_id)){
        temp2 = as.matrix(db[match(rowData[,column_hmdb_id],db$HMDB.ID),])
    }
    if(!is.na(column_chebi_id)){
        temp3 = as.matrix(db[match(rowData[,column_chebi_id],db$ChEBI.ID),])
    }
    if(exists("temp1")){
        res[is.na(res[,1]),] = temp1[is.na(res[,1]),]
    }
    if(exists("temp2")){
        res[is.na(res[,1]),] = temp2[is.na(res[,1]),]
    }
    if(exists("temp3")){
        res[is.na(res[,1]),] = temp3[is.na(res[,1]),]
    }
    rowData(met[["raw"]]) = data.frame(rowData,res)
    met
}

fill_experiments <- function(met) {
  # Fill normalized and imputed experiments.
  # Since we are using imputed and scaled data directly, raw, imputed and normalized will be the same set of data
  met_temp = met[["raw"]]
  sampleMap_imp = sampleMap(met)
  sampleMap_norm = sampleMap(met)
  sampleMap_imp_norm = sampleMap(met)
  
  levels(sampleMap_imp$assay) <- "imputed"
  levels(sampleMap_norm$assay) <- "norm"
  levels(sampleMap_imp_norm$assay) <- "norm_imputed"
  
  met2 = c(x=met,
           imputed=met_temp,
           sampleMap=sampleMap_imp)
  
  met3 = c(x=met2,
           norm=met_temp,
           sampleMap=sampleMap_norm)
  
  met4 = c(x=met3,
           norm_imputed=met_temp,
           sampleMap=sampleMap_imp_norm)
  
  met4
}
diff_test <- function(met, group_factors) {
    metadata(met) = vector("list",0)

    for (i in 1:length(group_factors)){
        if(length(levels(as.factor(colData(met)[[group_factors[i]]])))>2){
        coeff = sapply(1:nrow(assays(met)[["norm_imputed"]]),
                     function(x) aov(assays(met)[["norm_imputed"]][x,]~as.factor(colData(met)[[group_factors[i]]]))$coefficient)
        xlevels = unlist(aov(assays(met)[["norm_imputed"]][1,]~as.factor(colData(met)[[group_factors[i]]]))$xlevels)
        res = sapply(1:nrow(assays(met)[["norm_imputed"]]),
                     function(x) summary(aov(assays(met)[["norm_imputed"]][x,]~as.factor(colData(met)[[group_factors[i]]]))))
        res_df = data.frame(pval=as.vector(sapply(sapply(res,"[",i=5),"[",i=1)),
                            adj_pval=p.adjust(as.vector(sapply(sapply(res,"[",i=5),"[",i=1)),method = "fdr"),
                            dm=coeff[2,])
        metadata(met)[[paste0("anova_",group_factors[i],"_",paste(xlevels,collapse = "_vs_"))]] = res_df
        } else {
        xlev = levels(as.factor(colData(met)[[group_factors[i]]]))
        df = genefilter::rowttests(assays(met)[["norm_imputed"]],
                                   fac = as.factor(colData(met)[[group_factors[i]]]))
        res_df = data.frame(metabolite=rownames(df),
                            pval=df$p.value,
                            adj_pval=p.adjust(df$p.value,method="fdr"),
                            dm=df$dm)
        metadata(met)[[paste0("ttest_",group_factors[i],"_",xlev[2],"_vs_",xlev[1])]] = res_df
        }
    }

    met
}

volcano_plot <- function(met, group_factor, label_colors, dm_cutoff=0.5, p_adjust=TRUE, ...) {

    id = grep(group_factor,names(metadata(met)))[1]
    df=metadata(met)[[id]]
    name = names(metadata(met))[id]
    lv = levels(as.factor(colData(met)[[group_factor]]))
    xlabel = paste0("difference in means"," [",paste(lv,collapse = "-"),"]")

    if(p_adjust==TRUE){
        plot(df$dm, -1*log10(df$adj_pval),bty="n",pch=20,
             xlab=xlabel, ylab="adjusted -log10 p-value", ...)
        abline(h=-log10(0.05),lty=2)
        abline(v=dm_cutoff,lty=2)
        abline(v=-dm_cutoff,lty=2)
        Tsig = df$dm<(-dm_cutoff)&df$adj_pval<0.05
        points(x = df$dm[Tsig],y=-1*log10(df$adj_pval)[Tsig],col=label_colors[1],pch=20)
        Nsig = df$dm>(dm_cutoff)&df$adj_pval<0.05
        points(x = df$dm[Nsig],y=-1*log10(df$adj_pval)[Nsig],col=label_colors[2],pch=20)
    } else {
        plot(df$dm, -1*log10(df$pval),bty="n",pch=20,
             xlab=xlabel, ylab="-log10 p-value", ...)
        abline(h=-log10(0.05),lty=2)
        abline(v=dm_cutoff,lty=2)
        abline(v=-dm_cutoff,lty=2)
        Tsig = df$dm<(-dm_cutoff)&df$pval<0.05
        points(x = df$dm[Tsig],y=-1*log10(df$pval)[Tsig],col=label_colors[1],pch=20)
        Nsig = df$dm>(dm_cutoff)&df$pval<0.05
        points(x = df$dm[Nsig],y=-1*log10(df$pval)[Nsig],col=label_colors[2],pch=20)
    }
}
pca_plot <- function(met, group_factor, label_colors) {
    pca = prcomp(t(assay(met[["norm_imputed"]])))
    df = data.frame(pca$x)[,1:2]
    eigs <- pca$sdev^2
    vars = round(eigs[1:2] / sum(eigs),digits=2)*100
    df$grouping = as.vector(colData(met)[[group_factor]])
    p = ggplot(df, aes(x=PC1,y=PC2,colour=grouping)) +
        geom_point() + theme_classic() +
        scale_colour_manual(values=label_colors) +
        xlab(paste0("PC1 ","(",vars[1],"% of variance)")) +
        ylab(paste0("PC2 ","(",vars[2],"% of variance)"))
    return(p)
}



#-----------------------------------------------------
# READ METABOLITES ABUNDANCE VALUES FROM METABOLON INC (SCALED)
#-----------------------------------------------------

#curr_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
#`setwd(curr_dir)
input <- function(inputfile) {
	pfix <- prefix()
	parameters <- readParameters(inputfile)
assayT <<- read.table(paste(pfix, parameters["assay", 2], sep="/"), header = TRUE, sep=',')
rownames(assayT) <<- assayT$CHEM_ID
assayT$CHEM_ID <<- NULL
colData2 <<- read.table(paste(pfix, parameters["metadata", 2], sep="/"), header = TRUE, sep=',')
rownames(colData2) <<- colData2$PARENT_SAMPLE_NAME
colData2$PARENT_SAMPLE_NAME <<- NULL
colData2 <<- as.data.frame(colData2)

rowData2 <<- read.table(paste(pfix, parameters["annotations", 2], sep="/"), header=TRUE, sep='\t')
rownames(rowData2) <<- rowData2$CHEM_ID
rowData2$CHEM_ID <<- NULL
rowData2 <<- as.data.frame(rowData2)

}

run <- function() {
	met <<- create_mae(assayT,rowData2,colData2)


# Get annotations
met <<- get_SMPDBanno(met,
                     column_kegg_id=4,
                     column_hmdb_id=5,
                     column_chebi_id=NA)


# The data is already imputed and scaled
met <<- fill_experiments(met)
}

output <- function(outputfile) {
#-----------------------------------------------------
# PLOT DISTRIBUTION VALUES FROM METABOLON INC SCALING
#-----------------------------------------------------
group_factor="GROUP_NAME"
mdata = as.data.frame(longFormat(met[,,1],colDataCols=group_factor))
#mdata$value = log2(mdata$value)
plot1 = ggplot(mdata,
               mapping=aes(x=colname,
                           y=value,
                           fill=get(group_factor))) +
  geom_boxplot(lwd=0.2) + xlab("") + ggtitle("raw") +
  theme(axis.text.x=element_blank()) + ylab("Scaled abundance") 
  #scale_fill_manual(values=label_colors, name="group factor")
plot1

png(paste(outputfile, "quality_control.png", sep="/"))
#quartz.save('figures/quality_control.png', type = "png", device = dev.cur(), dpi = 300)

#-----------------------------------------------------
# VOLCANO PLOT
#-----------------------------------------------------
# Hypothesis testing
met = diff_test(met,
                group_factors = c(group_factor))
str(metadata(met), max.level=2)

# Volcano Plot
par(mfrow=c(1,2))
volcano_plot(met, 
             group_factor=group_factor,
             label_colors=c("darkseagreen","dodgerblue"),
             dm_cutoff=0.5,
             p_adjust = FALSE)

png(paste(outputfile, "volcano_fib4.png", sep="/"))
#quartz.save('figures/scaled_metabolon/volcano_fib4.png', type = "png", device = dev.cur(), dpi = 300)


# PCA and TSNE
source("http://peterhaschke.com/Code/multiplot.R")
pca_plot(met,
         group_factor=group_factor,
         label_colors=c("orange","darkseagreen","dodgerblue","red"))
#quartz.save('figures/pca.png', type = "png", device = dev.cur(), dpi = 300)
png(paste(outputfile, "pca.png", sep="/"))
}
