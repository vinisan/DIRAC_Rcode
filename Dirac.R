library(dplyr)
library(ggplot2)
library(doParallel)
library(pheatmap)




# this function reads the pathway  file, removes all pathways which has less than 3 genes and stored it as a list
readPathway<-function(path.File, exp.data){
    pathways<-readLines(path.File)
    path<-strsplit(pathways, "\t")
    
    return.pathways = list()
    for ( i in 1:length(path)){
        namesDIR<-as.matrix(path[[i]])
       
        if (length(path[[i]]) > 5 && length(which(rownames(exp.data) %in% namesDIR[,1])) > 5){
            return.pathways <-append(return.pathways, list(path[[i]]))
        }
    }
    
    names(return.pathways)<-sapply(return.pathways,'[',1)
    final.pathways = lapply(return.pathways,'[',-(1:2))
    return(final.pathways)
}





# does the pairwise comparison for the pathway; takes the rank ordered matrix for the pathway and samples of one class;
# for each sample in the class, counterrankwise pair comparison is performed through "rankvector" function and then appended to the class.Matrix,
# which represents the pairwise comparison for the whole class in which each column is the pairwise comparison for each sample.
doPairWise<-function(pathwayOrderMatrix, gene.pairs){
    
    pathCond.cols<-nrow(pathwayOrderMatrix)
    class.Matrix<-matrix(, nrow = gene.pairs)
   
    for ( sample.Count in 1:ncol(pathwayOrderMatrix)){
        sample<-vector('numeric', gene.pairs)
        pathCond<-as.matrix(pathwayOrderMatrix[,sample.Count])
        sample<-rankVector(pathCond, pathCond.cols, gene.pairs)
        class.Matrix<-cbind(class.Matrix,sample)
    }
    
    class.Matrix<-class.Matrix[,-1]
    mode.Vector<-vector('numeric',nrow(class.Matrix))
    nSamples <- ncol(class.Matrix)
        
    for ( gene.path in 1:nrow(class.Matrix) ){
        mode.gene = 0
        mode.gene = sum(class.Matrix[gene.path,])/nSamples
        
        if(mode.gene < 0.5){
            mode.Vector[gene.path] = 0
        }
        
        else{
            mode.Vector[gene.path] = 1
        }
        
    }
    return(list(class.Matrix, mode.Vector))
}


#pairwise comparison of the gene ranks generates counterbinary vector which is returned as an output to the main class vector which isreturned to the doPairWise function
rankVector<-function(pathCond, pathCond.cols, gene.pairs){
    counter= 0
   
    rank.compared<-vector('numeric', gene.pairs)
    for(i in 1:pathCond.cols){
        
        for (j in i:pathCond.cols){
            
            if ( i != j){
                counter= counter+1
                
                if(pathCond[i,1] >= pathCond[j,1]){
                   rank.compared[counter] = 1
                    
                }
                
                else{
                    rank.compared[counter] = 0
                }
                
            }
        }
    }
    
    return(rank.compared)
    
}

#this function generates the rank conservation score by comparing the expression matrix (ranked) with the template for that class.

calcRankMatching<-function(expr.Matrix,template){
    
    pairs<-nrow(expr.Matrix)
    rank.conservation<-vector("numeric", ncol(expr.Matrix))
    
    for ( columns in 1:ncol(expr.Matrix)){
        
        different<-which(expr.Matrix[,columns] != template)
        percentage.different = 1-(length(different)/pairs)
        rank.conservation[columns] = percentage.different
        
    }
    return(rank.conservation)
}


#Calculate the average accuracy using rank.difference matrix returns the average accuracy of the pathway as a classifier: James says "if i calculate accuracy (eta) as the average of sensitivity and specificity; sensitivity is TP / N1 and specificity is TN / N2. when i sum up TP and TN, i count 0.5 for ties" Average  = (TP/N1)*0.5 + (TN/N2)*0.5

calculateAccuracy<-function(rank.difference,nCond1,nCond2){
    
    true.positives = 0
    ties = 0
    specificity = 0
    sensitivity = 0
    accuracy = 0
    true.positives <- length(which(rank.difference[1:nCond1] >= 0))
    sensitivity = (true.positives)/nCond1
    true.negatives = 0
    true.negatives <- length(which(rank.difference[(nCond1 +1):(nCond1 + nCond2)] < 0))
    ties = length(which(rank.difference[(nCond1 +1):(nCond1 + nCond2)] == 0))
    specificity = (true.negatives )/nCond2
    accuracy = (specificity * 0.5) + (sensitivity * 0.5)
    
    return(accuracy)
}





#This function take the expression data, pathway list, indices for two conditions and a vector for pathway difference an inputs.
# for each pathway, data is subset and then for pathways with more than 3 genes, expression of genes are rank ordered and and after paiwise comparison
# through doPairwise function the mean changes in the ranks are calculated and difference between the pathway means are stored in pathway.Difference

runDirac<-function(data.Cond1, data.Cond2, pathway.List, pathway.accuracy){
    
    
    for (pathwayN in 1:length(pathway.List)){
        namesDIR<-as.matrix(pathway.List[[pathwayN]])
        
        pathwayNdata.Cond1<-data.Cond1[which(rownames(data.Cond1) %in% namesDIR[,1]),]
        pathwayNdata.Cond2<-data.Cond2[which(rownames(data.Cond2) %in% namesDIR[,1]),]
        
        pathwayNDataCond1.order = apply(pathwayNdata.Cond1,2, rank)
        pathwayNDataCond2.order = apply(pathwayNdata.Cond2,2, rank)
        
        pathCond.cols<-nrow(pathwayNDataCond1.order)
        gene.pairs <-(pathCond.cols*(pathCond.cols-1))/2
            
        cond1.matrix<-matrix(, nrow = gene.pairs)
        cond2.matrix<-matrix(, nrow = gene.pairs)
        
        cond1.list<-doPairWise(pathwayNDataCond1.order,gene.pairs)
        cond2.list<-doPairWise(pathwayNDataCond2.order,gene.pairs)
        
        total.Matrix<-cbind(cond1.list[[1]],cond2.list[[1]])
        rank.matching1<-calcRankMatching(total.Matrix,cond1.list[[2]])
        rank.matching2<-calcRankMatching(total.Matrix,cond2.list[[2]])
        rank.difference = rank.matching1 - rank.matching2
        
        nCond1 = ncol(data.Cond1)
        nCond2 = ncol(data.Cond2)
        pathway.accuracy[pathwayN] = calculateAccuracy(rank.difference,nCond1,nCond2)
    }
    return(pathway.accuracy)
}

#this function performs permutations and finds the if the permuted mean difference in more than pathway difference.
#each time permuted mean difference (permuted.Difference) is higher than pathway difference(pathway.Difference) "1" is added to a vector named
#"calculated.Difference". run.Dirac is run with permuted data and mean difference between the two conditions is computed for each pathway.



permuteDirac<-function(dirac.data,pathway.List, nCond1, nCond2, pathway.Difference, calculated.Accuracy){
    
    for (run in 1:permutations){
        
        print (run)
        permuted.Data<-dirac.data[,sample(ncol(dirac.data), replace = TRUE)]
        data.Cond1<-permuted.Data[,1:nCond1]
        data.Cond2<-permuted.Data[,(nCond1 +1):(nCond1 + nCond2)]
        
        pathway.Accuracy<-vector('numeric',length(pathway.List))
        permuted.Difference <- runDirac(data.Cond1, data.Cond2,pathway.List, pathway.Accuracy)
        more.Difference<-which(permuted.Difference >= pathway.Difference)
        
        if (length(more.Difference) > 0)
        {
            calculated.Accuracy[more.Difference]<-(calculated.Accuracy[more.Difference] + 1)   
        }
    }
    return(calculated.Accuracy)
}





# Cross Validation using leave one out crossValidation
doCrossValidation<-function(cond1.data,cond2.data,CV.pathwayList){
    
    CVaccuracy<-vector('numeric',length(CV.pathwayList) )
    index<-vector('character',length = (ncol(cond1.data) + ncol(cond2.data)))
    index[1:ncol(cond1.data)] = "1"
    index[(ncol(cond1.data)+1):(ncol(cond1.data) + ncol(cond2.data))] = "0"
    
    CVaccuracy[1:length(CV.pathwayList)] = 0
    rank.difference<-vector('numeric',length(CV.pathwayList))
    dirac.data<-cbind(cond1.data,cond2.data)
    
    for (samples in 1:ncol(dirac.data)){
        
        pathway.accuracy<-vector('numeric',length(CV.pathwayList))
        cvData<-dirac.data[,-samples]
        cvIndex<-index[-samples]
        holdData<-as.matrix(dirac.data[,samples])
        holdIndex<-index[samples]
        rownames(holdData)<-rownames(dirac.data)
        
        data1<-which(cvIndex == "1")
        data2<-which(cvIndex == "0")
       
        data.Cond1<-cvData[,data1]
        data.Cond2<-cvData[,data2]
        cvPathway.Accuracy<-vector('numeric',length(CV.pathwayList))
        
        nCond1<-ncol(data.Cond1)
        nCond2<-ncol(data.Cond2)
        
        for (pathwayN in 1:length(CV.pathwayList)){
            
            namesDIR<-as.matrix(CV.pathwayList[[pathwayN]])    
            pathwayNdata.Cond1<-data.Cond1[which(rownames(data.Cond1) %in% namesDIR[,1]),]
            pathwayNdata.Cond2<-data.Cond2[which(rownames(data.Cond2) %in% namesDIR[,1]),]
            
            holdDataGenes<-holdData[which(rownames(holdData) %in% namesDIR[,1]),]
            
            pathwayNDataCond1.order = apply(pathwayNdata.Cond1,2, rank)
            pathwayNDataCond2.order = apply(pathwayNdata.Cond2,2, rank)
            holdDataGenes.order     = as.matrix(rank(holdDataGenes))
            
            pathCond.cols<-nrow(pathwayNDataCond1.order)
            gene.pairs <-(pathCond.cols*(pathCond.cols-1))/2
            
            cond1.matrix<-matrix(, nrow = gene.pairs)
            cond2.matrix<-matrix(, nrow = gene.pairs)
            hold.matrix<-matrix(, nrow = gene.pairs)
            
            cond1.list<-doPairWise(pathwayNDataCond1.order,gene.pairs)
            cond2.list<-doPairWise(pathwayNDataCond2.order,gene.pairs)
            hold.list<-rankVector(holdDataGenes.order, nrow(holdDataGenes.order), gene.pairs)
            
            
            total.Matrix<-cbind(cond1.list[[1]],cond2.list[[1]])
            template1<-cond1.list[[2]]
            template2<-cond2.list[[2]]
            
            
            different1<-which(hold.list != template1)
            percentage.different1 = 1-(length(different1)/gene.pairs)
            rank.conservation1 = percentage.different1
            
            different2<-which(hold.list != template2)
            percentage.different2= 1-(length(different2)/gene.pairs)
            rank.conservation2 = percentage.different2
            
            rank.difference = rank.conservation1 - rank.conservation2
            if(rank.difference > 0 & holdIndex == 1)
            {
                CVaccuracy[pathwayN] = CVaccuracy[pathwayN] +1
            }
            else if(rank.difference < 0 & holdIndex == 0)
            {
                CVaccuracy[pathwayN] = CVaccuracy[pathwayN] +1
            }
            
            
        }
        
    }
    return((CVaccuracy/(ncol(dirac.data))))
}


#Dirac permutations using parallelization
doPermutations<-function(expr.data, pathway.list, nCond1, nCond2,  PathwayAccuracy, calculated.Accuracy, CVaccuracy, permutations, cores){
    
    calculated.Accuracy<-vector('numeric',length(pathway.list))
    registerDoParallel(cores)
    
    parallelPermuted <- foreach(i=1:cores, .combine='cbind') %dopar% permuteDirac(expr.data, pathway.list, nCond1, nCond2,  PathwayAccuracy, calculated.Accuracy)
    PermutedAccuracy =apply(parallelPermuted, 1, sum)
    Pvalue<-PermutedAccuracy/(permutations*16)
    q<-p.adjust(Pvalue,method = "BH", n = length(Pvalue))
    
    
    #Storing the values in a data matrix
    results<-cbind(names(pathway.list), PathwayAccuracy)
    results<-cbind(results, as.numeric(Pvalue))
    results<-cbind(results, as.numeric(q))
    results<-cbind(results, as.numeric(CVaccuracy))
    
    return(results)
}










#Main run

TindexDir<-read.csv("AggHumanDataforEVA.csv", nrow = 1, header = FALSE,row.names = 1)
Texpr.data<-read.csv("AggHumanDataforEVA.csv", skip = 1, header = TRUE, row.names = 1)


cond1.DIR<-which(TindexDir == "T9")
cond2.DIR<-which(TindexDir == "C9")

print(cond1.DIR)
print(cond2.DIR)
Tdata.Cond1<-NULL
Tdata.Cond2<-NULL
Tdata.Cond1<-Texpr.data[,cond1.DIR]
Tdata.Cond2<-Texpr.data[,cond2.DIR]
Tdirac.data<-cbind(Tdata.Cond1,Tdata.Cond2)
dim(Tdirac.data)
nCond1<-ncol(Tdata.Cond1)
nCond2<-ncol(Tdata.Cond2)




#path.File<-"c2.biocarta.v2.5.symbols.gmt"
path.File<-"c2.cp.kegg.v5.0.symbols.gmt"
Tpathway.list<-readPathway(path.File, Tdirac.data)
print(length(Tpathway.list))

print(length(Tpathway.list))
Tpathway.Accuracy<-vector('numeric',length(Tpathway.list))

## this runs the dirac to calculate the differences betweed two classes for each pathway
TPathwayAccuracy<-runDirac(Tdata.Cond1, Tdata.Cond2, Tpathway.list,  Tpathway.Accuracy)
TCVaccuracy<-doCrossValidation(Tdata.Cond1, Tdata.Cond2,Tpathway.list)

#Permutations with parallelization using 16 cores
set.seed(2334)
permutations = 1
cores= 16

storeResultsT<-doPermutations(Tdirac.data, Tpathway.list, nCond1, nCond2,  TPathwayAccuracy, Tcalculated.Accuracy, TCVaccuracy, permutations, cores)
hist(storeResultsT[,4], breaks = 10)

#Proteomic DIRAC analysis
PindexDir<-read.csv("AHumanProteeome.csv", nrow = 1, header = FALSE, row.names = 1)
Pexpr.data<-read.csv("AHumanProteeome.csv", skip = 1, header = FALSE, row.names = 1)

Pcond1.DIR<-which(PindexDir == "9T")
Pcond2.DIR<-which(PindexDir == "9C")

print(cond1.DIR)
print(cond2.DIR)
Pdata.Cond1<-NULL
Pdata.Cond2<-NULL
Pdata.Cond1<-Pexpr.data[,Pcond1.DIR]

Pdata.Cond2<-Pexpr.data[,Pcond2.DIR]
Pdirac.data<-cbind(Pdata.Cond1,Pdata.Cond2)
dim(Pdirac.data)
nCond1<-ncol(Pdata.Cond1)
nCond2<-ncol(Pdata.Cond2)



Ppathway.list<-readPathway(path.File, Pdirac.data)

print(length(Ppathway.list))
Ppathway.Accuracy<-vector('numeric',length(Ppathway.list))

## this runs the dirac to calculate the differences betweed two classes for each pathway
PPathwayAccuracy<-runDirac(Pdata.Cond1, Pdata.Cond2, Ppathway.list,  Ppathway.Accuracy)
PCVaccuracy<-doCrossValidation(Pdata.Cond1, Pdata.Cond2,Ppathway.list)

#Permutations
Pcalculated.Accuracy<-vector('numeric',length(Ppathway.list))
set.seed(2334)
storeResultsP<-doPermutations(Pdirac.data, Ppathway.list, nCond1, nCond2,  PPathwayAccuracy, Pcalculated.Accuracy, PCVaccuracy, permutations, cores)
hist(storeResultsP[,4], breaks = 10)


