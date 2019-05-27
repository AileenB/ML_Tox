require(tools)
require(xlsx)
require(class)
require(randomForest)

ML_Tox <- function(file, sheetInd=1, classLabelPos="active", classLabelNeg="passive", pca=TRUE, knn=TRUE, k=1, usePC=c(1:2), rf=TRUE, ntree=5000, rfe=TRUE, importance="MDA", setSeed=T, seedNumber=123){
  
  if (setSeed == T){
    set.seed(seedNumber)
  }
  # read input file

  if(file_ext(file) == "xlsx"){
    myData <- read.xlsx2(file, sheetInd)
    myData <- myData[which(myData[,1] != ""),]
    myData <- myData[,c(1:which(colnames(myData) == "Toxicity"))]
  } else if (file_ext(file) == "csv"){
    myData <- read.csv(file)
    myData <- myData[which(myData[,1] != ""),]
    myData <- myData[,c(1:which(colnames(myData) == "Toxicity"))]
  } else{
    print("Error: Please provide input as .xlsx or .csv")
  }
  
  myDescriptors <- myData[-c(1,dim(myData)[2])]
  rownames(myDescriptors) <- myData[,1]
  myDescriptors_numeric <- apply(myDescriptors,2,as.numeric)
  myTox <- myData[,dim(myData)[2]]
  myTox <- droplevels(myTox,"")
  
  # perform pca
  
  if(pca==T){
    myPCAResult <- prcomp(myDescriptors_numeric, scale = T)
  }
  
  # perform knn
  
  if(knn==T){
    
    myChosenPCs <- myPCAResult$x[,usePC]
    
    predictKNN <- function(i){
      myKNN <- knn(myChosenPCs[-i,],myChosenPCs[i,], cl=myTox[-i], k=k, prob=T)
      return(myKNN)
    }
  
    myPredictions <- sapply(c(1:dim(myChosenPCs)[1]),predictKNN)
  
    myComparisonOfClasses <- data.frame(PredictedClass=myPredictions,TrueClass=myTox)
    rownames(myComparisonOfClasses) <- rownames(myDescriptors)
    print("PCA + knn:")
    print(myComparisonOfClasses)
    
    myCorrectPrediction <- sum(myPredictions==myTox)
    print("Number of correct predictions:")
    print(myCorrectPrediction)
    
    mySensitivity <- sum(myPredictions == classLabelPos & myTox == classLabelPos) / sum(myPredictions == classLabelPos)
    print("Sensitivity:")
    print(mySensitivity)
    
    mySpecificity <- sum(myPredictions == classLabelNeg & myTox == classLabelNeg) / sum(myPredictions == classLabelNeg)
    print("Specificity:")
    print(mySpecificity)
    
    myBalancedAccuracy <- (mySensitivity + mySpecificity)/2 
    print("Balanced accuracy:")
    print(myBalancedAccuracy)
    cat("\n")
    
  }
  
  # perform rf
  
  if(rf==T){
    
    myDescriptors_numeric <- as.data.frame(myDescriptors_numeric)
    rownames(myDescriptors_numeric) <- myData[,1]
    
    predictRF <- function(i){
      myRF <- randomForest(myDescriptors_numeric[-c(i),],myTox[-c(i)], myDescriptors_numeric[c(i),],myTox[c(i)], ntree=ntree, importance=T)
      myImportance <- round(importance(myRF)[,4], 2)
      myPredictedClass <- myRF$test$predicted
      myList=list(prediction=myPredictedClass,importance=myImportance)
      return(myList)
    }
    
    myPredictions <- sapply(c(1:dim(myDescriptors_numeric)[1]),predictRF)
    
    myComparisonOfClasses <- data.frame(PredictedClass=unlist(myPredictions[1,]),TrueClass=myTox)
    rownames(myComparisonOfClasses) <- rownames(myDescriptors_numeric)
    print("RF - full model:")
    
    myImportanceMatrix <- do.call(cbind, myPredictions[2,])
    myParameterImportance <- sort(rowMeans(myImportanceMatrix))
    
    print(myParameterImportance)
    print(myComparisonOfClasses)
    
    myCorrectPrediction <- sum(unlist(myPredictions[1,])==myTox)
    print("Number of correct predictions:")
    print(myCorrectPrediction)
    
    mySensitivity <- sum(unlist(myPredictions[1,]) == classLabelPos & myTox == classLabelPos) / sum(unlist(myPredictions[1,]) == classLabelPos)
    print("Sensitivity:")
    print(mySensitivity)
    
    mySpecificity <- sum(unlist(myPredictions[1,]) == classLabelNeg & myTox == classLabelNeg) / sum(unlist(myPredictions[1,]) == classLabelNeg)
    print("Specificity:")
    print(mySpecificity)
    
    myBalancedAccuracy <- (mySensitivity + mySpecificity)/2 
    print("Balanced accuracy:")
    print(myBalancedAccuracy)
    cat("\n")
    
    if(rfe == T){
    
      minValue <- names(myParameterImportance[1])
      a <- which(colnames(myDescriptors_numeric)==minValue)
      
      print("RF - reduced models")
    
      while(length(a) < (dim(myDescriptors)[2]-1)){
      
        predictRF <- function(i){
          myRF <- randomForest(myDescriptors_numeric[,-c(a)][-c(i),],myTox[-c(i)], myDescriptors_numeric[,-c(a)][c(i),],myTox[c(i)], ntree=ntree, importance=T)
          if (importance=="Gini"){
            myImportance <- round(importance(myRF)[,4], 2)
          } else if (importance=="MDA"){
            myImportance <- round(importance(myRF)[,3], 2)
          }
          myPredictedClass <- myRF$test$predicted
          myList=list(prediction=myPredictedClass,importance=myImportance)
          return(myList)
        }
        
        myPredictions <- sapply(c(1:dim(myDescriptors_numeric)[1]),predictRF)
        
        myComparisonOfClasses <- data.frame(PredictedClass=unlist(myPredictions[1,]),TrueClass=myTox)
        rownames(myComparisonOfClasses) <- rownames(myDescriptors_numeric)
        
        myImportanceMatrix <- do.call(cbind, myPredictions[2,])
        myParameterImportance <- sort(rowMeans(myImportanceMatrix))
        
        print(myParameterImportance)
        print(myComparisonOfClasses)
        
        myCorrectPrediction <- sum(unlist(myPredictions[1,])==myTox)
        print("Number of correct predictions:")
        print(myCorrectPrediction)
        
        mySensitivity <- sum(unlist(myPredictions[1,]) == classLabelPos & myTox == classLabelPos) / sum(unlist(myPredictions[1,]) == classLabelPos)
        print("Sensitivity:")
        print(mySensitivity)
        
        mySpecificity <- sum(unlist(myPredictions[1,]) == classLabelNeg & myTox == classLabelNeg) / sum(unlist(myPredictions[1,]) == classLabelNeg)
        print("Specificity:")
        print(mySpecificity)
        
        myBalancedAccuracy <- (mySensitivity + mySpecificity)/2 
        print("Balanced accuracy:")
        print(myBalancedAccuracy)
        cat("\n")
      
        minValue <- names(sort(rowMeans(myImportanceMatrix))[1])
        a_tmp <- which(colnames(myDescriptors_numeric)==minValue)
        a <- c(a,a_tmp)
    
        
        }
    }    
  }
}



