##=============================================================
## 펫박스 추천시스템 개발
##=============================================================

##------------------------
## 패키지 인스톨
##------------------------
# installing and loading readxl package 
install.packages("dbscan")
install.packages('readxl')
install.packages('proxy')
install.packages("ggdendro")
install.packages("dplyr")
install.packages("ggplot2")
install.packages('stylo')
install.packages("philentropy")
install.packages("microbenchmark")

library(recommenderlab)
library(microbenchmark)
library(philentropy)
library(stylo)
library(proxy)


##----------------------------------------------------
# 표준화 및 정규화 함ㅅ
##----------------------------------------------------
scale<-function(x) {
  return(x-min(x))/(max(x)-min(x))
}
normalize <- function(x) {
  return((x-mean(x))/sd(x))
}


##----------------------------------------------------
# cosine distance function
##----------------------------------------------------
cosine_Dist <- function(x,y){
  jaccard_dst <- distance(as.matrix(y), method = "jaccard")
  as.dist(((1-(x%*%t(x))/(sqrt(rowSums(x^2)%*%t(rowSums(x^2)))))+jaccard_dst)/2)
}


##----------------------------------------------------
## 클러스터링 및 결과 출력 함수(Hierachical)
##----------------------------------------------------
display_cluster_h <- function(x, method){
  hc <- hclust(x, method)
  cut <- as.data.frame(cutree(hc, k=4))
  cut$names <- rownames(cut)
  names(cut) <- c("cut", "names")
  
  library(ggplot2)
  library(ggdendro)
  library(dplyr)
  
  hcdata <- dendro_data(hc, type="rectangle")
  hcdata$labels <- left_join(hcdata$labels, cut, by=c("label"="names"))
  
  ggplot(hcdata$segments) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    geom_text(data = hcdata$labels, aes(x, y, label = label, colour=factor(cut)), 
              hjust = 1, size = 4) +
    scale_color_manual(values=c("red", "blue", "green", "black"), guide_legend(title="clusters")) +
    labs(x="", y="") + coord_flip() + theme_bw()
}

##-------------------------------
## 입력데이터 분석 및 정제
##-------------------------------
analyzeData <- function (x, answer) {
  
  # 데이터 분포도 출력
  boxplot(x)
  head(as.matrix(x))
  summary(x)
  
  # 산점도 출력
  plot(x) 
  text(2.5,4,"abc",cex=6)
  pairs(x, panel=panel.smooth, cex.labels=2)
  
  # 산점도 Class 적용하여 출력
  library(mclust)
  clPairs(x, answer)
  
  # 단백질이 0.2이하인 것 제거
  #subset(st_all_main, Protein < 0.2) 
  
  #PCA(주성분 분석) 차원축소
  
  #
}


##-------------------------------
## 훈련 및 테스트 데이터로 나누기
##-------------------------------
divideTrainningSet <- function(st_all_main) {
  library(caret)
  set.seed(123)
  inTrain <- createDataPartition(st_all_main$Protein, p=0.7, list=FALSE)
  trainings <- st_all_main[inTrain,]
  testings <- st_all_main[-inTrain,]
  return(list(training=trainings, testing=testings)) 
}


##------------------------
## 거리 계산
##------------------------
calculateDistance <- function(x, mthd) {
  return(distance(as.matrix(x), method = mthd))
}


##------------------------
# 하이라키컬 클러스터링
##------------------------
hClustering <- function(distance, method, value, m) {
  hc <- hclust(as.dist(distance),method)
  cluster <- cutree(hc, h=value)
  plot(hc, hang=-1, cex=.7, main=m)
  
  #library(ape)
  #plot(as.phylo(hc),cex = .7)
  
  abline(h=value, col="blue")
  rect.hclust(hc, h=value, border = "red")
  return(cluster)
}



##----------------------------
# K-Means 클러스터링
##----------------------------
kmeansClustering <- function(x, numOfCluster) {
  wssPlot(x)
  result.kmeans <- kmeans(x, centers = numOfCluster, iter.max = 10000)
  result.kmeans$centers
  cluster<-as.factor(result.kmeans$cluster)
  return(cluster)
}



##-------------------------------------
## K-Means 클러스터링에서 K값 구하기
##-------------------------------------
wssPlot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)
  }
  plot(1:nc, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
}


##----------------------------
## 아이템 추천하기
##----------------------------
recommend_Contents <- function(x, cluter, item_number) {
  main_m = as.matrix(x)
  cat("요청 아이템 : ", main_m[item_number,1], sep="\n\n")
  cat("추천 아이템 : ", main_m[cluster==cluster[item_number],1], sep = "\n\n")
}


##----------------------------
## Knn Classification
##----------------------------
knnClassification <- function(trainData, train_label, testData, test_label, K ) {
  pre <- knn(train=trainData,
             cl = train_label,
             k = K)
  #결과보기
  pre
  table(pred=pre, true=test_label)
}


GetMax <-function(x) {
  for (m in 1:nrow(x)) {
    y <- replace(x[m,], m, values=9999)
    y <- format(y,scientific=FALSE)
    cat(m, which.min(y), min(y), "\n")
  }
}

GetCosMax <-function(x) {
  for (m in 1:nrow(x)) {
    y <- replace(x[m,], which(x[m,]==1), values=0)
    y <- format(y,scientific=FALSE)
    cat(which.max(y), max(y), "\n")
  }
}



######################################################
##------------------------
## Main
##------------------------
######################################################
library(readxl)
data_main <- read_excel("d:/test.xlsx", sheet = "all", range = "C1:J100", col_names = TRUE, col_types = "guess",na = "NA")
data_sub <- read_excel("d:/test.xlsx", sheet = "all", range = "W1:Aq100", col_names = TRUE, col_types = "guess",na = "NA")
data_class <- read_excel("d:/test.xlsx", sheet = "all", range = "K1:K100", col_names = TRUE, col_types = "guess",na = "NA")
data_main <- cbind(data_main, data_sub)

analyzeData(data_main[,2:8], data_class)          #데이터 분석
nrow(data_main)
data_main <- unique(data_main)
nrow(data_main)

tra_datas <- divideTrainningSet(data_main)               #Trainning Set 분리

tra_datas$training.data <- scale(tra_datas$training[,2:8]) #정규화
tra_datas$training.data <- cbind(tra_datas$training.data,tra_datas$training[,9:29])

library(philentropy)
getDistMethods()

myDistance <- as.matrix(cosine_Dist(as.matrix(tra_datas$training.data[,1:7]), tra_datas$training.data[,8:28]))
#myDistance <- calculateDistance(tra_datas$training.data, method) #거리계산
myDistance <- distance(tra_datas$training.data[,1:7], method = "euclidean")
GetMax(myDistance)

install.package('wordspace')
library(wordspace)
dist_mat=dist.matrix(tra_datas$training.data[,1:7],method='euclidean')
dist_mht <- dist(tra_datas$training.data[,1:7], method = "euclidean")



#rownames(myDistance) = t(tra_datas$training[,1])
#colnames(myDistance) = t(tra_datas$training[,1])

#Classification
#knnClassification()

#Clustering
hCluster <- hClustering(myDistance, "complete", 0.1, "fdf") #single linkage, complet linkage, centroid, Ward

kmCluster <- kmeansClustering(tra_datas$training.data, 2)

#Clustering 결과 보기
qplot(Protein,Ash,colour=kmCluster,data=tra_datas$training)
table(tra_datas$training$Class, kmCluster)

library(cluster)
clusplot(tra_datas$training.data[,1:7], kmCluster, color=TRUE, shade=FALSE, labels=13, lines=FALSE)

#추천하기
recommend_Contents(tra_datas$training, hCluster, 1)

#파일에 저장하기
write.csv(as.matrix(dist_euc),file="d:/dog euclidian distance results.csv", row.names=FALSE)
write.csv(as.matrix(myDistance),file="d:/myDistance.csv", row.names=FALSE)
write.csv(as.matrix(tra_datas$training),file="d:/Trainning Set.csv", row.names=FALSE)














##----------------------------
## 그룹핑
##----------------------------
grouping <- cutree(hc, h = 0.1)
grouping
sort(table(grouping), decreasing = TRUE)







#Use k-nearest Neighbors
library("dbscan")기
nn <- kNN(distance,5)
main_m = as.matrix(data_main)
cat("요청 아이템 : ", main_m[15,1], sep="\n\n")
cat("추천 아이템 : ", main_m[nn$id[15,],1], sep = "\n\n")





##-----------------------------
## 막대다이어그램으로 출력하기
##-----------------------------

#plot(st_all, xlab = "linear algebra", ylab = "machine learning", xlim = c(50, 100), ylim = c(50, 100), main = "scatter plot of final scores")
#text(st_all[,1], st_all[,2], labels = abbreviate(rownames(st_all)), cex = 0.8, pos = 1, col = "blue")
#points(st_all[1,1], st_all[1,2], col = "red", pch = 19) 

data <- dist_matrix_euc[1,]
xx <- barplot(data, col = "red", border = NA, xlab = "number of data", ylab = "value of distance")
#text(xx, data, labels = round(data,0), pos = 3, cex = 0.8, col = "black")

#image(as.matrix(dist_max), main = "Maximum Similarity") // ?Ÿ? ��?絵 ?? ????



