dataframe19 <- read.table("hist.BFT_mer_out")
plot(dataframe19[1:200,], type="l")
plot(dataframe19[20:200,], type="l")
points(dataframe19[20:200,])

sum(as.numeric(dataframe19[20:10001,1]*dataframe19[20:10001,2]))
# 88,974,720,109

max(dataframe19[115:140,2])
# value 115

gsize <- sum(as.numeric(dataframe19[20:10001,1]*dataframe19[20:10001,2]))/115
# 773,693,218

singlecopy<- sum(as.numeric(dataframe19[20:200,1]*dataframe19[20:200,2]))/115
# 641,239,001

(sum(as.numeric(dataframe19[20:200,1]*dataframe19[20:200,2]))) / (sum(as.numeric(dataframe19[20:10001,1]*dataframe19[20:10001,2])))
# 82.9%
