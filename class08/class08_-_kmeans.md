class08
================
Reshma Pyala
2/5/2019

K-means clustering
------------------

``` r
tmp <- c(rnorm(30,-3), rnorm(30,3))
x <- cbind(x=tmp, y=rev(tmp))
plot(x)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-1-1.png)

Use the kmeans() function setting k to 2 and nstart=20

Inspect/print the results

Q. How many points are in each cluster? Q. What ‘component’ of your result object details - cluster size? - cluster assignment/membership? - cluster center?

Plot x colored by the kmeans cluster assignment and add cluster centers as blue points

``` r
km <- kmeans(x, centers = 2, nstart = 20)
km
```

    ## K-means clustering with 2 clusters of sizes 30, 30
    ## 
    ## Cluster means:
    ##           x         y
    ## 1 -2.637778  2.795078
    ## 2  2.795078 -2.637778
    ## 
    ## Clustering vector:
    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ## 
    ## Within cluster sum of squares by cluster:
    ## [1] 53.61634 53.61634
    ##  (between_SS / total_SS =  89.2 %)
    ## 
    ## Available components:
    ## 
    ## [1] "cluster"      "centers"      "totss"        "withinss"    
    ## [5] "tot.withinss" "betweenss"    "size"         "iter"        
    ## [9] "ifault"

``` r
km$size
```

    ## [1] 30 30

Clusster membership assignment vector (i.e. which cluster group my data lies in!)

``` r
km$cluster
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
plot(x, col=km$cluster)
points(km$centers, col="blue", pch=16, cex=1.5)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
km$totss
```

    ## [1] 992.7104

Hierarchial clustering in R
---------------------------

``` r
# First we need to calculate point (dis)similarity
#   as the Euclidean distance between observations
dist_matrix <- dist(x)
# The hclust() function returns a hierarchical
#  clustering model
hc <- hclust(d = dist_matrix)
# the print method is not so useful here
hc
```

    ## 
    ## Call:
    ## hclust(d = dist_matrix)
    ## 
    ## Cluster method   : complete 
    ## Distance         : euclidean 
    ## Number of objects: 60

``` r
View(as.matrix(dist_matrix))
dim(as.matrix(dist_matrix))
```

    ## [1] 60 60

``` r
plot(hc)
abline(h=6, col="red")
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
cutree(hc, h=6) #Cut by height h
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
#Can also do this function as a way to divide it into groups without picking a specific height yourself
cutree(hc, k=2)
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
# Using different hierarchical clustering methods
d <- dist_matrix
hc.complete <- hclust(d, method="complete")
hc.average  <- hclust(d, method="average")
hc.single   <- hclust(d, method="single")

plot(hc.single)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
# Step 1. Generate some example data for clustering
x <- rbind(
  matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2),   # c1
  matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
  matrix(c(rnorm(50, mean = 1, sd = 0.3),           # c3
           rnorm(50, mean = 0, sd = 0.3)), ncol = 2))
colnames(x) <- c("x", "y")
# Step 2. Plot the data without clustering
plot(x)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
# Step 3. Generate colors for known clusters
#         (just so we can compare to hclust results)
col <- as.factor( rep(c("c1","c2","c3"), each=50) )
plot(x, col=col, pch=16, cex=.5)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-12-2.png)

Q. Use the dist(), hclust(), plot() and cutree() functions to return 2 and 3 clusters Q. How does this compare to your known 'col' groups?

``` r
d <- dist(x)
hc2 <- hclust(d)
plot(hc2)
abline(h=2, col="red")
abline(h=2.5, col="blue")
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
gp2 <- cutree(hc2, k=2)
gp3 <- cutree(hc2, k=3)
plot(x, col=gp2)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
plot(x, col=gp3)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-14-2.png)

General code for kmeans and hclust **kmeans(x, centers=2, nstart=20)** **hclust(dist(x))**

Could also do it all in one run like this: **plot(hclust(dist(x)))** but once you run it, you don't have the output saved anywhere so it is better to save things separately to an object and then run more lines of code so that you can use the object name later to run it instead of recalculating for it

Now trying to do PCA (Principal Component Analysis) in R
--------------------------------------------------------

``` r
mydata <- read.csv("https://tinyurl.com/expression-CSV",
row.names=1)

head(mydata)
```

    ##        wt1 wt2  wt3  wt4 wt5 ko1 ko2 ko3 ko4 ko5
    ## gene1  439 458  408  429 420  90  88  86  90  93
    ## gene2  219 200  204  210 187 427 423 434 433 426
    ## gene3 1006 989 1030 1017 973 252 237 238 226 210
    ## gene4  783 792  829  856 760 849 856 835 885 894
    ## gene5  181 249  204  244 225 277 305 272 270 279
    ## gene6  460 502  491  491 493 612 594 577 618 638

``` r
#Let's do PCA
pca <- prcomp(t(mydata), scale=TRUE)

#See what is returned by the procomp() function
attributes(pca)
```

    ## $names
    ## [1] "sdev"     "rotation" "center"   "scale"    "x"       
    ## 
    ## $class
    ## [1] "prcomp"

``` r
pca$names
```

    ## NULL

``` r
##A basic PC1 vs PC2 2D plot
plot(pca$x[,1], pca$x[,2])
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
##Variance captured per PC
pca.var <- pca$sdev^2
##Percent variance is often more informative to look at
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.var.per
```

    ##  [1] 92.6  2.3  1.1  1.1  0.8  0.7  0.6  0.4  0.4  0.0

``` r
#Now let's plot it
 pca.var <- pca$sdev^2
 pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
 barplot(pca.var.per, main="Scree Plot",
 xlab="Principal Component", ylab="Percent Variation")
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
## A vector of colors for wt and ko samples
colvec <- colnames(mydata)
colvec[grep("wt", colvec)] <- "red"
colvec[grep("ko", colvec)] <- "blue"

plot(pca$x[,1], pca$x[,2], col=colvec, pch=16,
xlab=paste0("PC1 (", pca.var.per[1], "%)"),
ylab=paste0("PC2 (", pca.var.per[2], "%)"))
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
plot(pca$x[,1], pca$x[,2], col=colvec, pch=16,xlab=paste0("PC1 (", pca.var.per[1], "%)"),ylab=paste0("PC2 (", pca.var.per[2], "%)"))
## Click to identify which sample is which
identify(pca$x[,1], pca$x[,2], labels=colnames(mydata))
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-21-1.png)

    ## integer(0)

Next hands on section
---------------------

``` r
x <- read.csv("UK_foods.csv", row.names =1)
head(x)
```

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

``` r
#*Came back and added the row.names and head so that the file is read with proper formatting and this way we don't have to adjust it later and complicate things even more*

#How many rows and columns are in your new data frame named *x*? What R functions could you use to answer this queestions?
#I opened up the data saved in the environment panel and saw that there are 17 rows and 5 columns
```

``` r
## Complete the following code to find out how many rows and columns are in x?
dim(x)
```

    ## [1] 17  4

``` r
ncol(x)
```

    ## [1] 4

``` r
nrow(x)
```

    ## [1] 17

``` r
#This is reading as 5 columns because it is counting x as a label but that space is meant to be blank - as a cross section between the two variables
```

``` r
## Preview the first 6 rows
head(x)
```

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

``` r
#Use head(x) to preview the first couple data points in this data. and use tail(x) to preview the last couple data points in this data.
tail(x)
```

    ##                   England Wales Scotland N.Ireland
    ## Fresh_fruit          1102  1137      957       674
    ## Cereals              1472  1582     1462      1494
    ## Beverages              57    73       53        47
    ## Soft_drinks          1374  1256     1572      1506
    ## Alcoholic_drinks      375   475      458       135
    ## Confectionery          54    64       62        41

``` r
# Note how the minus indexing works
rownames(x) <- x[,1]
x <- x[,-1]
head(x)
```

    ##     Wales Scotland N.Ireland
    ## 105   103      103        66
    ## 245   227      242       267
    ## 685   803      750       586
    ## 147   160      122        93
    ## 193   235      184       209
    ## 156   175      147       139

``` r
#Using -1 means to show all the data but column 1
```

``` r
dim(x)
```

    ## [1] 17  3

``` r
#This should now have 4 columns since we condensed the data a bit
```

``` r
x <- read.csv("UK_foods.csv", row.names=1)
head(x)
```

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

``` r
barplot(as.matrix(x), beside=T, col=rainbow(nrow(x)))
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
barplot(as.matrix(x), beside=F, col=rainbow(nrow(x)))
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-29-1.png)

``` r
pairs(x, col=rainbow(10), pch=16)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-30-1.png)

Let's run pca on this to get some insight

``` r
# Use the prcomp() PCA function 
pca <- prcomp( t(x) )
summary(pca)
```

    ## Importance of components:
    ##                             PC1      PC2      PC3       PC4
    ## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
    ## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    ## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

``` r
# Plot PC1 vs PC2
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500))
text(pca$x[,1], pca$x[,2], colnames(x),col=c("orange", "red", "blue", "darkgreen"))
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-32-1.png)

``` r
v <- round( pca$sdev^2/sum(pca$sdev^2) * 100 )
v
```

    ## [1] 67 29  4  0

``` r
## or the second row here...
z <- summary(pca)
z$importance
```

    ##                              PC1       PC2      PC3          PC4
    ## Standard deviation     324.15019 212.74780 73.87622 4.188568e-14
    ## Proportion of Variance   0.67444   0.29052  0.03503 0.000000e+00
    ## Cumulative Proportion    0.67444   0.96497  1.00000 1.000000e+00

``` r
barplot(v, xlab="Principal Component", ylab="Percent Variation")
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-35-1.png)

### Examine the "loadings"

This will help us determine how the original variables (dimensions) contribute to our new PCs

``` r
## Lets focus on PC1 as it accounts for > 90% of variance 
par(mar=c(10, 3, 0.35, 0))
barplot( pca$rotation[,1], las=2 )
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-36-1.png)

``` r
## The inbuilt biplot() can be useful for small datasets 
biplot(pca)
```

![](class08_-_kmeans_files/figure-markdown_github/unnamed-chunk-37-1.png)

``` r
##prcomp function general set up

# prcomp(t(x))

##Main functions from today

#kmeans(x, centers = 2, nstart = 20)
#hclust(dist(x))
#prcomp(t(x))
```
