# ============================================================================
cat('\n *** loading dependences and native functions ***\n\n')
require(sp)
require(rgdal)
require(raster)
require(gstat)
require(ggplot2)
require(gridExtra)
require(ggdendro)
require(rrcov)
require(GGally)
require(cluster)
# ============================================================================
# gstat variogram function
.gstat.semvar <- function(spdf, treshold=3, np='NA') {
    maxDist <- sqrt((max(spdf@coords[,1]) - min(spdf@coords[,1]))^2 + (max(spdf@coords[,2]) - min(spdf@coords[,2]))^2)
    if (np == 'NA') {
        semvar <- variogram(spdf@data[, 1] ~ 1, locations = spdf, data = spdf, cutoff = maxDist/treshold)
    } else {
        semvar <- variogram(spdf@data[, 1] ~ 1, locations = spdf, data = spdf, cutoff = maxDist/treshold, width = maxDist/np)
    }
    return(semvar)
}
# ============================================================================
# selection model
.selection.model <- function(semvar, spdf) {
    stats <- data.frame()
    for (i in c('Cir', 'Sph', 'Exp', 'Gau')) {
        fit <- .fit(semvar, i, spdf)
        fit.stats <- .fit.stats(semvar, spdf, i, fit$nugg, fit$sill, fit$rang)
        stats <- rbind(stats, cbind(fit, fit.stats))
    }
    return(stats)
}

# gstat Ordinary kriging
.gstat.Ok <- function (spdf, grd, model, nugg=0, sill=0, rang=0) {
    gstat.model <- vgm(psill = sill, model = model, nugget = nugg, range = rang)
    return(krige(formula = spdf@data[, 1] ~ 1, locations = spdf, newdata = grd, model = gstat.model, debug.level=0))
}

# gstat lave one out 
.gstat.loo <- function(spdf, model, nugg, sill, rang) {
    gstat.model <- vgm(psill = sill, model = model, nugget = nugg, range = rang)
    estim <- data.frame()
    for (i in 1:nrow(spdf@data)) {
        ok <- krige(formula = spdf[-i,]@data[, 1] ~ 1, locations = spdf[-i,], newdata = spdf[i,], model = gstat.model, debug.level=0)
        estim <- rbind(estim, ok@data)
    }
    return(estim)
}

# plot models 
.plot.model <- function(Object, model, nugg, sill, rang, Color='blue') {
    self.model <- .vgm(model, nugg, sill, rang)
    # Object: spdf / ggplot
    if (class(Object)[1] == 'gg') {
        self.plot <- Object
        self.plot <- self.plot + stat_function(fun=self.model, colour = Color)
    } else {
        semvar <- Object
        self.plot <- ggplot(data=semvar) + geom_point(aes(x=dist, y=gamma, alpha=np))
        self.plot <- self.plot + stat_function(fun=self.model, colour = Color)
    }
    return(self.plot)
}

# ============================================================================
#                           PCA

pca <- prcomp(dataset, scale=F, center=T)
pcaplot <- ggbiplot(pca, obs.scale = 1, var.scale = 1, circle = TRUE)
pca$rotation
# ============================================================================
#                           DENDOGRAM

# multiple dendogram ---------------------------------------------------------
.exploreMethod <- function(dataset, metDist = 'euclidean') 
{
    Plots <- list()
    for(i in c('ward', 'single', 'complete', 'average', 'median', 'centroid')) 
    {
        hc <- hclust(dist(dataset, metDist), method = i)
        Plots <- append(Plots, list(ggplot(segment(dendro_data(as.dendrogram(hc)))) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + labs(title = i, x = '', y = '')))
    }
    return(Plots)
}
.cutDendro <- function(dataset, Method = 'ward', metDist = 'euclidean') 
{
    hc <- hclust(dist(dataset, method = metDist), method = Method)
    h <- min(hc$height[hc$height > mean(hc$height) + 1.25*sd(hc$height)]) # Mojena (1977)
    k <- 1 + length(hc$height[hc$height > mean(hc$height) + 1.25*sd(hc$height)]) 
    return(list(k = k, h = h, group = cutree(hc, k = k)))
}

grp <- .cutDendro(dataset)$group

# ============================================================================

# CENTROIDS ------------------------------------------------------------------  
.centroids <- function (dataset, group) {
    df <- data.frame(group = group, dataset)
    centroids <- data.frame()
    for (i in 1:length(levels(factor(group)))) {
        centroids <- rbind(centroids, colMeans(df[df$group == i,]))
    }
    names(centroids) <- names(df)
    return(centroids)
}

cent <- .centroids(dataset, grp)

# ============================================================================
#                       K-means

K <- kmeans(dataset, centers = as.matrix(cent[, 2:length(cent)]))

# ----------------------------------------------------------------------------

tb <- table(K$cluster, grp)
sum(diag(tb)/sum(tb))

# PLOTS

clusplot(dataset, K$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

scores <- data.frame(pca$x, 
                     label = row.names(dataset), 
                     # categories = jan.phyto$Area, 
                     grDendro = grp, 
                     grKmeans = K$cluster)

p <- ggplot(scores) +
    geom_point(aes(x = PC1, y = PC2), size = 2) +
    stat_ellipse(aes(x = PC1,y = PC2, fill = factor(grDendro)), geom = "polygon", level = 0.95, alpha = 0.2) +
    geom_text(aes( label = label, x = PC1, y = PC2), size = 2, angle = 0, vjust = 3)

# ============================================================================

# Dendogram
hc <- hclust(dist(dataset), method='ward')
hcData <- dendro_data(as.dendrogram(hc))
scores <- merge(scores, hcData$labels, by = 'label')

dendroplot <- ggplot(segment(hcData))
dendroplot <- dendroplot + geom_segment(aes(x = x, y = y, xend = xend, yend = yend))
dendroplot <- dendroplot + geom_text(data = scores, aes(label= label, x = x, y = -1, color = factor(grKmeans)), size = 5, angle = 90)