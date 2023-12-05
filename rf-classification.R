#!/usr/bin/env Rscript

# This script creates a classification model using 3 training regions, without band 10, adding NDVI and NDWI, and uses two
# passes of neighborhood information.

# Author: Bruno de Oliveira Schneider - UFLA, 2023

cat('\n=====================================================================\n')
cat('3 TRs, B10 drop, +NDVI +NDWI, 2 neighborhood passes.\n')

library(terra)

# Load training regions delimiters
trainingRegion1 <- vect('regiao-treinamento.shp')
trainingRegion2 <- vect('regiao-treinamento2.shp')
trainingRegion3 <- vect('regiao-treinamento3.shp')
# Load classification region delimiter
classifRegion <- vect('regiao-classificacao.shp')

# Load Sentinel-2 data (higher/reference resolution first)
b2  <- rast('T23KMS_20220625T131301_B02.jp2')
b2T <- crop(b2, trainingRegion1)
# Other data with the same higher resolution
b3  <- rast('T23KMS_20220625T131301_B03.jp2')
b4  <- rast('T23KMS_20220625T131301_B04.jp2')
b8  <- rast('T23KMS_20220625T131301_B08.jp2')

# Load bands that need to be resampled
b1  <- rast('T23KMS_20220625T131301_B01.jp2')
b5  <- rast('T23KMS_20220625T131301_B05.jp2')
b6  <- rast('T23KMS_20220625T131301_B06.jp2')
b7  <- rast('T23KMS_20220625T131301_B07.jp2')
b8a <- rast('T23KMS_20220625T131301_B8A.jp2')
b9  <- rast('T23KMS_20220625T131301_B09.jp2')
b10 <- rast('T23KMS_20220625T131301_B10.jp2')
b11 <- rast('T23KMS_20220625T131301_B11.jp2')
b12 <- rast('T23KMS_20220625T131301_B12.jp2')

# Create multi spectral image
bNIR <- crop(b8, trainingRegion1)
bRed <- crop(b4, trainingRegion1)
bMIR <- resample(b12, b2T)
bNDVI <- (bNIR-bRed) / (bNIR+bRed)
names(bNDVI) <- c('NDVI')
bNDWI <- (bNIR-bMIR) / (bNIR+bMIR)
names(bNDWI) <- c('NDWI')
multiTraining1 <- c(resample(b1, b2T),
                b2T,
                crop(b3, trainingRegion1),
                bRed,
                resample(b5, b2T),
                resample(b6, b2T),
                resample(b7, b2T),
                bNIR,
                resample(b8a, b2T),
                resample(b9, b2T),
                #resample(b10, b2T),
                resample(b11, b2T),
                bMIR)
names(multiTraining1) <- c('coastal', 'blue', 'green', 'red', 'RedE5', 'RedE6', 'RedE7', 'NIR', 'NNIR', 'wVapor',
                           'SWIR_11', 'SWIR_12')

# Create vegetation indices NDVI and NDWI
b2t2 <- crop(b2, trainingRegion2)
bNIR2 <- crop(b8, trainingRegion2)
bRed2 <- crop(b4, trainingRegion2)
bMIR2 <- resample(b12, b2t2)
bNDVI2 <- (bNIR2-bRed2) / (bNIR2+bRed2)
names(bNDVI2) <- c('NDVI')
bNDWI2 <- (bNIR2-bMIR2) / (bNIR2+bMIR2)
names(bNDWI2) <- c('NDWI')
multiTraining2 <- c(resample(b1, b2t2),
                b2t2,
                crop(b3, trainingRegion2),
                bRed2,
                resample(b5, b2t2),
                resample(b6, b2t2),
                resample(b7, b2t2),
                bNIR2,
                resample(b8a, b2t2),
                resample(b9, b2t2),
                #resample(b10, b2T),
                resample(b11, b2t2),
                bMIR2)
names(multiTraining2) <- c('coastal', 'blue', 'green', 'red', 'RedE5', 'RedE6', 'RedE7', 'NIR', 'NNIR', 'wVapor',
                           'SWIR_11', 'SWIR_12')
b2t3 <- crop(b2, trainingRegion3)
bNIR3 <- crop(b8, trainingRegion3)
bRed3 <- crop(b4, trainingRegion3)
bMIR3 <- resample(b12, b2t3)
bNDVI3 <- (bNIR3-bRed3) / (bNIR3+bRed3)
names(bNDVI3) <- c('NDVI')
bNDWI3 <- (bNIR3-bMIR3) / (bNIR3+bMIR3)
names(bNDWI3) <- c('NDWI')
multiTraining3 <- c(resample(b1, b2t3),
                b2t3,
                crop(b3, trainingRegion3),
                bRed3,
                resample(b5, b2t3),
                resample(b6, b2t3),
                resample(b7, b2t3),
                bNIR3,
                resample(b8a, b2t3),
                resample(b9, b2t3),
                #resample(b10, b2t3),
                resample(b11, b2t3),
                bMIR3)
names(multiTraining3) <- c('coastal', 'blue', 'green', 'red', 'RedE5', 'RedE6', 'RedE7', 'NIR', 'NNIR', 'wVapor',
                           'SWIR_11', 'SWIR_12')


# Classification region data
b2C <- crop(b2, classifRegion)
bNIRC <- crop(b8, classifRegion)
bRedC <- crop(b4, classifRegion)
bMIRC <- resample(b12, b2C)
bNDVIC <- (bNIRC-bRedC) / (bNIRC+bRedC)
names(bNDVIC) <- c('NDVI')
bNDWIC <- (bNIRC-bMIRC) / (bNIRC+bMIRC)
names(bNDWIC) <- c('NDWI')
multiClassif <- c(resample(b1, b2C),
                b2C,
                crop(b3, classifRegion),
                bRedC,
                resample(b5, b2C),
                resample(b6, b2C),
                resample(b7, b2C),
                bNIRC,
                resample(b8a, b2C),
                resample(b9, b2C),
                #resample(b10, b2C),
                resample(b11, b2C),
                bMIRC)
names(multiClassif) <- c('coastal', 'blue', 'green', 'red', 'RedE5', 'RedE6', 'RedE7', 'NIR', 'NNIR', 'wVapor',
                         'SWIR_11', 'SWIR_12')

# Plot RGB image from classification image if intended
#plotRGB(multiTraining1, 4,3,2)
#img <- multiClassif[[c(4,3,2)]] # copy from the original set

# Load verified EMPAMIG/EMATER coffee polygons
cropPolys <- vect('cafe-lavras-edt.shp')
trainingPolys1 <- crop(cropPolys, trainingRegion1)
trainingPolys1$AREA_ha <- NULL # remove area column
trainingPolys2 <- crop(cropPolys, trainingRegion2)
trainingPolys2$AREA_ha <- NULL # remove area column
trainingPolys3 <- crop(cropPolys, trainingRegion3)
trainingPolys3$AREA_ha <- NULL # remove area column
classifPolys <- crop(cropPolys, classifRegion)
classifPolys$AREA_ha <- NULL # remove area column

# Rasterize crop polygons / balancing coffee samples
rasterClasses <- rasterize(trainingPolys1, b2T, field='Categoria', background=0)
names(rasterClasses) <- c('class')
trainingData1 <- as.data.frame(c(rasterClasses, multiTraining1, bNDVI, bNDWI))
rasterClasses2 <- rasterize(trainingPolys2, b2t2, field='Categoria', background=0)
names(rasterClasses2) <- c('class')
trainingData2 <- as.data.frame(c(rasterClasses2, multiTraining2, bNDVI2, bNDWI2))
trainingData2 <- trainingData2[trainingData2$class > 0, ] # filter non-coffee samples
rasterClasses3 <- rasterize(trainingPolys3, b2t3, field='Categoria', background=0)
names(rasterClasses3) <- c('class')
trainingData3 <- as.data.frame(c(rasterClasses3, multiTraining3, bNDVI3, bNDWI3))
trainingData3 <- trainingData3[trainingData3$class > 0, ] # filter non-coffee samples
fullTrainingData <- rbind(trainingData1, rbind(trainingData2, trainingData3))
cat(paste(nrow(fullTrainingData), 'samples used in training.\n'))
cat('Number of samples for each class:')
print(table(fullTrainingData$class))

# Create dataframe for classification
classifData <- as.data.frame(c(multiClassif, bNDVIC, bNDWIC))

# Define inputs and output
predictors <- c(names(multiTraining1), 'NDVI', 'NDWI')
form <- as.formula(paste('class ~', paste(predictors, collapse='+'))) #warning: names must not contain spaces

# Create classification model (this may take some time...)
library("ranger")
classifModel <- ranger(form, data=fullTrainingData, importance='impurity', classification=TRUE)
cat('Band importance ranking:\n')
print(classifModel$variable.importance[order(classifModel$variable.importance)])
cat(paste('Classification error (OOB):', classifModel$prediction.error, '\n'))

# Classify!
pred <- predict(classifModel, data=classifData, type='response')

# Accuracy function
ovAcc <- function(conmat) {
    # number of total cases/samples
    n = sum(conmat)
    # number of correctly classified cases per class
    diag = diag(conmat)
    # Overall Accuracy
    OA = sum(diag) / n
    # observed (true) cases per class
    rowsums = apply(conmat, 1, sum)
    p = rowsums / n
    # predicted cases per class
    colsums = apply(conmat, 2, sum)
    q = colsums / n
    expAccuracy = sum(p*q)
    kappa = (OA - expAccuracy) / (1 - expAccuracy)
    # Producer accuracy
    PA <- diag / colsums
    # User accuracy
    UA <- diag / rowsums
    outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
    #print(outAcc)

    global_acc = data.frame(overallAccuracy=OA, overallKappa=kappa)
    #print(global_acc)
    cat('Producer accuracy for each class:\n')
    print(PA)
    cat(paste('Overall accuracy:', format(OA, digits=4), 'Kappa:', format(kappa, digits=4), '\n'))
    # Based from: http://gsp.humboldt.edu/olm/Courses/GSP_216/lessons/accuracy/metrics.html
}


# Create image from classification
rastClassif1 <- rast(ncols=ncol(b2C), nrows=nrow(b2C), nlyrs=1, crs=crs(b2C), extent=ext(b2C))
values(rastClassif1) <- pred$predictions
names(rastClassif1) <- c('class')

# Save to file, if intended. Supposing 2 classes, generates values up to 200.
#writeRaster(rastClassif1*100, 'classificacao.png', overwrite=TRUE)
#plot(rastClassif1)
#plot(classifPolys, add=TRUE)

# Rasterize coffee polygons to generation a true classification
trueClassif <- terra::rasterize(classifPolys, b2C, field='Categoria', background=0)

# Confusion Matrix
realFactors <- factor(values(trueClassif), levels=c(0:2))
compFactors <- factor(pred$predictions, levels=c(0:2))
conMat <- table(compFactors, realFactors)
ovAcc(conMat)

# Generate raster from first neighborhood information. TR1.
cat('Computing neighborhood data\n')
neighborhood1 <- terra::rasterize(trainingPolys1[trainingPolys1$Categoria==1], b2T, field=1, background=0)
neighborhood1 <- terra::focal(neighborhood1, w=5, fun='sum', fillvalue=0)
names(neighborhood1) <- c('nbh1')
neighborhood2 <- terra::rasterize(trainingPolys1[trainingPolys1$Categoria==2], b2T, field=1, background=0)
neighborhood2 <- terra::focal(neighborhood2, w=5, fun='sum', fillvalue=0)
names(neighborhood2) <- c('nbh2')
multiTraining1 <- c(multiTraining1, neighborhood1, neighborhood2)
#plot(neighborhood1)
#print("Neighborhood histogram:")
#print(table(as.vector(neighborhood1)))

# Generate raster from first neighborhood information. TR2.
neighborhood1 <- terra::rasterize(trainingPolys2[trainingPolys2$Categoria==1], b2t2, field=1, background=0)
neighborhood1 <- terra::focal(neighborhood1, w=5, fun='sum', fillvalue=0)
names(neighborhood1) <- c('nbh1')
neighborhood2 <- terra::rasterize(trainingPolys2[trainingPolys2$Categoria==2], b2t2, field=1, background=0)
neighborhood2 <- terra::focal(neighborhood2, w=5, fun='sum', fillvalue=0)
names(neighborhood2) <- c('nbh2')
multiTraining2 <- c(multiTraining2, neighborhood1, neighborhood2)
#plot(neighborhood1)

# Generate raster from first neighborhood information. TR3.
neighborhood1 <- terra::rasterize(trainingPolys3[trainingPolys3$Categoria==1], b2t3, field=1, background=0)
neighborhood1 <- terra::focal(neighborhood1, w=5, fun='sum', fillvalue=0)
names(neighborhood1) <- c('nbh1')
neighborhood2 <- terra::rasterize(trainingPolys3[trainingPolys3$Categoria==2], b2t3, field=1, background=0)
neighborhood2 <- terra::focal(neighborhood2, w=5, fun='sum', fillvalue=0)
names(neighborhood2) <- c('nbh2')
multiTraining3 <- c(multiTraining3, neighborhood1, neighborhood2)
#plot(neighborhood1)

# Create new classification model using neighborhood data (this may take some time...)
trainingData1 <- as.data.frame(c(rasterClasses, multiTraining1, bNDVI, bNDWI))
trainingData2 <- as.data.frame(c(rasterClasses2, multiTraining2, bNDVI2, bNDWI2))
trainingData2 <- trainingData2[trainingData2$class>0, ] # filtrar a classe 0
trainingData3 <- as.data.frame(c(rasterClasses3, multiTraining3, bNDVI3, bNDWI3))
trainingData3 <- trainingData3[trainingData3$class>0, ] # filtrar a classe 0
# Join all data
fullTrainingData <- rbind(trainingData1, rbind(trainingData2, trainingData3))
predictors <- c(names(multiTraining1), 'NDVI', 'NDWI')
form <- as.formula(paste('class ~', paste(predictors, collapse='+')))
cat('Creating new classification model\n')
classifModel <- ranger(form, data=fullTrainingData, importance="impurity", classification=TRUE)
cat(paste('Classification error (OOB):', classifModel$prediction.error, '\n'))

# Create neighborhood data for CR
neighborhood1 <- rastClassif1
neighborhood1[neighborhood1$class==2] = 0 # filter class 2
neighborhood1 <- terra::focal(neighborhood1, w=5, fun='sum', fillvalue=0)
names(neighborhood1) <- c('nbh1')
neighborhood2 <- rastClassif1
neighborhood2[neighborhood2$class==1] = 0 # filter class 1
neighborhood2[neighborhood2$class==2] = 1
neighborhood2 <- terra::focal(neighborhood2, w=5, fun='sum', fillvalue=0)
names(neighborhood2) <- c('nbh2')

classifData <- as.data.frame(c(multiClassif, neighborhood1, neighborhood2, bNDVIC, bNDWIC))
cat('Second mapping:\n')
pred <- predict(classifModel, data=classifData, type='response')
rastClassif2 <- rastClassif1
values(rastClassif2) <- pred$prediction
#writeRaster(rastClassif2*51, 'classificacao.png', overwrite=TRUE)

# Update confusion matrix
compFactors <- factor(pred$predictions, levels=c(0:2))
conMat <- table(compFactors, realFactors)
ovAcc(conMat)
#writeRaster(rastClassif2*100, 'classificacao-it1.png', overwrite=TRUE)

cat('\nSecond neighborhood iteration\n')
neighborhood1 <- rastClassif2
neighborhood1[neighborhood1$class==2] = 0 # filter class 2
neighborhood1 <- terra::focal(neighborhood1, w=5, fun='sum', fillvalue=0)
names(neighborhood1) <- c('nbh1')
neighborhood2 <- rastClassif2
neighborhood2[neighborhood2$class==1] = 0 # filter class 1
neighborhood2[neighborhood2$class==2] = 1
neighborhood2 <- terra::focal(neighborhood2, w=5, fun='sum', fillvalue=0)
names(neighborhood2) <- c('nbh2')

classifData <- as.data.frame(c(multiClassif, neighborhood1, neighborhood2, bNDVIC, bNDWIC))
cat('Third mapping:\n')
pred <- predict(classifModel, data=classifData, type='response')
rastClassif3 <- rastClassif1
values(rastClassif3) <- pred$prediction
compFactors <- factor(pred$predictions, levels=c(0:2))
conMat <- table(compFactors, realFactors)
ovAcc(conMat)
#writeRaster(rastClassif3*100, 'classificacao-it2.png', overwrite=TRUE)

