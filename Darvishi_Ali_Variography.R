# Variography
# Geomodelling - Geostatistics for Natural Resource Evaluation
# Name: Ali Darvishi Koushali
# Dataset: testData123.csv

# libraries
install.packages("sp", dependencies=TRUE)
install.packages("gstat", dependencies=TRUE)

# load the libraries
library(sp)
library(gstat)
library(MASS)

# inking
cf = function(x) {
  crange = range(x)
  hsv(0.7*(crange[2]-x)/(crange[2]-crange[1]),0.7,0.7)
}

## Load data

testdata <- read.table("testData123.csv", header=TRUE, sep=",")

# check your data import
colnames(testdata)
sapply(testdata, class)
head(testdata)
summary(testdata)

# convert to spatial dataset
ds = SpatialPointsDataFrame(testdata[,c("x","y")], testdata)

# spatial distribution
plot(ds, main="spatial distribution")

# variogram cloud of Co - complete area
vgmCloud_Co = variogram(Co~1, ds, cloud=TRUE)
plot(vgmCloud_Co, main="Variogram cloud of Co / complete area")

# variogram cloud of Ni - complete area
vgmCloud_Ni = variogram(Ni~1, ds, cloud=TRUE)
plot(vgmCloud_Ni, main="Variogram cloud of Ni / complete area")


# empirical variogram of Co 
vgmEmp_Co = variogram(Co~1, ds)
plot(vgmEmp_Co, main="empirical variogram of Co / complete area")
vgmEmp_Co
var(testdata$Co)
tail(vgmEmp_Co)

# empirical variogram of Ni
vgmEmp_Ni = variogram(Ni~1, ds)
plot(vgmEmp_Ni, main=" empirical variogram of Ni / complete area")
vgmEmp_Ni
var(testdata$Ni)
tail(vgmEmp_Ni)


# Variogram modelling for complete area

# Co
model_Co_initial = vgm(psill=16.5, model="Sph", range=1400, nugget=1.5)

#fit model for Co
model_Co_fit = fit.variogram(vgmEmp_Co, model_Co_initial)

# to find fitted parameters
model_Co_fit

# plot empirical variogram and fitted model for Co
plot(vgmEmp_Co, model_Co_fit, main="Fitted variogram model of Co - complete area")


# Ni
model_Ni_initial = vgm(psill=31, model="Sph", range=1500, nugget=2.5)

# fit model for Ni
model_Ni_fit = fit.variogram(vgmEmp_Ni, model_Ni_initial)

# to find fitted parameters
model_Ni_fit

# plot Empirical variogram and fitted model for Ni
plot(vgmEmp_Ni, model_Ni_fit, main="Fitted variogram model of Ni - complete area")



# Split the area into two parts
# western part
testdata_west = testdata[testdata$x <= 1000, ]

# eastern part
testdata_east = testdata[testdata$x > 1000, ]

# number of points
dim(testdata_west)
dim(testdata_east)

# spatial datasets
ds_west = SpatialPointsDataFrame(testdata_west[,c("x","y")], testdata_west)
ds_east = SpatialPointsDataFrame(testdata_east[,c("x","y")], testdata_east)

# plot both parts
plot(ds_west, main="Spatial distribution - western part")
plot(ds_east, main="Spatial distribution - eastern part")



# Empirical variograms for Co for both parts
vgmEmp_Co_west = variogram(Co~1, ds_west)
plot(vgmEmp_Co_west, plot.numbers=TRUE,
     main="Empirical variogram of Co - western part")

vgmEmp_Co_east = variogram(Co~1, ds_east)
plot(vgmEmp_Co_east, plot.numbers=TRUE,
     main="Empircal variogram of Co - eastern part")

vgmEmp_Co_west
vgmEmp_Co_east

# Empirical variograms for ni for both parts

vgmEmp_Ni_west = variogram(Ni~1, ds_west)
plot(vgmEmp_Ni_west, plot.numbers=TRUE,
     main="Empirical variogram of Ni - western part")

vgmEmp_Ni_east = variogram(Ni~1, ds_east)
plot(vgmEmp_Ni_east, plot.numbers=TRUE,
     main="Empirical variogram of Ni - eastern part")

vgmEmp_Ni_west
vgmEmp_Ni_east


# variogram models for Co - western and eastern parts
model_Co_west_initial = vgm(psill=6, model="Sph", range=800, nugget=9)
model_Co_west_fit = fit.variogram(vgmEmp_Co_west, model_Co_west_initial)
model_Co_west_fit
plot(vgmEmp_Co_west, model_Co_west_fit,
     main="Fitted variogram model of Co - western part")


model_Co_east_initial = vgm(psill=7, model="Sph", range=700, nugget=2)
model_Co_east_fit = fit.variogram(vgmEmp_Co_east, model_Co_east_initial)
model_Co_east_fit
plot(vgmEmp_Co_east, model_Co_east_fit,
     main="Fitted variogram model of Co - eastern part")


# variogram models for Ni - western and eastern parts
model_Ni_west_initial = vgm(psill=14, model="Sph", range=1000, nugget=13)
model_Ni_west_fit = fit.variogram(vgmEmp_Ni_west, model_Ni_west_initial)
model_Ni_west_fit
plot(vgmEmp_Ni_west, model_Ni_west_fit,
     main="Fitted variogram model of Ni - western part")


model_Ni_east_initial = vgm(psill=10, model="Sph", range=900, nugget=2.5)
model_Ni_east_fit = fit.variogram(vgmEmp_Ni_east, model_Ni_east_initial)
model_Ni_east_fit
plot(vgmEmp_Ni_east, model_Ni_east_fit,
     main="Fitted variogram model of Ni - eastern part")



# anisotrophy

# variogram map for Co - complete area
vgmMap_Co = variogram(Co~1, ds, cutoff=1500, width=100, map=TRUE)
plot(vgmMap_Co, main="Variogram map of Co - complete area")


# variogram map for Ni - complete area
vgmMap_Ni = variogram(Ni~1, ds, cutoff=1500, width=100, map=TRUE)
plot(vgmMap_Ni, main="Variogram map of Ni - complete area")


# direction variograms for Co
vgmDir_Co = variogram(Co~1, ds, alpha=c(0,45,90,135), tol.hor=22.5)
plot(vgmDir_Co, plot.numbers=TRUE,
     main="Directional variograms of Co - complete area")


# direction variograms for Ni
vgmDir_Ni = variogram(Ni~1, ds, alpha=c(0,45,90,135), tol.hor=22.5)
plot(vgmDir_Ni, plot.numbers=TRUE,
     main="Direction variograms of Ni - complete area")




# Directional variogram models for Co

vgmDir_Co_0 = vgmDir_Co[vgmDir_Co$dir.hor == 0, ]
vgmDir_Co_90 = vgmDir_Co[vgmDir_Co$dir.hor == 90, ]

model_Co_0_initial = vgm(psill=10, model="Sph", range=1300, nugget=2)
model_Co_0_fit = fit.variogram(vgmDir_Co_0, model_Co_0_initial)
model_Co_0_fit
plot(vgmDir_Co_0, model_Co_0_fit,
     main="Directional fitted variogram model of Co - 0 degree")

model_Co_90_initial = vgm(psill=20, model="Sph", range=1300, nugget=2)
model_Co_90_fit = fit.variogram(vgmDir_Co_90, model_Co_90_initial)
model_Co_90_fit
plot(vgmDir_Co_90, model_Co_90_fit,
     main="Directional fitted variogram model of Co - 90 degree")


# Directional variogram models for Ni

vgmDir_Ni_0 = vgmDir_Ni[vgmDir_Ni$dir.hor == 0, ]
vgmDir_Ni_90 = vgmDir_Ni[vgmDir_Ni$dir.hor == 90, ]

model_Ni_0_initial = vgm(psill=16, model="Sph", range=1300, nugget=3)
model_Ni_0_fit = fit.variogram(vgmDir_Ni_0, model_Ni_0_initial)
model_Ni_0_fit
plot(vgmDir_Ni_0, model_Ni_0_fit,
     main="Directional fitted variogram model of Ni - 0 degree")

model_Ni_90_initial = vgm(psill=36, model="Sph", range=1500, nugget=3.5)
model_Ni_90_fit = fit.variogram(vgmDir_Ni_90, model_Ni_90_initial)
model_Ni_90_fit
plot(vgmDir_Ni_90, model_Ni_90_fit,
     main="Directional fitted variogram model of Ni - 90 degree")
