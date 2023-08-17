# Analysis of PM10 data using the informed partition model
# March 3rd 2023 BYU office

setwd("~/Research/BYU/InformedPriorPartitions/analysis/")
library(salso)
library(drpm)
library(MCMCpack)
library(rmutil, warn.conflicts=FALSE)
library(pdfCluster)
library(mclust)

#
# This data is found in the gstat package
library(gstat)
data(DE_RB_2005)
dat <- DE_RB_2005

# Create ymat with columns corresponding to time rows stations
N <- length(dat@sp)
Tm <- 365
y <- matrix(NA, nrow=N, ncol=Tm)
for(i in 1:Tm){
  y[dat@index[dat@index[,2]==i,1], i] <- dat@data[dat@index[,2]==i,1]
}
# Try to create an average PM10 per month
year <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),
          rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),
          rep(11,30),rep(12,31))
week <- rep(1:52, each=7)
ymn <- t(apply(y, 1, function(x) tapply(x,year,mean, na.rm=TRUE)))
ymnwk <- t(apply(y[,-365], 1, function(x) tapply(x, week,mean, na.rm=TRUE)))

# Keep those times that are only missing one station
keep <- which(apply(apply(y, 1, is.na),1,sum) == 1)
c(1:365)[keep]
ysub <- y[, keep]

## Keep those that don't have any missing values when averageing over a month
rmve <- -c(4,16,25,27,30,43,52,59,69)
ysub1 <- ymnwk[rmve,]
ysub2 <- ymn[rmve,]

library(terra)
v <- vect(dat@sp@coords, crs="+init=epsg:32632 +proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
vt <- project(v, "+proj=longlat +datum=WGS84")
mon_loc <- data.frame(geom(vt)[, c("x", "y")][rmve,],
                      lab=as.character(1:60))

cex_val <- t((t(ysub2) - apply(ysub2, 2, min))/
               (apply(apply(ysub2,2,range),2,diff)))


# spatial coordinates
s_coords <- (mon_loc[,c(1,2)])

smn <- apply(s_coords,2,mean)
ssd <- apply(s_coords,2,sd)
s_std <- t((t(s_coords) - smn)/ssd)



ip_mat <- matrix(1, nrow=nrow(ysub2), ncol=5)

# initial partition 1 puts units into left and right regions
# based on coordinate 1 (x-axis)
ip_mat[s_std[,1] > median(s_std[,1]),1] <- 2
plot(s_std, col=ip_mat[,1])

# initial partition 2 puts units into up and down regions
# based on coordinate 2 (y-axis)
ip_mat[s_std[,2] > median(s_std[,2]), 2] <- 2
plot(s_std, col=ip_mat[,2])

# initial partition 3  puts units into four cartician quadrants
ip_mat[s_std[,1] > 0 & s_std[,2] > 0, 3] <- 2
ip_mat[s_std[,1] > 0 & s_std[,2] < 0, 3] <- 3
ip_mat[s_std[,1] < 0 & s_std[,2] < 0, 3] <- 4
plot(s_std, col=ip_mat[,3])

# fit model with initial partition that puts units into
# groups of nine using  x and y
ip_mat[s_std[,1] > 0.75 & s_std[,2] > -0.75 & s_std[,2] < 0.75, 4] <- 2
ip_mat[s_std[,1] > 0.75 & s_std[,2] < -0.75, 4] <- 3
ip_mat[s_std[,1] > -0.75 & s_std[,1] < 0.75 & s_std[,2] > 0.75, 4] <- 4
ip_mat[s_std[,1] > -0.75 & s_std[,1] < 0.75 & s_std[,2] > -0.75 & s_std[,2] < 0.75, 4] <- 5
ip_mat[s_std[,1] > -0.75 & s_std[,1] < 0.75 & s_std[,2] < -0.75, 4] <- 6
ip_mat[s_std[,1] < -0.75 & s_std[,2] > 0.75, 4] <- 7
ip_mat[s_std[,1] < -0.75 & s_std[,2] > -0.75 & s_std[,2] < 0.75, 4] <- 8
ip_mat[s_std[,1] < -0.75 & s_std[,2] < -0.75, 4] <- 9
plot(s_std, col=rainbow(9)[ip_mat[,4]])


# Plot of monitoring stations using
library(maps)
library(ggmap)
g_map <- get_map( 'germany', zoom=6, maptype="satellite")
gg0 <- ggmap(g_map) +
  #  geom_point(data=mon_loc, aes(x, y), col='orange') +
  geom_text(data=mon_loc, aes(x=x,y=y,label=lab), size=3,
            col=rainbow(9)[ip_mat[,4]]) +
  ylab("latitude") + xlab("longitude")
#ggsave("~/Research/BYU/InformedPriorPartitions/latex/figures/PM10_initial_partition.pdf",
#       width=7, height=7)s




#             # m0, s20, A,    At,  Al,  be
modelPriors <- c(0, 100, 2.5, 100,  5, 1)


niter <- 150000; nburn <- 50000; nthin <- 50
nout <- (niter - nburn)/nthin

aprior <- rbind(c(1,9))

set.seed(303)
# Fit model with no initial partition
drpm0 <- drpm_fit(y=ysub2[,1, drop=FALSE],
                  M=1,
                  initial_partition = NULL,
                  starting_alpha = 0.5,
                  unit_specific_alpha=FALSE,
                  time_specific_alpha=FALSE,
                  alpha_0 = FALSE, # FALSE -> update alpha
                  eta1_0 = TRUE,
                  phi1_0 = TRUE,
                  modelPriors=modelPriors,
                  alphaPriors = aprior, # This is not used here
                  draws=niter, burn=nburn, thin=nthin,verbose=TRUE)


# fit model with initial partition and global prior for alpha
drpm1 <- drpm_fit(y=ysub2[,1, drop=FALSE],
                  initial_partition = ip_mat[,4],
                  starting_alpha = 0.5,
                  unit_specific_alpha=FALSE,
                  time_specific_alpha=FALSE,
                  alpha_0 = FALSE, # FALSE -> update alpha
                  eta1_0 = TRUE,
                  phi1_0 = TRUE,
                  modelPriors=modelPriors,
                  alphaPriors = aprior,
                  draws=niter, burn=nburn, thin=nthin,verbose=TRUE)


# fit model with initial partition and local prior for alpha
aprior <- matrix(1, nrow=nrow(s_coords), ncol=2)
aprior[c(22,32,39,36),1] <- 9
aprior[-c(22,32,39,36),2] <- 9
drpm2 <- drpm_fit(y=ysub2[,1, drop=FALSE],
                  initial_partition = ip_mat[,4],
                  starting_alpha = 0.5,
                  unit_specific_alpha=TRUE,
                  time_specific_alpha=FALSE, # only one time point so doesn't matter
                  alpha_0 = FALSE, # FALSE -> update alpha
                  eta1_0 = TRUE,
                  phi1_0 = TRUE,
                  modelPriors=modelPriors,
                  alphaPriors = aprior,
                  draws=niter, burn=nburn, thin=nthin,verbose=TRUE)


lpml1 <- c(drpm0$lpml,
           drpm1$lpml,
           drpm2$lpml)

waic1 <- c(drpm0$waic,
           drpm1$waic,
           drpm2$waic)

# co-clustering probabilities
ccprob1 <- rbind(psm(t(drpm0$Si[1,,]))[c(22, 32, 36, 39), c(22, 32, 36, 39)][c(2,12)],
                 psm(t(drpm1$Si[1,,]))[c(22, 32, 36, 39), c(22, 32, 36, 39)][c(2,12)],
                 psm(t(drpm2$Si[1,,]))[c(22, 32, 36, 39), c(22, 32, 36, 39)][c(2,12)])

library(xtable)
xtable(cbind(lpml1, waic1, ccprob1), digits=2)



pest0 <- salso(t(drpm0$Si[1,,]))
pest1 <- salso(t(drpm1$Si[1,,]))
pest2 <- salso(t(drpm2$Si[1,,]))


adjustedRandIndex(pest0, pest1)
adjustedRandIndex(pest0, pest2)
adjustedRandIndex(pest1, pest2)



# Can I plot the coclustering probabilities
library(fields)

nmon <- nrow(ysub2)



# Figure that displays the coclustering probabilities for each model fit
pdf("~/Research/BYU/InformedPriorPartitions/latex/figures/co_clustering2.pdf",
    height=7, width=18)

#par(mfrow=c(1,4),oma=c( 0,0,0,0))

layout(matrix(c(1,2,3,4,4,4), nrow = 2, byrow = TRUE),
       height=c(10,1))
#layout.show(4)
fields.style()
#par(mar=c( 0.7,1.7,0.6,0.1)) # bottome, left, top, right
par(mar=c( 2.0, 2.0, 0.5, 2.0)) # bottome, left, top, right
ord <- order(salso(t(drpm0$Si[1,,])))
image(psm(t(drpm0$Si[1,,]))[ord, ord],
           axes=FALSE, xlab="", ylab="", col = tim.colors())
axis(1, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.70, las=3)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.70, las=2)


#par(mar=c( 0.7,1.7,0.6,0.1)) # bottome, left, top, right
par(mar=c( 2.0, 2.0, 0.5, 2.0)) # bottome, left, top, right
ord <- order(salso(t(drpm1$Si[1,,])))
image(psm(t(drpm1$Si[1,,]))[ord, ord],
           axes=FALSE, xlab="", ylab="", col = tim.colors())
axis(1, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.70, las=3)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.70, las=2)



#par(mar=c( 1.75,1.7,1.1,0.1)) # bottome, left, top, right
par(mar=c( 2.0, 2.0, 0.5, 2.0)) # bottome, left, top, right
ord <- order(salso(t(drpm2$Si[1,,])))
image(psm(t(drpm2$Si[1,,]))[ord, ord],
           axes=FALSE, xlab="", ylab="", col = tim.colors())
axis(1, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.70, las=3)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.70, las=2)

#par(mar=c( 2.0,10.0,0.5,1.1)) # bottome, left, top, right
image(psm(t(drpm2$Si[1,,]))[ord, ord], col="white", xaxt = "n", yaxt = "n", axes=FALSE)
image.plot(legend.only=TRUE, zlim=c(0,1), horizontal=TRUE, graphics.reset=TRUE,
           legend.mar=25, legend.width=6)

dev.off()




# Plot of four maps that display the initial partition,
# partition estimate from baseline model, that from a global alpha model
# and that from a unit local model.
mon_loc_pest <- rbind(mon_loc, mon_loc, mon_loc, mon_loc)
mon_loc_pest$model <- factor(rep(c('initial partition','baseline model','global model','unit local model'), each=nrow(mon_loc)),
                             levels=c('initial partition','baseline model','global model','unit local model'))
mon_loc_pest$pest <- c(ip_mat[,4], pest0, pest1, pest2)
mon_loc_pest$size <- rep((cex_val[,1]+0.5)*0.1, times=4)

  ggmap(g_map,
      maprange = FALSE,
      base_layer = ggplot(data = mon_loc_pest,
                                aes(x = x, y = y, col=as.factor(pest), size=(size+1)))) +
  #  geom_point(data=mon_loc, aes(x, y), col='orange') +
  geom_text(data=mon_loc_pest, aes(x=x,y=y,label=lab)) +
  ylab("latitude") + xlab("longitude") +
  facet_wrap(~model) + scale_size(range = c(1.75,4), guide = F) +
  theme(legend.position = "none")

  theme(plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))

ggsave("~/Research/BYU/InformedPriorPartitions/latex/figures/PM10_partition_estimates2.pdf",
       height=10, width=10)





# Now fit models with all twelve time points rather than one
set.seed(202)
# Fit model with no initial partition and global alpha, alpha ~ beta(a,b)
aprior <- rbind(c(1,9))
drpm01 <- drpm_fit(y=ysub2,
                  initial_partition = NULL, # no initial partition
                  starting_alpha = 0.5,
                  unit_specific_alpha=FALSE,
                  time_specific_alpha=FALSE,
                  alpha_0 = FALSE, # FALSE -> update alpha
                  eta1_0 = TRUE,
                  phi1_0 = TRUE,
                  modelPriors=modelPriors,
                  alphaPriors = aprior,
                  draws=niter, burn=nburn, thin=nthin,verbose=TRUE)



# Fit model with no initial partition and time-specific alpha, alpha_t ~ beta(a,b)
drpm02 <- drpm_fit(y=ysub2,
                   initial_partition = NULL, # no initial partition
                   starting_alpha = 0.5,
                   unit_specific_alpha=FALSE,
                   time_specific_alpha=TRUE,
                   alpha_0 = FALSE, # FALSE -> update alpha
                   eta1_0 = TRUE,
                   phi1_0 = TRUE,
                   modelPriors=modelPriors,
                   alphaPriors = aprior,
                   draws=niter, burn=nburn, thin=nthin,verbose=TRUE)




# fit model with initial partition and alpha \sim beta(a,b)
drpm11 <- drpm_fit(y=ysub2,
                  initial_partition = ip_mat[,4],
                  starting_alpha = 0.5,
                  unit_specific_alpha=FALSE,
                  time_specific_alpha=FALSE,
                  alpha_0 = FALSE, # FALSE -> update alpha
                  eta1_0 = TRUE,
                  phi1_0 = TRUE,
                  modelPriors=modelPriors,
                  alphaPriors = aprior,
                  draws=niter, burn=nburn, thin=nthin,verbose=TRUE)

# fit model with initial partition and alpha_t \sim beta(a,b)
drpm21 <- drpm_fit(y=ysub2,
                   initial_partition = ip_mat[,4],
                   starting_alpha = 0.5,
                   unit_specific_alpha=FALSE,
                   time_specific_alpha=TRUE,
                   alpha_0 = FALSE, # FALSE -> update alpha
                   eta1_0 = TRUE,
                   phi1_0 = TRUE,
                   modelPriors=modelPriors,
                   alphaPriors = aprior,
                   draws=niter, burn=nburn, thin=nthin,verbose=TRUE)


# fit model with initial partition and alpha_i \sim beta(a_i, b_i)
aprior <- matrix(1, nrow=nrow(s_coords), ncol=2)
aprior[c(22,32,39,36),1] <- 9
aprior[-c(22,32,39,36),2] <- 9
drpm31 <- drpm_fit(y=ysub2,
                  initial_partition = ip_mat[,4],
                  starting_alpha = 0.5,
                  unit_specific_alpha=TRUE,
                  time_specific_alpha=FALSE,
                  alpha_0 = FALSE, # FALSE -> update alpha
                  eta1_0 = TRUE,
                  phi1_0 = TRUE,
                  modelPriors=modelPriors,
                  alphaPriors = aprior,
                  draws=niter, burn=nburn, thin=nthin,verbose=TRUE)



# fit model with initial partition and alpha_it \sim beta(a_i, b_i)
set.seed(101)
drpm41 <- drpm_fit(y=ysub2,
                   initial_partition = ip_mat[,4],
                   starting_alpha = 0.5,
                   unit_specific_alpha=TRUE,
                   time_specific_alpha=TRUE,
                   alpha_0 = FALSE, # FALSE -> update alpha
                   eta1_0 = TRUE,
                   phi1_0 = TRUE,
                   modelPriors=modelPriors,
                   alphaPriors = aprior,
                   draws=niter, burn=nburn, thin=nthin,verbose=TRUE)

load("~/Research/BYU/InformedPriorPartitions/analysis/modelFits_PM10.RData")


lpml2 <- c(drpm01$lpml,
drpm02$lpml,
drpm11$lpml,
drpm21$lpml,
drpm31$lpml,
drpm41$lpml)

waic2 <- c(drpm01$waic,
drpm02$waic,
drpm11$waic,
drpm21$waic,
drpm31$waic,
drpm41$waic)

library(xtable)
xtable(cbind(lpml2[-2], waic2[-2]), digits=2)


# MSE values (not included in the document)
tmp01 <- t(apply(drpm01$fitted,c(1,2),mean))
tmp02 <- t(apply(drpm02$fitted,c(1,2),mean))
tmp11 <- t(apply(drpm11$fitted,c(1,2),mean))
tmp21 <- t(apply(drpm21$fitted,c(1,2),mean))
tmp31 <- t(apply(drpm31$fitted,c(1,2),mean))
tmp41 <- t(apply(drpm41$fitted,c(1,2),mean))

c(mean((tmp01 - ysub2)^2),
mean((tmp02 - ysub2)^2),
mean((tmp11 - ysub2)^2),
mean((tmp21 - ysub2)^2),
mean((tmp31 - ysub2)^2),
mean((tmp41 - ysub2)^2))

# Plot fitted vs observed values (not included in the document)
plot(ysub2[,1], tmp01[,1], ylim=range(ysub2))
points(ysub2[,1], tmp02[,1], col='red')
points(ysub2[,1], tmp11[,1], col='brown')
points(ysub2[,1], tmp21[,1], col='orange')
points(ysub2[,1], tmp31[,1], col='blue')
points(ysub2[,1], tmp41[,1], col='green')


# Compare partition estimates from model with one time period to those with
# 12.
adjustedRandIndex(salso(t(drpm1$Si[1,,])), salso(t(drpm11$Si[1,,])))
adjustedRandIndex(salso(t(drpm1$Si[1,,])), salso(t(drpm01$Si[1,,])))

adjustedRandIndex(salso(t(drpm1$Si[1,,])), ip_mat[,4])
adjustedRandIndex(salso(t(drpm11$Si[1,,])), ip_mat[,4])
adjustedRandIndex(salso(t(drpm01$Si[1,,])), ip_mat[,4])


plot(s_coords, col=salso(t(drpm1$Si[1,,])), cex=cex_val[,12]+0.5)
plot(s_coords, col=salso(t(drpm11$Si[1,,])), cex=cex_val[,12]+0.5)
plot(s_coords, col=salso(t(drpm01$Si[1,,])), cex=cex_val[,12]+0.5)

# Plot of partition estimates time period 1 and 12 (not included in document)
par(mfrow=c(2,3))
plot(mon_loc[,-3], cex=cex_val[,1]+0.5, pch=19, col=salso(t(drpm01$S[1,,])))
plot(mon_loc[,-3], cex=cex_val[,1]+0.5, pch=19, col=salso(t(drpm02$S[1,,])))
plot(mon_loc[,-3], cex=cex_val[,1]+0.5, pch=19, col=salso(t(drpm11$S[1,,])))
plot(mon_loc[,-3], cex=cex_val[,1]+0.5, pch=19, col=salso(t(drpm21$S[1,,])))
plot(mon_loc[,-3], cex=cex_val[,1]+0.5, pch=19, col=salso(t(drpm31$S[1,,])))
plot(mon_loc[,-3], cex=cex_val[,1]+0.5, pch=19, col=salso(t(drpm41$S[1,,])))

plot(mon_loc[,-3], cex=cex_val[,12]+0.5, pch=19, col=salso(t(drpm01$S[12,,])))
plot(mon_loc[,-3], cex=cex_val[,12]+0.5, pch=19, col=salso(t(drpm02$S[12,,])))
plot(mon_loc[,-3], cex=cex_val[,12]+0.5, pch=19, col=salso(t(drpm11$S[12,,])))
plot(mon_loc[,-3], cex=cex_val[,12]+0.5, pch=19, col=salso(t(drpm21$S[12,,])))
plot(mon_loc[,-3], cex=cex_val[,12]+0.5, pch=19, col=salso(t(drpm31$S[12,,])))
plot(mon_loc[,-3], cex=cex_val[,12]+0.5, pch=19, col=salso(t(drpm41$S[12,,])))

# plot alpha
alpha01 <- apply(drpm01$alpha,2,mean)[2]
alpha02 <- apply(drpm02$alpha,2,mean)[-1]
alpha11 <- apply(drpm11$alpha,2,mean)[1]
alpha21 <- apply(drpm21$alpha,2,mean)
alpha31 <- apply(drpm31$alpha[1,,],1,mean )
alpha41 <- apply(drpm41$alpha, c(1,2), mean)

alpha_mat <- cbind(rbind(alpha41, alpha31), c(alpha21,NA))
image.plot(alpha_mat, xlab="time", ylab='station',axes=FALSE, zlim=c(0,1))
axis(1, at = c(0:(13-1))/(13-1), labels=c(1:12,NA), tick=TRUE, cex.axis=0.70, las=1)
axis(2, at = c(0:(nmon))/(nmon), labels=c(1:nmon, NA), tick=TRUE, cex.axis=0.70, las=2)

image.plot(alpha41, xlab="time", ylab='station',axes=FALSE, zlim=c(0,1))
axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon), tick=TRUE, cex.axis=0.70, las=2)



# Plot co-clustering probabilities for a stations 22, 32 and 36, 39
ccp <- matrix(NA, 12, 10)
for(jj in 1:12){
  ccp[jj, c(1,2)] <- diag(psm(t(drpm01$Si[jj,,]))[c(22, 36), c(32,39)])
  ccp[jj, c(3,4)] <- diag(psm(t(drpm11$Si[jj,,]))[c(22, 36), c(32,39)])
  ccp[jj, c(5,6)] <- diag(psm(t(drpm21$Si[jj,,]))[c(22, 36), c(32,39)])
  ccp[jj, c(7,8)] <- diag(psm(t(drpm31$Si[jj,,]))[c(22, 36), c(32,39)])
  ccp[jj, c(9,10)] <- diag(psm(t(drpm41$Si[jj,,]))[c(22, 36), c(32,39)])
}


pdf("~/Research/BYU/InformedPriorPartitions/latex/figures/co_clustering_all_time.pdf",
      height=7, width=9)
  par(mfrow=c(1,1))
  par(mar = c(3.5, 6, 2, 6.5))
  image.plot(ccp, xlab="", ylab='',axes=FALSE, zlim=c(0,1), horizontal=TRUE)
  axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.9, las=1)
  axis(2, at = c(0:(10-1))/(10-1),
        labels=rep(c(expression(Pr(c[22]==c[32])),expression(Pr(c[36]==c[39]))), times=5),
        tick=TRUE, cex.axis=0.9, las=2)
  axis(4, at = c(0.06, 0.28, 0.5, 0.72, 0.94),
       labels=c("baseline","global","unit local","time local","unit x time local"),
       tick=TRUE, cex.axis=0.9, las=2)
  mtext("time",side=1, line=2)
dev.off()



relabel_labels <- function(z) {
  # Get the unique labels in the input vector
  unique_labels <- unique(z)

  # Create a mapping from old labels to new labels
  old_to_new_labels <- 1:length(unique_labels)
  names(old_to_new_labels) <- unique_labels

  # Use the match() function to update labels
  new_labels <- old_to_new_labels[match(z, unique_labels)]

  return(new_labels)
}

library(salso)
pdf("~/Research/BYU/InformedPriorPartitions/latex/figures/station22_32_36_39.pdf",
    height=8, width=8)
clrs <- c("red","blue", "green","orange","purple","turquoise","magenta","orangered2", "cyan")

                      # b,l,t,r
par(mfrow=c(2,2), mar=c(4,4,1,1))

matplot(t(ysub2), type='p', col='gray95', lty=1, pch=1, cex=0.5,
        xlab="month", ylab=expression('PM'[10]))

for(jj in 1:12){
  ppp <- relabel_labels(salso(t(drpm01$Si[jj,,])))
  text(rep(jj,4), ysub2[c(22,32,36,39),jj], labels=c(22,32,36,39),cex=0.75,
       col=clrs[ppp[c(22,32,36,39)]])
  points(rep(jj,length(ppp)-4), ysub2[-c(22,32,36,39),jj],
         col=clrs[ppp[-c(22,32,36,39)]], cex=0.25)
}
#legend(x="topleft", legend=c("cluster1","cluster2"), col=c("red","blue"),
#       pch=1,ncol=2)

matplot(t(ysub2), type='p', col='gray95', lty=1, pch=1, cex=0.5,
        xlab="month", ylab=expression('PM'[10]))
for(jj in 1:12){
  ppp <- relabel_labels(salso(t(drpm11$Si[jj,,])))
  text(rep(jj,4), ysub2[c(22,32,36,39),jj], labels=c(22,32,36,39),cex=0.75,
       col=clrs[ppp[c(22,32,36,39)]])
  points(rep(jj,length(ppp)-4), ysub2[-c(22,32,36,39),jj],
         col=clrs[ppp[-c(22,32,36,39)]], cex=0.25)

}

matplot(t(ysub2), type='p', col='gray95', lty=1, pch=1, cex=0.5,
        xlab="month", ylab=expression('PM'[10]))
for(jj in 1:12){
  ppp <- relabel_labels(salso(t(drpm21$Si[jj,,])))
  text(rep(jj,4), ysub2[c(22,32,36,39),jj], labels=c(22,32,36,39),cex=0.75,
       col=clrs[ppp[c(22,32,36,39)]])
  points(rep(jj,length(ppp)-4), ysub2[-c(22,32,36,39),jj],
         col=clrs[ppp[-c(22,32,36,39)]], cex=0.25)

}

matplot(t(ysub2), type='p', col='gray95', lty=1, pch=1, cex=0.5,
        xlab="month", ylab=expression('PM'[10]))
for(jj in 1:12){
  ppp <- relabel_labels(salso(t(drpm31$Si[jj,,])))
  text(rep(jj,4), ysub2[c(22,32,36,39),jj], labels=c(22,32,36,39),cex=0.75,
       col=clrs[ppp[c(22,32,36,39)]])
  points(rep(jj,length(ppp)-4), ysub2[-c(22,32,36,39),jj],
         col=clrs[ppp[-c(22,32,36,39)]], cex=0.25)

}

#legend(x="topleft", legend=c("cluster1","cluster2","cluster3"), col=c("red","blue","green"),
#       pch=1,ncol=3)
dev.off()



par(mfrow=c(1,2))
plot(1:12,ccp[,2], type='b', ylim=c(0,1))
points(1:12,ccp[,4], type='b', col='brown')
points(1:12,ccp[,6], type='b', col='orange')
points(1:12,ccp[,8], type='b', col='green')
points(1:12,ccp[,10], type='b', col='blue')

plot(1:12,ccp[,1], type='b', ylim=c(0,1))
points(1:12,ccp[,3], type='b', col='brown')
points(1:12,ccp[,5], type='b', col='orange')
points(1:12,ccp[,7], type='b', col='green')
points(1:12,ccp[,9], type='b', col='blue')



par(mfrow=c(4,3))
for(jj in 1:12){
  ord <- order(salso(t(drpm11$Si[jj,,])))
  image(psm(t(drpm11$Si[jj,,]))[ord, ord],
        axes=FALSE, xlab="", ylab="", col = tim.colors())
  axis(1, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.60, las=2)
  axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.60, las=2)

}


par(mfrow=c(4,3))
for(jj in 1:12){
  ord <- order(ip_mat[,4])
  image(psm(t(drpm01$Si[jj,,]))[ord, ord],
        axes=FALSE, xlab="", ylab="", col = tim.colors())
  axis(1, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.60, las=2)
  axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon)[ord], tick=TRUE, cex.axis=0.60, las=2)

}


cp590 <- matrix(NA, 12, nrow(ysub2))
cp591 <- matrix(NA, 12, nrow(ysub2))
cp592 <- matrix(NA, 12, nrow(ysub2))
cp593 <- matrix(NA, 12, nrow(ysub2))
cp594 <- matrix(NA, 12, nrow(ysub2))
station <- 25
for(jj in 1:12){
  cp590[jj,] <- psm(t(drpm01$Si[jj,,]))[station,]
  cp591[jj,] <- psm(t(drpm11$Si[jj,,]))[station,]
  cp592[jj,] <- psm(t(drpm21$Si[jj,,]))[station,]
  cp593[jj,] <- psm(t(drpm31$Si[jj,,]))[station,]
  cp594[jj,] <- psm(t(drpm41$Si[jj,,]))[station,]
}
par(mfrow=c(2,3))
image.plot(cp590, xlab="time", ylab='station',axes=FALSE, zlim=c(0,1), main="baseline model")
axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon), tick=TRUE, cex.axis=0.70, las=2)

image.plot(cp591, xlab="time", ylab='station',axes=FALSE, zlim=c(0,1), main="global model")
axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon), tick=TRUE, cex.axis=0.70, las=2)

image.plot(cp592, xlab="time", ylab='station',axes=FALSE, zlim=c(0,1), main="time local model")
axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon), tick=TRUE, cex.axis=0.70, las=2)

image.plot(cp593, xlab="time", ylab='station',axes=FALSE, zlim=c(0,1), main="unit local model")
axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon), tick=TRUE, cex.axis=0.70, las=2)

image.plot(cp594, xlab="time", ylab='station',axes=FALSE, zlim=c(0,1), main="unit x time local model")
axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
axis(2, at = c(0:(nmon-1))/(nmon-1), labels=c(1:nmon), tick=TRUE, cex.axis=0.70, las=2)


# Find the adjusted rand index based on the initial partitions
ari_mat <- matrix(NA, 12, 6)
for(jj in 1:12){
  cat("jj = ", jj, "\n")
  mn01 <- mn02 <- mn1 <- mn2 <- mn3 <- mn4 <- 0
  for(kk in 1:nout){
    mn01 <- mn01 + adjustedRandIndex(ip_mat[,4], drpm01$Si[jj,,kk])/nout
    mn02 <- mn02 + adjustedRandIndex(ip_mat[,4], drpm02$Si[jj,,kk])/nout
    mn1 <- mn1 + adjustedRandIndex(ip_mat[,4], drpm11$Si[jj,,kk])/nout
    mn2 <- mn2 + adjustedRandIndex(ip_mat[,4], drpm21$Si[jj,,kk])/nout
    mn3 <- mn3 + adjustedRandIndex(ip_mat[,4], drpm31$Si[jj,,kk])/nout
    mn4 <- mn4 + adjustedRandIndex(ip_mat[,4], drpm41$Si[jj,,kk])/nout
  }
  ari_mat[jj,1] <- mn01
  ari_mat[jj,2] <- mn02
  ari_mat[jj,3] <- mn1
  ari_mat[jj,4] <- mn2
  ari_mat[jj,5] <- mn3
  ari_mat[jj,6] <- mn4
}


par(mfrow=c(2,3))
plot(ari_mat[,1], pch=19, ylim=range(ari_mat), ylab="ARI", xlab="time")
plot(ari_mat[,2], pch=19, ylim=range(ari_mat), ylab="ARI", xlab="time")
plot(ari_mat[,3], pch=19, ylim=range(ari_mat), ylab="ARI", xlab="time")
plot(ari_mat[,4], pch=19, ylim=range(ari_mat), ylab="ARI", xlab="time")
plot(ari_mat[,5], pch=19, ylim=range(ari_mat), ylab="ARI", xlab="time")
plot(ari_mat[,6], pch=19, ylim=range(ari_mat), ylab="ARI", xlab="time")


# Find the lagged adjusted rand index
ari_mat01 <- matrix(0, 12, 12)
ari_mat02 <- matrix(0, 12, 12)
ari_mat1 <- matrix(0, 12, 12)
ari_mat2 <- matrix(0, 12, 12)
ari_mat3 <- matrix(0, 12, 12)
ari_mat4 <- matrix(0, 12, 12)
for(jj in 1:12){
  cat("jj = ", jj, "\n")
  for(ii in 1:12){
    mn01 <- mn02 <- mn1 <- mn2 <- mn3 <- mn4 <- 0
    for(kk in 1:nout){
      mn01 <- mn01 + adjustedRandIndex(drpm01$Si[jj,,kk], drpm01$Si[ii,,kk])/nout
      mn02 <- mn02 + adjustedRandIndex(drpm02$Si[jj,,kk], drpm02$Si[ii,,kk])/nout
      mn1 <- mn1 + adjustedRandIndex(drpm11$Si[jj,,kk], drpm11$Si[ii,,kk])/nout
      mn2 <- mn2 + adjustedRandIndex(drpm21$Si[jj,,kk], drpm21$Si[ii,,kk])/nout
      mn3 <- mn3 + adjustedRandIndex(drpm31$Si[jj,,kk], drpm31$Si[ii,,kk])/nout
      mn4 <- mn4 + adjustedRandIndex(drpm41$Si[jj,,kk], drpm41$Si[ii,,kk])/nout
    }
    ari_mat01[jj,ii] <- mn01
    ari_mat02[jj,ii] <- mn02
    ari_mat1[jj,ii] <- mn1
    ari_mat2[jj,ii] <- mn2
    ari_mat3[jj,ii] <- mn3
    ari_mat4[jj,ii] <- mn4
  }
}

library(fields)
pdf("~/Research/BYU/InformedPriorPartitions/latex/figures/ARI.pdf", height=6, width=9.5)
  par(mfrow=c(2,3))
  par(mar = c(3, 3, 3, 5))

  image.plot(ari_mat01, xlab="", ylab='',axes=FALSE, zlim=c(-0.1,1), main="baseline")
#  mtext("time", side=c(2), line=1.5);  mtext("time", side=c(1), line=2.25);
  axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
  axis(2, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=2)

  image.plot(ari_mat1, xlab="", ylab='',axes=FALSE, zlim=c(-0.1,1), main="global")
#    mtext("time", side=c(2), line=1.5);  mtext("time", side=c(1), line=2.25);
  axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
  axis(2, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=2)

  image.plot(ari_mat2, xlab="", ylab='',axes=FALSE, zlim=c(-0.1,1), main="unit local")
#  mtext("time", side=c(2), line=1.5);  mtext("time", side=c(1), line=2.25);
  axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
  axis(2, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=2)

  image.plot(ari_mat3, xlab="", ylab='',axes=FALSE, zlim=c(-0.1,1), main="time local")
#  mtext("time", side=c(2), line=1.5);  mtext("time", side=c(1), line=2.25);
  axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
  axis(2, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=2)

  image.plot(ari_mat4, xlab="", ylab='',axes=FALSE, zlim=c(-0.1,1), main="unit x time local")
#  mtext("time", side=c(2), line=1.5);  mtext("time", side=c(1), line=2.25);
  axis(1, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=1)
  axis(2, at = c(0:(12-1))/(12-1), labels=c(1:12), tick=TRUE, cex.axis=0.70, las=2)
dev.off()


# install.packages("reshape")
library(reshape)

# Data
set.seed(8)
m <- matrix(round(rnorm(200), 2), 10, 10)
colnames(m) <- paste("Col", 1:10)
rownames(m) <- paste("Row", 1:10)

# Transform the matrix in long format
ari_mat01_df <- melt(ari_mat01)
gg1 <- ggplot(ari_mat01_df, aes(x = X1, y = X2, fill = value)) +
  geom_tile()

ari_mat1_df <- melt(ari_mat1)
gg2 <- ggplot(ari_mat1_df, aes(x = X1, y = X2)) +
  geom_tile(aes(fill=value)) + scale_colour_gradientn(colours=tim.colors(64))

ari_mat2_df <- melt(ari_mat2)
gg3 <- ggplot(ari_mat2_df, aes(x = X1, y = X2, fill = value)) +
  geom_tile()

apply(drpm01$alpha,2,mean)[1]
apply(drpm02$alpha,2,mean)[1]
apply(drpm11$alpha,2,mean)[2]
apply(drpm21$alpha,2,mean)
apply(drpm31$alpha,c(1,2),mean)[1,]
apply(drpm41$alpha,c(1,2),mean)



# number of clusters over time.  This is included in the model.
par(mfrow=c(2,3))
plot(apply(sapply(1:12, function(x) apply(t(drpm01$Si[x,,]), 1, function(y) length(unique(y)))),2,mean), type='b', ylab="number of clusters", xlab="time", main="model0", ylim=c(1,9))
plot(apply(sapply(1:12, function(x) apply(t(drpm02$Si[x,,]), 1, function(y) length(unique(y)))),2,mean), type='b', ylab="number of clusters", xlab="time", main="model0", ylim=c(1,9))
plot(apply(sapply(1:12, function(x) apply(t(drpm11$Si[x,,]), 1, function(y) length(unique(y)))),2,mean), type='b', ylab="number of clusters", xlab="time", main="model3", ylim=c(1,9))
plot(apply(sapply(1:12, function(x) apply(t(drpm21$Si[x,,]), 1, function(y) length(unique(y)))),2,mean), type='b', ylab="number of clusters", xlab="time", main="model4", ylim=c(1,9))
plot(apply(sapply(1:12, function(x) apply(t(drpm31$Si[x,,]), 1, function(y) length(unique(y)))),2,mean), type='b', ylab="number of clusters", xlab="time", main="model5", ylim=c(1,9))
plot(apply(sapply(1:12, function(x) apply(t(drpm41$Si[x,,]), 1, function(y) length(unique(y)))),2,mean), type='b', ylab="number of clusters", xlab="time", main="model6", ylim=c(1,9))

pdf("~/Research/BYU/InformedPriorPartitions/latex/figures/number_of_clusters.pdf", height=8, width=8)
  par(mfrow=c(1,1))
  par(mar = c(4, 4.75, 2, 2))
  plot(apply(sapply(1:12, function(x) apply(t(drpm01$Si[x,,]), 1, function(y) length(unique(y)))),2,mean),
       type='b', ylab="number of clusters", xlab="time",  ylim=c(1,9), pch=19, cex.lab=2, cex.axis=1.5)
#points(apply(sapply(1:12, function(x) apply(t(drpm02$Si[x,,]), 1, function(y) length(unique(y)))),2,mean),
#       type='b', col='red')
  points(apply(sapply(1:12, function(x) apply(t(drpm11$Si[x,,]), 1, function(y) length(unique(y)))),2,mean),
         type='b', col='brown', pch=19)
  points(apply(sapply(1:12, function(x) apply(t(drpm21$Si[x,,]), 1, function(y) length(unique(y)))),2,mean),
         type='b', col='orange', pch=19)
  points(apply(sapply(1:12, function(x) apply(t(drpm31$Si[x,,]), 1, function(y) length(unique(y)))),2,mean),
         type='b', col='green', pch=19)
  points(apply(sapply(1:12, function(x) apply(t(drpm41$Si[x,,]), 1, function(y) length(unique(y)))),2,mean),
         type='b', col='blue', pch=19)
  legend(x='topright', legend=c("baseline", "global", "unit local","time local","unit x time local"),
         col=c("black","brown","orange","green","blue"), lty=1,  pch=19, cex=1.75)
dev.off()


# Old stuff that isn't included in the main document
# Plot the estimated partitions using germany

gg1 <- ggmap(g_map) +
  #  geom_point(data=mon_loc, aes(x, y), col='orange') +
  geom_text(data=mon_loc, aes(x=x,y=y,label=lab), size=3,
            col=rainbow(length(unique(pest0)))[pest0]) +
  ylab("latitude") + xlab("longitude") +
  theme(plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))

gg2 <- ggmap(g_map) +
  #  geom_point(data=mon_loc, aes(x, y), col='orange') +
  geom_text(data=mon_loc, aes(x=x,y=y,label=lab), size=3,
            col=rainbow(length(unique(pest1)))[pest1]) +
  ylab("latitude") + xlab("longitude") +
  theme(plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))

gg3 <- ggmap(g_map) +
  #  geom_point(data=mon_loc, aes(x, y), col='orange') +
  geom_text(data=mon_loc, aes(x=x,y=y,label=lab), size=3,
            col=rainbow(length(unique(pest2)))[pest2]) +
  ylab("latitude") + xlab("longitude") +
  theme(plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))

library(patchwork)
(gg0 + gg1)/(gg2 + gg3)


