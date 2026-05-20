#Wilcoxon Rank test

w <- all_merged3_scenerio0[(all_merged3_scenerio0$Name=="WISpR"),]
spot <- all_merged3_scenerio0[(all_merged3_scenerio0$Name=="SPOTlight"),]
ster <- all_merged3_scenerio0[(all_merged3_scenerio0$Name=="Stereoscope"),]
dw <- all_merged3_scenerio0[(all_merged3_scenerio0$Name=="DWLS"),]
sdw <- all_merged3_scenerio0[(all_merged3_scenerio0$Name=="S-DWLS"),]
rc <- all_merged3_scenerio0[(all_merged3_scenerio0$Name=="RCTD"),]


w2 <- all_merged3_scenerio1[(all_merged3_scenerio1$Name=="WISpR"),]
spot2 <- all_merged3_scenerio1[(all_merged3_scenerio1$Name=="SPOTlight"),]
ster2 <- all_merged3_scenerio1[(all_merged3_scenerio1$Name=="Stereoscope"),]
dw2 <- all_merged3_scenerio1[(all_merged3_scenerio1$Name=="DWLS"),]
sdw2 <- all_merged3_scenerio1[(all_merged3_scenerio1$Name=="S-DWLS"),]
rc2 <- all_merged3_scenerio1[(all_merged3_scenerio1$Name=="RCTD"),]


w3 <- all_merged3_scenerio2[(all_merged3_scenerio2$Name=="WISpR"),]
spot3 <- all_merged3_scenerio2[(all_merged3_scenerio2$Name=="SPOTlight"),]
ster3 <- all_merged3_scenerio2[(all_merged3_scenerio2$Name=="Stereoscope"),]
dw3 <- all_merged3_scenerio2[(all_merged3_scenerio2$Name=="DWLS"),]
sdw3 <- all_merged3_scenerio2[(all_merged3_scenerio2$Name=="S-DWLS"),]
rc3 <- all_merged3_scenerio2[(all_merged3_scenerio2$Name=="RCTD"),]


w4 <- all_merged3_scenerio3[(all_merged3_scenerio3$Name=="WISpR"),]
spot4 <- all_merged3_scenerio3[(all_merged3_scenerio3$Name=="SPOTlight"),]
ster4 <- all_merged3_scenerio3[(all_merged3_scenerio3$Name=="Stereoscope"),]
dw4 <- all_merged3_scenerio3[(all_merged3_scenerio3$Name=="DWLS"),]
sdw4 <- all_merged3_scenerio3[(all_merged3_scenerio3$Name=="S-DWLS"),]
rc4 <- all_merged3_scenerio3[(all_merged3_scenerio3$Name=="RCTD"),]

w_spot <- wilcox.test(x= w$RMSE, y =spot$RMSE,
                      alternative = c("less"),
                      mu = 0, paired = TRUE,
                      conf.int = TRUE)

w_ster <- wilcox.test(x= w$RMSE, y =ster$RMSE,
                      alternative = c("less"),
                      mu = 0, paired = TRUE,
                      conf.int = TRUE)

w_dw <- wilcox.test(x= w$RMSE, y =dw$RMSE,
                    alternative = c("less"),
                    mu = 0, paired = TRUE,
                    conf.int = TRUE)

w_sdw <- wilcox.test(x= w$RMSE, y =sdw$RMSE,
                     alternative = c("less"),
                     mu = 0, paired = TRUE,
                     conf.int = TRUE)

w_rc <- wilcox.test(x= w$RMSE, y =rc$RMSE,
                    alternative = c("less"),
                    mu = 0, paired = TRUE,
                    conf.int = TRUE)
library(moments)

skewness(w$RMSE); skewness(dw$RMSE); skewness(rc$RMSE); skewness(sdw$RMSE); skewness(ster$RMSE); skewness(spot$RMSE)
skewness(w2$RMSE); skewness(dw2$RMSE); skewness(rc2$RMSE); skewness(sdw2$RMSE); skewness(ster2$RMSE); skewness(spot2$RMSE)
skewness(w3$RMSE); skewness(dw3$RMSE); skewness(rc3$RMSE); skewness(sdw3$RMSE); skewness(ster3$RMSE); skewness(spot3$RMSE)
skewness(w4$RMSE); skewness(dw4$RMSE); skewness(rc4$RMSE); skewness(sdw4$RMSE); skewness(ster4$RMSE); skewness(spot4$RMSE)

par(mfrow = c(6, 4), mar=c(1,1,1,1))
# Add multiple vertical lines at specific locations with different colors

p1 <- hist(w$RMSE, breaks = seq(min(w$RMSE), max(w$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.032), col = c('red'), lwd = 1.5, lty = 'dashed')
p7 <- hist(w2$RMSE, breaks = seq(min(w2$RMSE), max(w2$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.050), col = c('red'), lwd = 1.5, lty = 'dashed')
p13 <- hist(w3$RMSE, breaks = seq(min(w3$RMSE), max(w3$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.071), col = c('red'), lwd = 1.5, lty = 'dashed')
p19 <- hist(w4$RMSE, breaks = seq(min(w4$RMSE), max(w4$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.057), col = c('red'), lwd = 1.5, lty = 'dashed')

p2 <- hist(dw$RMSE, breaks = seq(min(dw$RMSE), max(dw$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.038), col = c('red'), lwd = 1.5, lty = 'dashed')
p8 <- hist(dw2$RMSE, breaks = seq(min(dw2$RMSE), max(dw2$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.058), col = c('red'), lwd = 1.5, lty = 'dashed')
p14 <- hist(dw3$RMSE, breaks = seq(min(dw3$RMSE), max(dw3$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.071), col = c('red'), lwd = 1.5, lty = 'dashed')
p20 <- hist(dw4$RMSE, breaks = seq(min(dw4$RMSE), max(dw4$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.062), col = c('red'), lwd = 1.5, lty = 'dashed')

p3 <- hist(rc$RMSE, breaks = seq(min(rc$RMSE), max(rc$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL) 
abline(v = c(0.052), col = c('red'), lwd = 1.5, lty = 'dashed')
p9 <- hist(rc2$RMSE, breaks = seq(min(rc2$RMSE), max(rc2$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL) 
abline(v = c(0.123), col = c('red'), lwd = 1.5, lty = 'dashed')
p15 <- hist(rc3$RMSE, breaks = seq(min(rc3$RMSE), max(rc3$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL) 
abline(v = c(0.068), col = c('red'), lwd = 1.5, lty = 'dashed')
p21 <- hist(rc4$RMSE, breaks = seq(min(rc4$RMSE), max(rc4$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL) 
abline(v = c(0.065), col = c('red'), lwd = 1.5, lty = 'dashed')

p4 <- hist(sdw$RMSE, breaks = seq(min(sdw$RMSE), max(sdw$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.046), col = c('red'), lwd = 1.5, lty = 'dashed')
p10 <- hist(sdw2$RMSE, breaks = seq(min(sdw2$RMSE), max(sdw2$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.059), col = c('red'), lwd = 1.5, lty = 'dashed')
p16 <- hist(sdw3$RMSE, breaks = seq(min(sdw3$RMSE), max(sdw3$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.073), col = c('red'), lwd = 1.5, lty = 'dashed')
p22 <- hist(sdw4$RMSE, breaks = seq(min(sdw4$RMSE), max(sdw4$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.067), col = c('red'), lwd = 1.5, lty = 'dashed')

p6 <- hist(spot$RMSE, breaks = seq(min(spot$RMSE), max(spot$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.233), col = c('red'), lwd = 1.5, lty = 'dashed')
p12 <- hist(spot2$RMSE, breaks = seq(min(spot2$RMSE), max(spot2$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.260), col = c('red'), lwd = 1.5, lty = 'dashed')
p18 <- hist(spot3$RMSE, breaks = seq(min(spot3$RMSE), max(spot3$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.272), col = c('red'), lwd = 1.5, lty = 'dashed')
p24 <- hist(spot4$RMSE, breaks = seq(min(spot4$RMSE), max(spot4$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.261), col = c('red'), lwd = 1.5, lty = 'dashed')

p5 <- hist(ster$RMSE, breaks = seq(min(ster$RMSE), max(ster$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.046), col = c('red'), lwd = 1.5, lty = 'dashed')
p11 <- hist(ster2$RMSE, breaks = seq(min(ster2$RMSE), max(ster2$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.201), col = c('red'), lwd = 1.5, lty = 'dashed')
p17 <- hist(ster3$RMSE, breaks = seq(min(ster3$RMSE), max(ster3$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.210), col = c('red'), lwd = 1.5, lty = 'dashed')
p23 <- hist(ster4$RMSE, breaks = seq(min(ster4$RMSE), max(ster4$RMSE), length.out = 21), col=rgb(0,0,1,1/4), xlim=c(0,0.4), ylim= c(0,300), main = NULL)
abline(v = c(0.208), col = c('red'), lwd = 1.5, lty = 'dashed')

