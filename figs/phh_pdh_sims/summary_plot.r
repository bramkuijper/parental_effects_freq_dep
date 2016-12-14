library("lattice")
library("colorRamps")

#the.data <- read.table("all_sims.csv",sep=";",header=T)
the.data <- read.table("summary_sims_phh_pdh.csv",sep=";",header=T)

the.data[,"vc"] <- the.data$v / the.data$c

the.data[,"vc"] <- ifelse(the.data[,"vc"] > 1.0, 1.0, ifelse(the.data[,"vc"] < 0, 0, the.data[,"vc"]))

pdf("summary_sims_phh.pdf",width=10,height=7)
print(wireframe(meanphh ~ c * v | d,
                #        data=the.data[the.data$v < the.data$c,],
                data=the.data,
                pch=21,
                cex=0.2,
                strip=function(strip.levels,...) { 
                    strip.default(strip.levels=T, ...) },
                default.scales=list(arrows=F)
                )
        )
dev.off()

pdf("summary_sims_pdh.pdf",width=10,height=7)
print(wireframe(meanpdh ~ c * v | d,
                data=the.data,
                #                data=the.data,
                pch=21,
                cex=0.2,
                strip=function(strip.levels,...) { 
                    strip.default(strip.levels=T, ...) },
                default.scales=list(arrows=F)
                )
        )
dev.off()

pdf("summary_sims_freq_hawk.pdf",width=10,height=7)
print(cloud((.5 * f_1 + (1-f_0-f_1)) + vc ~ c * v | d,
                data=the.data,
                #                data=the.data,
                pch=21,
                cex=0.2,
                strip=function(strip.levels,...) { 
                    strip.default(strip.levels=T, ...) },
                default.scales=list(arrows=F)
                )
        )
dev.off()

pdf("summary_sims_f0.pdf",width=10,height=7)
print(cloud(f_0 + f_1 + (1-f_0-f_1) ~ c * v | d,
                data=the.data,
                #                data=the.data,
                pch=21,
                cex=0.2,
                auto.key=T,
                strip=function(strip.levels,...) { 
                    strip.default(strip.levels=T, ...) },
                default.scales=list(arrows=F)
                )
        )
dev.off()

pdf("summary_sims_phh_min_pdh.pdf",width=10,height=7)
print(levelplot((meanphh - meanpdh) ~ c * v | d,
                data=the.data,
                #                data=the.data,
                #  pch=21,
                #cex=0.2,
                ,strip=function(strip.levels,...) { 
                    strip.default(strip.levels=T, ...) }
                ,default.scales=list(arrows=F)
                ,col.regions=matlab.like
                )
        )
dev.off()
