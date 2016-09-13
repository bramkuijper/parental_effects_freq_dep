library("lattice")
library("colorRamps")

the.data <- read.table("all_sims.csv",sep=";",header=T)

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
                data=the.data[the.data$v < the.data$c,],
                #                data=the.data,
                pch=21,
                cex=0.2,
                strip=function(strip.levels,...) { 
                    strip.default(strip.levels=T, ...) },
                default.scales=list(arrows=F)
                )
        )
dev.off()

pdf("summary_sims_phh_min_pdh.pdf",width=10,height=7)
print(wireframe((meanphh - meanpdh) ~ c * v | d,
                data=the.data[the.data$v < the.data$c,],
                #                data=the.data,
                pch=21,
                cex=0.2,
                strip=function(strip.levels,...) { 
                    strip.default(strip.levels=T, ...) },
                default.scales=list(arrows=F)
                )
        )
dev.off()
