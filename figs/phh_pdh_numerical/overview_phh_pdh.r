library("lattice")
library("colorRamps")

the.data <- read.table("summary_phh_pdh_new.csv",sep=";",header=T)


if (!("v" %in% names(the.data)))
{
    the.data[,"v"] <- round(1.0 - the.data[,"mh_1"],3)
    the.data[,"c"] <- round(2.0 * the.data[,"mh_2"] + the.data[,"v"] - 2,3)
} else 
{
    the.data[,"err"] <- 0
}




pdf("overview_phh_new.pdf")
print(wireframe(phh ~ c * v | d * err,
                    data=the.data,
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
#                    data=the.data[the.data$c >= the.data$v,],
                    zlim=c(0,1),
                    default.scales=list(arrows=F)
                    )
            )
dev.off()

pdf("overview_pdh_new.pdf")
print(wireframe(pdh ~ c * v | d * err,
                    data=the.data,
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
#                    data=the.data[the.data$c >= the.data$v,],
                    zlim=c(0,1),

                    default.scales=list(arrows=F)
                    )
            )
dev.off()

pdf("overview_phh_pdh_diff.pdf")
print(levelplot((phh - pdh) ~ c * v | d * err,
                    data=the.data[the.data$v < the.data$c,],
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
                    col.regions=matlab.like
#                    data=the.data[the.data$c >= the.data$v,],
                            #                    zlim=c(0,1),
                            #                    default.scales=list(arrows=F)
                    )
            )
dev.off()


pdf("freq_hawk.pdf")
print(wireframe((.5 * f_1 + (1-f_0-f_1)) ~ c * v | d * err,
                    data=the.data,
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
#                    data=the.data[the.data$c >= the.data$v,],
                    zlim=c(0,1),
                    default.scales=list(arrows=F)
                    )
            )
dev.off()


pdf("rvals.pdf")
print(cloud(vh_1 + vh_2 + vd_0 + vd_1 ~ c * v | d * err,
                    data=the.data,
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
#                    data=the.data[the.data$c >= the.data$v,],
                    auto.key=T,
                    pch=".",
                    default.scales=list(arrows=F)
                    )
            )
dev.off()

pdf("rvals_vd0_vd1.pdf")
print(cloud(vd_0 + vd_1 ~ c * v | d * err,
                    data=the.data,
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
#                    data=the.data[the.data$c >= the.data$v,],
                    auto.key=T,
                    pch=".",
                    default.scales=list(arrows=F)
                    )
            )
dev.off()

pdf("rvals_vh1_vd0.pdf")
print(cloud(vh_1 + vd_0 ~ c * v | d * err,
                    data=the.data,
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
#                    data=the.data[the.data$c >= the.data$v,],
                    auto.key=T,
                    pch=".",
                    default.scales=list(arrows=F)
                    )
            )
dev.off()
