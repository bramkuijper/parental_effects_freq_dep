library("lattice")
library("colorRamps")

the.data <- read.table("summary_p_only.csv",sep=";",header=T)
the.data <- read.table("summary_p_only_new.csv",sep=";",header=T)

if (!("v" %in% names(the.data)))
{
    the.data[,"v"] <- round(1.0 - the.data[,"mh_1"],3)
    the.data[,"c"] <- round(2.0 * the.data[,"mh_2"] + the.data[,"v"] - 2,3)
}


the.data[,"vc"] <- the.data$v / the.data$c

the.data[,"vc"] <- ifelse(the.data$vc < 0, 0, the.data$vc)
the.data[,"vc"] <- ifelse(the.data$vc > 1, 1, the.data$vc)



pdf("overview_p_new.pdf")
print(wireframe(p ~ c * v | d,
                    data=the.data,
                    lwd=0.1,
                    zlim=c(0,1),
                    xlim=c(0,1),
                    ylim=c(0,1),
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()
pdf("overview_fi_new.pdf")
print(cloud(f_0 + f_1 + (1-f_0-f_1) ~ c * v | d,
                    data=the.data,
                    pch=21,
                    cex=0.1,
                    zlim=c(0,1),
                    xlim=c(0,1),
                    ylim=c(0,1),
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()


pdf("overview_vi_new.pdf")
print(cloud(vh_1 + vd_1 + vd_0 + vh_2 ~ c * v | d,
                    data=the.data,
                    pch=21,
                    cex=0.1,
                    xlim=c(0,1),
                    ylim=c(0,1),
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()


pdf("overview_vh1_vd1_new.pdf")
print(cloud(vh_1 + vd_1 ~ c * v | d,
                    data=the.data,
                    pch=21,
                    cex=0.1,
                    xlim=c(0,1),
                    ylim=c(0,1),
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()

pdf("overview_vh2_vd0_new.pdf")
print(cloud(vh_2 + vd_0 ~ c * v | d,
                    data=the.data,
                    pch=21,
                    cex=0.1,
                    xlim=c(0,1),
                    ylim=c(0,1),
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()

pdf("overview_f0_new.pdf")
print(cloud(f_0 ~ c * v | d,
                    data=the.data[the.data$f_0 > 0.01,],
                    pch=21,
                    cex=0.1,
                    zlim=c(0,1),
                    xlim=c(0,1),
                    ylim=c(0,1),
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()


pdf("overview_f0_new.pdf")
print(cloud(f_0 ~ c * v | d,
                    data=the.data[the.data$f_0 > 0.01,],
                    pch=21,
                    cex=0.1,
                    zlim=c(0,1),
                    xlim=c(0,1),
                    ylim=c(0,1),
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()
