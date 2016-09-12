library("lattice")
library("colorRamps")

the.data <- read.table("summary_p_only.csv",sep=";",header=T)

pdf("overview_p.pdf")
print(wireframe(p ~ c * v | d,
                    data=the.data,
                    lwd=0.1,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()
pdf("overview_fi.pdf")
print(cloud(f_0 + f_1 + (1-f_0-f_1) ~ c * v | d,
                    data=the.data,
                    pch=21,
                    cex=0.1,
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()

pdf("overview_f0.pdf")
print(cloud(f_0 ~ c * v | d,
                    data=the.data[the.data$f_0 > 0.01,],
                    pch=21,
                    cex=0.1,
                    auto.key=T,
                    strip=function(strip.levels, ...) { strip.default(strip.levels=T, ...) },
                    default.scales=list(arrows=F)
                    )
            )
dev.off()
