library("lattice")
library("colorRamps")

the.data <- read.table("summary_phh_pdh.csv",sep=";",header=T)

print(levelplot((phh - pdh) ~ c * v | d,
                    strip=function(strip.levels,...) 
                            { strip.default(strip.levels=T,...) },
                    data=the.data[the.data$c >= the.data$v,],
                    col.regions=matlab.like,
                    zlim=c(-1,0),
                    default.scales=list(arrows=F)
                    )
            )
