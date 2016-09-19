library("lattice")

if (!exists("the.data"))
{
the.data <- read.table("summary_allp.csv",sep=";",header=T)
the.data[,"f_2"] <- 1.0 - the.data[,"f_0"] - the.data[,"f_1"]

# make plots of migrant rv
the.data[,"rv.hawk.m"] <- (
            the.data[,"f_0"] * 2 * (1 - the.data[,"v"]/2) * the.data[,"vh_1"]
            +
            # dove dies on a hd patch
            the.data[,"f_1"] * 1.0 * the.data[,"vh_2"]
            +
            # hawk dies on a hd patch
            the.data[,"f_1"] * (1-the.data[,"v"]) * the.data[,"vh_1"]
            +
            # hawk dies on a dd patch
            the.data[,"f_2"] * (1 - (the.data[,"v"]/2 - the.data[,"c"]/2)) * the.data[,"vh_2"])

the.data[,"rv.dove.m"] <- (
            the.data[,"f_0"] * 2 * (1 - the.data[,"v"]/2) * the.data[,"vd_0"]
            +
            # dove dies on a hd patch
            the.data[,"f_1"] * 1.0 * the.data[,"vd_1"]
            +
            # hawk dies on a hd patch
            the.data[,"f_1"] * (1-the.data[,"v"]) * the.data[,"vd_0"]
            +
            # hawk dies on a dd patch
            the.data[,"f_2"] * (1 - (the.data[,"v"]/2 - the.data[,"c"]/2)) * the.data[,"vd_1"])
the.data[,"diff.hd"] <- the.data[,"vd_1"] + .5 * the.data[,"vh_1"] - 1.5 * the.data[,"vh_2"]
}

pdf("overview_freq_hawk.pdf",width=15,height=10)
print(wireframe((f_2 + f_1) ~ c * v | d,data=the.data,
        zlim=c(0,1),
        zlab=list(rot=100,label="Frequency hawk"),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()

pdf("overview_vd0.pdf",width=15,height=10)
print(wireframe(vd_0 ~ c * v | d,data=the.data,
        zlab=list(rot=100,label="RV Dove on DD patch"),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
pdf("overview_vh2_vd1.pdf",width=15,height=10)
print(wireframe(vh_2 - vd_1 ~ c * v | d,data=the.data,
        zlim=c(0,0.6),
        zlab=list(rot=100,label="RV Dove on DD patch"),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
pdf("overview_vd1.pdf",width=15,height=10)
print(wireframe(vd_1 ~ c * v | d,data=the.data,
        zlab=list(rot=100,label="RV Dove on HD patch"),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
pdf("overview_vh1.pdf",width=15,height=10)
print(wireframe(vh_1 ~ c * v | d,data=the.data,
        zlab=list(rot=100,label="RV Hawk on HD patch"),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
pdf("overview_vh2.pdf",width=15,height=10)
print(wireframe(vh_2 ~ c * v | d,data=the.data,
        zlab=list(rot=100,label="RV Hawk on HH patch"),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()

pdf("overview_pdh0.pdf",width=15,height=10)
print(wireframe(pdh0 ~ c * v | d,data=the.data,
        zlim=c(0,1),
        zlab=list(rot=100,label="% H by dove on DD patch"),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
pdf("overview_pdh1.pdf",width=15,height=10)

print(wireframe(pdh1 ~ c * v | d,data=the.data,
        zlim=c(0,1),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
pdf("overview_phh1.pdf",width=15,height=10)

print(wireframe(phh1 ~ c * v | d,data=the.data,
        zlim=c(0,1),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
pdf("overview_phh2.pdf",width=15,height=10)

print(wireframe(phh2 ~ c * v | d,data=the.data,
        zlim=c(0,1),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()

pdf("overview_rv_hawk_mig.pdf",width=15,height=10)

print(wireframe(rv.hawk.m ~ c * v | d,data=the.data,
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()

pdf("overview_rv_dove_mig.pdf",width=15,height=10)

print(wireframe(rv.dove.m ~ c * v | d,data=the.data,
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()

pdf("overview_rvdiff_dove_hawk_mig.pdf",width=15,height=10)

print(wireframe(rv.hawk.m - rv.dove.m ~ c * v | d,data=the.data,
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()

pdf("overview_diffhd",width=15,height=10)

print(wireframe(diff.hd ~ c * v | d,data=the.data,
        zlim=c(0,1),
        strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
        default.scales=list(arrows=F)
        ))
dev.off()
