library(pgenlibr)

dir <- "/home/joshua/vcu/gwsem"
pvar <- NewPvar(paste(dir, "inst/extdata/example.pvar", sep='/'))
d1 <- NewPgen(paste(dir, "inst/extdata/example.pgen", sep='/'),
              pvar)
ct <- GetRawSampleCt(d1)
buf <- matrix(0, 2, ct)
ReadAlleles(d1, buf, 1)

GetVariantId(pvar, 3)
GetAlleleCode(pvar, 2, 2)
