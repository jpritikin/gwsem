# gdb:
# dir /opt/rbgen_v1.1.4/src/

library(rbgen)

ranges = data.frame(
  chromosome = "01",
  start = 1,
  end = 10000
)

dir <- "/home/joshua/vcu/gwsem"
d1 <- bgen.load(paste(dir, "inst/extdata/example.bgen", sep='/'), ranges)
colSums(d1$data[2,,])
d1$data[2,1:10,]
