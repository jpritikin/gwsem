set print thread-events off
dir src
set breakpoint pending on
# Rf_errorcall is used for some out of memory conditions
b Rf_errorcall
catch throw
# Hints:
# Use Rf_PrintValue(SEXP) to pretty print R data from gdb
