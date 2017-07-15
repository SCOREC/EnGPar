add_library(engpar_ftn engpar.f90)
set(HEADERS)
engpar_export_lib(engpar_ftn "${HEADERS}")
