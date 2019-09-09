#include "openmx.h"
#include "LoadDataAPI.h"

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

static const R_CallMethodDef CallEntries[] = {
  {NULL, NULL, 0}
};

void setup2(AddLoadDataProviderType);

extern "C" {
  void attribute_visible R_init_gwsem(DllInfo *dll) {
    // next line is necessary to avoid a NOTE from R CMD check
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, Rboolean::TRUE); // necessary for .onLoad() to work
    AddLoadDataProviderType AddLoadDataProvider =
      (AddLoadDataProviderType) R_GetCCallable("OpenMx", "AddLoadDataProvider");
    setup2(AddLoadDataProvider);
  }
}
