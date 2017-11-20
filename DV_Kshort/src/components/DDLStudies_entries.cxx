#include "GaudiKernel/DeclareFactoryEntries.h"
#include "../Ks.h"

#ifndef XAOD_ANALYSIS
DECLARE_NAMESPACE_ALGORITHM_FACTORY(DDL, Ks) 
#endif

DECLARE_FACTORY_ENTRIES(DDLStudies)
{
  #ifndef XAOD_ANALYSIS
  DECLARE_NAMESPACE_ALGORITHM(DDL, Ks);
  #endif
}
