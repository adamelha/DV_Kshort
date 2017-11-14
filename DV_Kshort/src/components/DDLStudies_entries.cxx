// Please place this file in DDLStudies/src/components/

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "../Kshort_DDL.h"// This was added  for Kshort_DDL

#ifndef XAOD_ANALYSIS
DECLARE_NAMESPACE_ALGORITHM_FACTORY(DDL, Kshort_DDL) // This was added  for Kshort_DDL
#endif

DECLARE_FACTORY_ENTRIES(DDLStudies)
{
  #ifndef XAOD_ANALYSIS
  DECLARE_NAMESPACE_ALGORITHM(DDL, Kshort_DDL);// This was added  for Kshort_DDL
  #endif
}
