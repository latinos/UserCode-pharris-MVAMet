#ifndef PH_UTILS_LINKDEF_H
#define PH_UTILS_LINKDEF_H

#include "pharris/MVAMet/interface/GBRTree.h"
#include "pharris/MVAMet/interface/GBRForest.h"

#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class GBRTree+; 
#pragma link C++ class GBRForest+; 


#endif
