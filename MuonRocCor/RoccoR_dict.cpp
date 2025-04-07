// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RoccoR_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "RoccoR.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *RoccoR_Dictionary();
   static void RoccoR_TClassManip(TClass*);
   static void *new_RoccoR(void *p = nullptr);
   static void *newArray_RoccoR(Long_t size, void *p);
   static void delete_RoccoR(void *p);
   static void deleteArray_RoccoR(void *p);
   static void destruct_RoccoR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RoccoR*)
   {
      ::RoccoR *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RoccoR));
      static ::ROOT::TGenericClassInfo 
         instance("RoccoR", "RoccoR.h", 122,
                  typeid(::RoccoR), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RoccoR_Dictionary, isa_proxy, 0,
                  sizeof(::RoccoR) );
      instance.SetNew(&new_RoccoR);
      instance.SetNewArray(&newArray_RoccoR);
      instance.SetDelete(&delete_RoccoR);
      instance.SetDeleteArray(&deleteArray_RoccoR);
      instance.SetDestructor(&destruct_RoccoR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RoccoR*)
   {
      return GenerateInitInstanceLocal(static_cast<::RoccoR*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::RoccoR*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RoccoR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::RoccoR*>(nullptr))->GetClass();
      RoccoR_TClassManip(theClass);
   return theClass;
   }

   static void RoccoR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_RoccoR(void *p) {
      return  p ? new(p) ::RoccoR : new ::RoccoR;
   }
   static void *newArray_RoccoR(Long_t nElements, void *p) {
      return p ? new(p) ::RoccoR[nElements] : new ::RoccoR[nElements];
   }
   // Wrapper around operator delete
   static void delete_RoccoR(void *p) {
      delete (static_cast<::RoccoR*>(p));
   }
   static void deleteArray_RoccoR(void *p) {
      delete [] (static_cast<::RoccoR*>(p));
   }
   static void destruct_RoccoR(void *p) {
      typedef ::RoccoR current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::RoccoR

namespace {
  void TriggerDictionaryInitialization_RoccoR_dict_Impl() {
    static const char* headers[] = {
"RoccoR.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/thayil/boost/include",
"/cvmfs/cms.cern.ch/el8_amd64_gcc12/lcg/root/6.30.03-723f04ba093d0553281d42c7b0f6eee1/include/",
"/cms/thayil/pseudoaxions/attoAODs/CMSSW_14_0_0/src/PhysicsTools/NanoAODTools/MuonRocCor/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RoccoR_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$RoccoR.h")))  RoccoR;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RoccoR_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "RoccoR.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"RoccoR", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RoccoR_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RoccoR_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RoccoR_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RoccoR_dict() {
  TriggerDictionaryInitialization_RoccoR_dict_Impl();
}
