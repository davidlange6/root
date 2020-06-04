// @(#)root/thread:$Id$
// // Author: Xavier Valls Pla   08/05/20
//
/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RTaskArena                                                           //
//                                                                      //
// This file implements the method to initialize and retrieve ROOT's    //
// global task arena, together with a method to check for active        //
// CPU bandwith control, and a class to wrap the tbb task arena with    //
// the purpose of keeping tbb off the installed headers                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_RTaskArena
#define ROOT_RTaskArena

#include "RConfigure.h"
#include <memory>

// exclude in case ROOT does not have IMT support
#ifndef R__USE_IMT
// No need to error out for dictionaries.
# if !defined(__ROOTCLING__) && !defined(G__DICTIONARY)
#  error "Cannot use ROOT::Internal::RTaskArenaWrapper without defining R__USE_IMT."
# endif
#else

/// tbb::task_arena is an alias of tbb::interface7::task_arena, which doesn't allow
/// to forward declare tbb::task_arena without forward declaring tbb::interface7
namespace tbb{
namespace interface7{class task_arena;}
using task_arena = interface7::task_arena;
}

namespace ROOT {
namespace Internal {

////////////////////////////////////////////////////////////////////////////////
/// Wrapper for tbb::task_arena.
///
/// Necessary in order to keep tbb away from ROOT headers.
/// This class is thought out to be used as a singleton.
////////////////////////////////////////////////////////////////////////////////
class RTaskArenaWrapper {
public:
   RTaskArenaWrapper(unsigned maxConcurrency = 0);
   ~RTaskArenaWrapper(); // necessary to set size back to zero
   static unsigned TaskArenaSize();
   std::unique_ptr<tbb::task_arena> &Access();
private:
   std::unique_ptr<tbb::task_arena> fTBBArena;
   static unsigned fNWorkers;
};


////////////////////////////////////////////////////////////////////////////////
// Factory function returning a shared pointer to the instance of the global
// RTaskArenaWrapper. The task arena still needs to be initialized with a
// call to RTaskArena::Initialize(maxConcurrency)
//
// Allows for reinstantiation of the global RTaskArenaWrapper once all the
// references to the previous one are gone and the object destroyed.
////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<ROOT::Internal::RTaskArenaWrapper> GetGlobalTaskArena(unsigned maxConcurrency = 0);

} // namespace Internal
} // namespace ROOT

#endif   // R__USE_IMT
#endif   // ROOT_RTaskArena