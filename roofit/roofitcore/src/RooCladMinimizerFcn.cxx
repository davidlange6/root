/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   VC, Vince Croft,     DIANA / NYU,        vincent.croft@cern.ch          *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// RooCladMinimizerFcn is am interface class to the ROOT::Math function
// for minimization. See RooGradMinimizer.cxx for more information.
//

#include "RooCladMinimizerFcn.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"

#include "Riostream.h"
#include "TIterator.h"
#include "TClass.h"
#include "Fit/Fitter.h"
#include "Math/Minimizer.h"

#include <algorithm> // std::equal
#include <iostream>

RooCladMinimizerFcn::RooCladMinimizerFcn(RooAbsReal *funct, RooMinimizer *context, bool verbose)
   : RooAbsMinimizerFcn(RooArgList( * std::unique_ptr<RooArgSet>(funct->getParameters(RooArgSet{})) ), context, verbose),
     _grad(getNDim()), _grad_params(getNDim()), _funct(funct),
     has_been_calculated(getNDim())
{
   // TODO: added "parameters" after rewrite in april 2020, check if correct
   auto parameters = _context->fitter()->Config().ParamsSettings();
   synchronizeParameterSettings(parameters, kTRUE, verbose);
   synchronizeGradientParameterSettings(parameters);
   setStrategy(ROOT::Math::MinimizerOptions::DefaultStrategy());
}

RooCladMinimizerFcn::RooCladMinimizerFcn(const RooCladMinimizerFcn &other)
   : RooAbsMinimizerFcn(other), _grad(other._grad), _grad_params(other._grad_params), _funct(other._funct),
     has_been_calculated(other.has_been_calculated), none_have_been_calculated(other.none_have_been_calculated)
{
}

////////////////////////////////////////////////////////////////////////////////

double RooCladMinimizerFcn::DoEval(const double *x) const
{
      Bool_t parameters_changed = kFALSE;

   // Set the parameter values for this iteration
   for (unsigned index = 0; index < NDim(); index++) {
      // also check whether the function was already evaluated for this set of parameters
      parameters_changed |= SetPdfParamVal(index, x[index]);
   }

   // Calculate the function for these parameters
   RooAbsReal::setHideOffset(kFALSE);
   double fvalue = _funct->getVal();
   RooAbsReal::setHideOffset(kTRUE);

   if (!parameters_changed) {
      return fvalue;
   }

   if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {

      if (_printEvalErrors >= 0) {

         if (_doEvalErrorWall) {
            oocoutW(static_cast<RooAbsArg *>(nullptr), Eval)
               << "RooCladMinimizerFcn: Minimized function has error status." << std::endl
               << "Returning maximum FCN so far (" << _maxFCN
               << ") to force MIGRAD to back out of this region. Error log follows" << std::endl;
         } else {
            oocoutW(static_cast<RooAbsArg *>(nullptr), Eval)
               << "RooCladMinimizerFcn: Minimized function has error status but is ignored" << std::endl;
         }

         TIterator *iter = _floatParamList->createIterator();
         RooRealVar *var;
         Bool_t first(kTRUE);
         ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << "Parameter values: ";
         while ((var = (RooRealVar *)iter->Next())) {
            if (first) {
               first = kFALSE;
            } else
               ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << ", ";
            ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << var->GetName() << "=" << var->getVal();
         }
         delete iter;
         ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << std::endl;

         RooAbsReal::printEvalErrors(ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval), _printEvalErrors);
         ooccoutW(static_cast<RooAbsArg *>(nullptr), Eval) << std::endl;
      }

      if (_doEvalErrorWall) {
         fvalue = _maxFCN + 1;
      }

      RooAbsReal::clearEvalErrorLog();
      _numBadNLL++;
   } else if (fvalue > _maxFCN) {
      _maxFCN = fvalue;
   }

   // Optional logging
   if (_verbose) {
      std::cout << "\nprevFCN" << (_funct->isOffsetting() ? "-offset" : "") << " = " << std::setprecision(10) << fvalue
                << std::setprecision(4) << "  ";
      std::cout.flush();
   }

   _evalCounter++;
   return fvalue;
}

void RooCladMinimizerFcn::resetHasBeenCalculatedFlags() const
{
   for (auto it = has_been_calculated.begin(); it != has_been_calculated.end(); ++it) {
      *it = false;
   }
   none_have_been_calculated = true;
}

bool RooCladMinimizerFcn::syncParameter(double x, std::size_t ix) const
{
   bool parameter_has_changed = (_grad_params[ix] != x);

   if (parameter_has_changed) {
      _grad_params[ix] = x;
      // Set the parameter values for this iteration
      // TODO: this is already done in DoEval as well; find efficient way to do only once
      SetPdfParamVal(ix, x);

      if (!none_have_been_calculated) {
         resetHasBeenCalculatedFlags();
      }
   }

   return parameter_has_changed;
}

bool RooCladMinimizerFcn::syncParameters(const double *x) const
{
   bool has_been_synced = false;

   for (std::size_t ix = 0; ix < NDim(); ++ix) {
      bool parameter_has_changed = (_grad_params[ix] != x[ix]);

      if (parameter_has_changed) {
         _grad_params[ix] = x[ix];
         // Set the parameter values for this iteration
         // TODO: this is already done in DoEval as well; find efficient way to do only once
         SetPdfParamVal(ix, x[ix]);
      }

      has_been_synced |= parameter_has_changed;
   }

   if (has_been_synced) {
      resetHasBeenCalculatedFlags();
   }

   return has_been_synced;
}

void RooCladMinimizerFcn::runDerivator(unsigned int i_component) const
{
   // check whether the derivative was already calculated for this set of parameters
   if (!has_been_calculated[i_component]) {
      // Calculate the derivative etc for these parameters
      std::vector<double> gradResult(getNDim());
      _funct->evaluateGradient(gradResult.data());
      for(std::size_t iComp = 0; iComp < getNDim(); ++iComp) {
         _grad[iComp] = {gradResult[iComp], 0.0, 0.0};
         has_been_calculated[iComp] = true;
      }
      none_have_been_calculated = false;
   }
}

double RooCladMinimizerFcn::DoDerivative(const double *x, unsigned int i_component) const
{
   syncParameters(x);
   runDerivator(i_component);
   return _grad[i_component].derivative;
}

double RooCladMinimizerFcn::DoDerivativeWithPrevResult(const double *x, unsigned int i_component,
                                                       double * previous_grad, double * /*previous_g2*/,
                                                       double * /*previous_gstep*/) const
{
   syncParameters(x);
   _grad[i_component] = {previous_grad[i_component], 0.0, 0.0};
   runDerivator(i_component);
   previous_grad[i_component] = _grad[i_component].derivative;
   return _grad[i_component].derivative;
}

////////////////////////////////////////////////////////////////////////////////

void RooCladMinimizerFcn::setStrategy(int istrat)
{
   assert(istrat >= 0);
   ROOT::Minuit2::MnStrategy strategy(static_cast<unsigned int>(istrat));
}

Bool_t
RooCladMinimizerFcn::Synchronize(std::vector<ROOT::Fit::ParameterSettings> &parameters, Bool_t optConst, Bool_t verbose)
{
   Bool_t returnee = synchronizeParameterSettings(parameters, optConst, verbose);
   synchronizeGradientParameterSettings(parameters);
   setStrategy(_context->fitter()->Config().MinimizerOptions().Strategy());
   return returnee;
}
