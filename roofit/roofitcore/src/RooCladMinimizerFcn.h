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

#ifndef ROO_CLAD_MINIMIZER_FCN
#define ROO_CLAD_MINIMIZER_FCN

#include "RooAbsMinimizerFcn.h"

#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnMatrix.h" // MnAlgebraicVector
#include "Math/IFunction.h"   // IMultiGradFunction
#include "Fit/ParameterSettings.h"
#include "Minuit2/NumericalDerivator.h"

#include <vector>

class RooCladMinimizerFcn : public ROOT::Math::IMultiGradFunction, public RooAbsMinimizerFcn {
public:
   RooCladMinimizerFcn(RooAbsReal *funct, RooMinimizer *context, bool verbose = false);
   RooCladMinimizerFcn(const RooCladMinimizerFcn &other);
   ROOT::Math::IMultiGradFunction *Clone() const override { return new RooCladMinimizerFcn(*this); }

   ROOT::Minuit2::MnStrategy getStrategy() const;
   double getErrorDef() const;
   void setStrategy(int istrat);

   Bool_t Synchronize(std::vector<ROOT::Fit::ParameterSettings> &parameter_settings, Bool_t optConst,
                      Bool_t verbose = kFALSE) override;

   void synchronizeGradientParameterSettings(std::vector<ROOT::Fit::ParameterSettings> & /*parameter_settings*/) const {}

   inline bool returnsInMinuit2ParameterSpace() const override { return false; }
   inline unsigned int NDim() const override { return getNDim(); }
   inline std::string getFunctionName() const override { return _funct->GetName(); }
   inline std::string getFunctionTitle() const override { return _funct->GetTitle(); }
   inline void setOffsetting(Bool_t flag) override { _funct->enableOffsetting(flag); }

   bool fit(ROOT::Fit::Fitter& fitter) const override { return fitter.FitFCN(*this); };
   ROOT::Math::IMultiGenFunction* getMultiGenFcn() override { return this; };

private:
   void runDerivator(unsigned int i_component) const;

   void resetHasBeenCalculatedFlags() const;
   bool syncParameter(double x, std::size_t ix) const;
   bool syncParameters(const double *x) const;

   inline void setOptimizeConstOnFunction(RooAbsArg::ConstOpCode opcode, Bool_t doAlsoTrackingOpt) override
   {
      _funct->constOptimizeTestStatistic(opcode, doAlsoTrackingOpt);
   }

public:
   enum class GradientCalculatorMode { ExactlyMinuit2, AlmostMinuit2 };

private:
   // IMultiGradFunction overrides
   double DoEval(const double *x) const override;
   double DoDerivative(const double *x, unsigned int icoord) const override;
   double DoDerivativeWithPrevResult(const double *x, unsigned int i_component, double *previous_grad,
                                     double *previous_g2, double *previous_gstep) const override;

   // members
   // mutable because ROOT::Math::IMultiGradFunction::DoDerivative is const
protected:
   mutable std::vector<ROOT::Minuit2::DerivatorElement> _grad;
   mutable std::vector<double> _grad_params;

private:
   RooAbsReal *_funct;
   mutable std::vector<bool> has_been_calculated;
   mutable bool none_have_been_calculated = false;
};
#endif
