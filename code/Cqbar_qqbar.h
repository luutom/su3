// -*- C++ -*-
// $Id$
/*! \file
 *  \brief Meson- definitions
 */

#include "chromabase.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <qdp.h>
#include "mom_project.h"
#include "readwrite.h"

using namespace QDP;
using namespace Chroma ;

#ifndef __cqbar_qqbar_system_h__
#define __cqbar_qqbar_system_h__

class Cqbar_qqbar_system
{
public:
  
  // Default constructor -- forward declaration
  Cqbar_qqbar_system() {
  };  

  // blocks for the exchange terms
  multi1d<Propagator> Mblock_qq;
  multi1d<Propagator> Mblock_Qq;

  // initialization routine
  int initialize(const multi1d<int>& src,const multi1d<int>& mom, const multi1d<int>& latsize,
                               std::string directory, std::string base);

  int calcDirect(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index);
  int calcExchange(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index);
  int calcExchange3bar(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index);

  int calcDirect_P0(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index);
  int calcExchange_P0(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index);
  
  int calc_qq(const LatticePropagator& propQ, const int Gamma_index);
  int calc_Qq(const LatticePropagator& propQ, const LatticePropagator& propQbar, const int Gamma_index);

  int calc_15_rep(void);
  int calc_6_rep(void);
  int calc_3bar_rep(void);
  
  int reset_momentum(const multi1d<int>& src,const multi1d<int>& mom);

  // meson mass routine
  int write_qq_and_Qq(const int g1, const int g2, const std::string p);
  int write_15_and_6(const int g1, const int g2, const std::string p);
  int write_15(const int g1, const int g2, const std::string p);
  int write_6(const int g1, const int g2, const std::string p);
  int write_3bar(const int g1, const int g2, const std::string p);
  int write_direct(const int g1, const int g2, const std::string p);
  // consistency test routines 
  bool doChecks;
  
  // free memory
  int free();

  // meson correlator
  multi1d<DComplex> direct;
  multi1d<DComplex> exchange;
  multi1d<DComplex> corr_Qq;
  multi1d<DComplex> corr_qq;
  multi1d<DComplex> rep_15;
  multi1d<DComplex> rep_6;
  multi1d<DComplex> rep_3bar;
  
private:
    
  std::string dirname;
  std::string basename;

  multi1d<int> lat_size;

  // zero momentum container
  mom_project P1xyz;
  mom_project P2xyz;

  // read/write container
  readwrite rw;

  // where the meson source is located
  multi1d<int> srce_pt;
  multi1d<int> p_mom;

};


#endif
