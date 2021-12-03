#include "chroma.h"
#include "permutations.h"
#include <iostream>
#include <string>
#include "mom_project.h"
#include "readwrite.h"

#include "Cqbar_qqbar.h"

using namespace QDP;

int Cqbar_qqbar_system::initialize(const multi1d<int>& src,const multi1d<int>& mom, const multi1d<int>& latsize,
                               std::string directory, std::string base)
{
   srce_pt.resize(4);
   lat_size.resize(4);
   p_mom.resize(4);  
    
  /* Set the meson source position */
  srce_pt[0]=src[0];
  srce_pt[1]=src[1];
  srce_pt[2]=src[2];
  srce_pt[3]=src[3];

  /* Set lattice size constants */
  lat_size[0]=latsize[0];
  lat_size[1]=latsize[1];
  lat_size[2]=latsize[2];
  lat_size[3]=latsize[3];
  
  // momentum at the sink
  p_mom[0]=mom[0];
  p_mom[1]=mom[1];
  p_mom[2]=mom[2];
  p_mom[3]=mom[3];
  P1xyz.initialize(p_mom,srce_pt);

  p_mom[0]=-mom[0];
  p_mom[1]=-mom[1];
  p_mom[2]=-mom[2];
  P2xyz.initialize(p_mom,srce_pt);

  
  dirname = directory;
  basename = base;
    


  return 0;

}

int Cqbar_qqbar_system::reset_momentum(const multi1d<int>& src,const multi1d<int>& mom)
{

  srce_pt.resize(4);
  p_mom.resize(4);  

  srce_pt[0]=src[0];
  srce_pt[1]=src[1];
  srce_pt[2]=src[2];
  srce_pt[3]=src[3];
  
  p_mom[0]=mom[0];
  p_mom[1]=mom[1];
  p_mom[2]=mom[2];
  p_mom[3]=mom[3];
  P1xyz.initialize(p_mom,srce_pt);

  p_mom[0]=-mom[0];
  p_mom[1]=-mom[1];
  p_mom[2]=-mom[2];
  P2xyz.initialize(p_mom,srce_pt);
  
  return 0;
}

int Cqbar_qqbar_system::free()
{
  srce_pt.resize(0);
  lat_size.resize(0);
  p_mom.resize(0);
    
  /* Resize the udd blocks accordingly */
  Mblock_qq.resize(0);
  Mblock_Qq.resize(0);

  /* Set size of correlators */
  corr_qq.resize( 0 );
  corr_Qq.resize( 0 );
  direct.resize( 0 );
  exchange.resize( 0 );
  rep_15.resize( 0 );
  rep_6.resize( 0 );
  rep_3bar.resize( 0 );
  
  return 0;
  
}

int Cqbar_qqbar_system::calc_qq(const LatticePropagator& propq, const int Gamma_index)
{
  LatticeComplex temp; 
  SpinMatrixD g5Gamma = 1.0;
  SpinMatrixD Gammag5 = 1.0;
  /* Set size of correlator */
  corr_qq.resize( lat_size[3] );

  /*
    Calculates correlator for light-light meson system
   */

  // now calculate g5 * Gamma
  g5Gamma = g5Gamma * Gamma(15);  // g_5
  Gammag5 = Gamma(Gamma_index) * g5Gamma; // Gamma * g_5
  g5Gamma = g5Gamma * Gamma(Gamma_index);  // g_5 * Gamma

  if( doChecks ) {
  };
  
  temp = trace(Gammag5 * adj(propq) * g5Gamma * propq); // Note! I remove overall minus sign (same as hadspec)
  corr_qq = P1xyz.project_momentum( temp );


  return 0;

}

int Cqbar_qqbar_system::calc_Qq(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma_index)
{
  LatticeComplex temp; 
  SpinMatrixD g5Gamma = 1.0;
  SpinMatrixD Gammag5 = 1.0;
  /* Set size of correlator */
  corr_Qq.resize( lat_size[3] );

  /*
    Calculates correlator for heavy-light meson system
    Here propq is the propagator for the light quark and propQ is the propagator for the heavy quark
   */
  
  // now calculate g5 * Gamma
  g5Gamma = g5Gamma * Gamma(15);  // g_5
  Gammag5 = Gamma(Gamma_index) * g5Gamma; // Gamma * g_5
  g5Gamma = g5Gamma * Gamma(Gamma_index);  // g_5 * Gamma

  if( doChecks ) {
  };
  
  temp = trace(Gammag5 * adj(propq) * g5Gamma * propQ); // Note! I remove overall minus sign (asame as hadspec)
  corr_Qq = P2xyz.project_momentum( temp );


  return 0;

}

int Cqbar_qqbar_system::calcDirect(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index)
{
  //  /* Set size of correlator */
  direct.resize( lat_size[3] );

  //  /*
  //    Here propq is the propagator for the light quark and propQ is the propagator for the heavy quark, and Gamma = gamma_5 * gamma_i
  //   */

  if( doChecks ) {
  };

  calc_qq(propq, Gamma1_index);  // calculates tr(Gamma1 gamma5 qq^dag gamma5 Gamma1 qq) and momentum projects
  calc_Qq(propq, propQ, Gamma2_index); // calculates tr(Gamma2 gamma5 qq^dag gamma5 Gamma2 QQ) and momentum projects

  direct = corr_qq * corr_Qq;

  return 0;

}

int Cqbar_qqbar_system::calcDirect_P0(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index)
{
  LatticeComplex temp, temp1, temp2;
  SpinMatrixD g5Gamma1 = 1.0;
  SpinMatrixD Gamma1g5 = 1.0;
  SpinMatrixD g5Gamma2 = 1.0;
  SpinMatrixD Gamma2g5 = 1.0;
  //  /* Set size of correlator */
  direct.resize( lat_size[3] );

  //  /*
  //    Here propq is the propagator for the light quark and propQ is the propagator for the heavy quark, and Gamma = gamma_5 * gamma_i
  //   */

  // now calculate g5 * Gamma1
  g5Gamma1 = g5Gamma1 * Gamma(15);  // g_5
  Gamma1g5 = Gamma(Gamma1_index) * g5Gamma1; // Gamma1 * g_5
  g5Gamma1 = g5Gamma1 * Gamma(Gamma1_index);  // g_5 * Gamma1

   // now calculate g5 * Gamma2
  g5Gamma2 = g5Gamma2 * Gamma(15);  // g_5
  Gamma2g5 = Gamma(Gamma2_index) * g5Gamma2; // Gamma2 * g_5
  g5Gamma2 = g5Gamma2 * Gamma(Gamma2_index);  // g_5 * Gamma2
  
  if( doChecks ) {
  };

  // the follwoing calculates tr(Gamma1 gamma5 qq^dag gamma5 Gamma1 qq)*tr(Gamma2 gamma5 qq^dag gamma5 Gamma2 QQ) and then momentum projects
  temp1 = trace(Gamma1g5 * adj(propq) * g5Gamma1 * propq);
  temp2 = trace(Gamma2g5 * adj(propq) * g5Gamma2 * propQ);
  temp = temp1*temp2;

  direct = P1xyz.project_momentum( temp );

  return 0;

}

int Cqbar_qqbar_system::calcExchange(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index)
{
  LatticePropagator temp; 
  exchange.resize( lat_size[3] );
  multi1d < Propagator > Mblock;
  Mblock.resize( lat_size[3] );
  SpinMatrixD g5Gamma1 = 1.0;
  SpinMatrixD Gamma1g5 = 1.0;
  SpinMatrixD g5Gamma2 = 1.0;
  SpinMatrixD Gamma2g5 = 1.0;

  //  /* Resize the meson block accordingly */
  Mblock_qq.resize(lat_size[3]);
  Mblock_Qq.resize(lat_size[3]);

  // make gamma matrices;
  g5Gamma1 = g5Gamma1 * Gamma(15); // g5
  Gamma1g5 = Gamma(Gamma1_index) * g5Gamma1; // Gamma1 * g5
  g5Gamma1 = g5Gamma1 * Gamma(Gamma1_index); // g5 * Gamma1

  g5Gamma2 = g5Gamma2 * Gamma(15); // g5
  Gamma2g5 = Gamma(Gamma2_index) * g5Gamma2; // Gamma2 * g5
  g5Gamma2 = g5Gamma2 * Gamma(Gamma2_index); // g5 * Gamma2

  if( doChecks ) {
  };
  
  temp = Gamma2g5 * adj(propq) * g5Gamma1 * propq;  // G2 g5 propq^\dag g5 G1 propq
  Mblock_qq = P1xyz.project_momentum(temp);

  temp = Gamma1g5 * adj(propq) * g5Gamma2 * propQ;  // G1 g5 propq^\dag g5 G2 propQ
  Mblock_Qq = P2xyz.project_momentum(temp);

  // CHECK!!!!!!!
  for(int t = 0; t < lat_size[3]; t++) {  // Here I do the trace timeslice by timeslice
    Mblock[t] = Mblock_qq[t] * Mblock_Qq[t];
    exchange[t] = trace(Mblock[t]);
  }

  Mblock.resize( 0 );

  return 0;

}

int Cqbar_qqbar_system::calcExchange3bar(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index)
{
  LatticePropagator temp; 
  exchange.resize( lat_size[3] );
  multi1d < Propagator > Mblock;
  Mblock.resize( lat_size[3] );
  SpinMatrixD g5Gamma1 = 1.0;
  //  SpinMatrixD Gamma1g5 = 1.0;
  SpinMatrixD g5Gamma2 = 1.0;
  SpinMatrixD Gamma2g5 = 1.0;

  //  /* Resize the meson block accordingly */
  Mblock_qq.resize(lat_size[3]);
  Mblock_Qq.resize(lat_size[3]);

  // make gamma matrices;
  g5Gamma1 = g5Gamma1 * Gamma(15); // g5
  //  Gamma1g5 = Gamma(Gamma1_index) * g5Gamma1; // Gamma1 * g5
  g5Gamma1 = g5Gamma1 * Gamma(Gamma1_index); // g5 * Gamma1

  g5Gamma2 = g5Gamma2 * Gamma(15); // g5
  Gamma2g5 = Gamma(Gamma2_index) * g5Gamma2; // Gamma2 * g5
  //  g5Gamma2 = g5Gamma2 * Gamma(Gamma2_index); // g5 * Gamma2

  if( doChecks ) {
  };
  
  temp = Gamma2g5 * adj(propq) * g5Gamma1 * Gamma(15) * adj(propq) * Gamma(15);  // G2 g5 propq^\dag g5 G1 g5 propq^\dag g5
  Mblock_qq = P1xyz.project_momentum(temp);

  temp = Gamma(Gamma1_index) * propq * Gamma(Gamma2_index) * propQ;  // G1 propq G2 propQ
  Mblock_Qq = P2xyz.project_momentum(temp);

  // CHECK!!!!!!!
  for(int t = 0; t < lat_size[3]; t++) {  // Here I do the trace timeslice by timeslice
    Mblock[t] = Mblock_qq[t] * Mblock_Qq[t];
    exchange[t] = trace(Mblock[t]);
  }

  Mblock.resize( 0 );

  return 0;

}

int Cqbar_qqbar_system::calcExchange_P0(const LatticePropagator& propq, const LatticePropagator& propQ, const int Gamma1_index, const int Gamma2_index)
{
  LatticePropagator temp1, temp2;
  LatticeComplex temp;

  exchange.resize( lat_size[3] );
  SpinMatrixD g5Gamma1 = 1.0;
  SpinMatrixD Gamma1g5 = 1.0;
  SpinMatrixD g5Gamma2 = 1.0;
  SpinMatrixD Gamma2g5 = 1.0;

  // make gamma matrices;
  g5Gamma1 = g5Gamma1 * Gamma(15); // g5
  Gamma1g5 = Gamma(Gamma1_index) * g5Gamma1; // Gamma1 * g5
  g5Gamma1 = g5Gamma1 * Gamma(Gamma1_index); // g5 * Gamma1

  g5Gamma2 = g5Gamma2 * Gamma(15); // g5
  Gamma2g5 = Gamma(Gamma2_index) * g5Gamma2; // Gamma2 * g5
  g5Gamma2 = g5Gamma2 * Gamma(Gamma2_index); // g5 * Gamma2

  if( doChecks ) {
  };

  temp1 = Gamma2g5 * adj(propq) * g5Gamma1 * propq;  // G2 g5 propq^\dag g5 G1 propq
  temp2 = Gamma1g5 * adj(propq) * g5Gamma2 * propQ;  // G1 g5 propq^\dag g5 G2 propQ
  temp = trace( temp1*temp2 );

  exchange = P1xyz.project_momentum( temp );
  
  return 0;

}

int Cqbar_qqbar_system:: calc_15_rep()
{
  rep_15.resize(lat_size[3]);
  
  rep_15 = direct - exchange;
  
  return 0;
}

int Cqbar_qqbar_system::calc_6_rep()
{
  rep_6.resize(lat_size[3]);
  
  rep_6 = direct + exchange;
  
  return 0;
}

int Cqbar_qqbar_system::calc_3bar_rep()
{
  rep_3bar.resize(lat_size[3]);

  rep_3bar = exchange; // + direct

  return 0;
}

int Cqbar_qqbar_system::write_qq_and_Qq(const int gamma1, const int gamma2, const std::string p)
{
  string g1, g2;

  g1 = to_string(gamma1); // for light-light system
  g2 = to_string(gamma2); // for heavy-light system

  rw.writeCorrelator(corr_Qq,dirname+"Q_g"+g2+"_qbar"+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(corr_Qq,dirname+"Q_g"+g2+"_qbar"+p+basename,1, lat_size[3], srce_pt[3]);

  rw.writeCorrelator(corr_qq,dirname+"q_g"+g1+"_qbar"+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(corr_qq,dirname+"q_g"+g1+"_qbar"+p+basename,1, lat_size[3], srce_pt[3]);

  return 0;
}

int Cqbar_qqbar_system::write_15(const int gamma1, const int gamma2, const std::string p)
{
  string g1, g2;

  g1 = to_string(gamma1); // for light-light system
  g2 = to_string(gamma2); // for heavy-light system

  rw.writeCorrelator(rep_15,dirname+"15rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(rep_15,dirname+"15rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt[3]);

  return 0;
}

int Cqbar_qqbar_system::write_6(const int gamma1, const int gamma2, const std::string p)
{
  string g1, g2;

  g1 = to_string(gamma1); // for light-light system
  g2 = to_string(gamma2); // for heavy-light system

  rw.writeCorrelator(rep_6,dirname+"6rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(rep_6,dirname+"6rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt[3]);
  
  return 0;
}

int Cqbar_qqbar_system::write_direct(const int gamma1, const int gamma2, const std::string p)
{
  string g1, g2;

  g1 = to_string(gamma1); // for light-light system
  g2 = to_string(gamma2); // for heavy-light system

  rw.writeCorrelator(direct,dirname+"direct_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(direct,dirname+"direct_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt[3]);
  return 0;
}

int Cqbar_qqbar_system::write_3bar(const int gamma1, const int gamma2, const std::string p)
{
  string g1, g2;

  g1 = to_string(gamma1); // for light-light system
  g2 = to_string(gamma2); // for heavy-light system

  rw.writeCorrelator(rep_3bar,dirname+"3brep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(rep_3bar,dirname+"3brep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt[3]);
  return 0;
}

int Cqbar_qqbar_system::write_15_and_6(const int gamma1, const int gamma2, const std::string p)
{
  string g1, g2;

  g1 = to_string(gamma1); // for light-light system
  g2 = to_string(gamma2); // for heavy-light system

  rw.writeCorrelator(rep_6,dirname+"6rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(rep_6,dirname+"6rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt[3]);
  
  rw.writeCorrelator(rep_15,dirname+"15rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(rep_15,dirname+"15rep_g"+g1+"g"+g2+p+basename,1, lat_size[3], srce_pt[3]);

  return 0;
}
