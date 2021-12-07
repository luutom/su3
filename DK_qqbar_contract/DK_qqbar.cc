// Contraction code for qqbar blocks generated from the QQQ_NUCNUC Chroma routine
// We write out all individual momentum options which can later be combined into the correct irreps.
// 2021-12-06 Andre Walker-Loud
//

#include <iostream>
#include <cstdio>
#include <libgen.h>
#include <math.h>
#include <iomanip>

#include "chroma.h"
#include "meas/hadron/qqbar_w.h"
#include "read_meson_block.h"
#include "DK_qqbar.h"

using namespace std;
using namespace QDP;
using namespace Chroma ;

const double PI = 2*acos(0.0);

//
//  Contraction routine
//  15: Tr[A]Tr[B] - Tr[A.B]
//   6: Tr[A]Tr[B] + Tr[A.B]
//
void DK_15(multi1d<DComplex>& corr,
    multi1d<DPropagator>& meson_block_p1,
    multi1d<DPropagator>& meson_block_p2)
    {
    int Nt = meson_block_p1.size() ;
    corr.resize(Nt) ;

    for(int t(0);t<Nt;t++){
        DPropagator Phi1 = meson_block_p1[t]*Gamma(Ns*Ns - 1) ; // Gamma15 = gamma_5
        DPropagator Phi2 = meson_block_p2[t]*Gamma(Ns*Ns - 1) ;
        corr[t] = trace(Phi1)*trace(Phi2) - trace(Phi1*Phi2) ;
        }
    }

void DK_6(multi1d<DComplex>& corr,
    multi1d<DPropagator>& meson_block_p1,
    multi1d<DPropagator>& meson_block_p2)
    {
    int Nt = meson_block_p1.size() ;
    corr.resize(Nt) ;

    for(int t(0);t<Nt;t++){
        DPropagator Phi1 = meson_block_p1[t]*Gamma(Ns*Ns - 1) ; // Gamma15 = gamma_5
        DPropagator Phi2 = meson_block_p2[t]*Gamma(Ns*Ns - 1) ;
        corr[t] = trace(Phi1)*trace(Phi2) + trace(Phi1*Phi2) ;
        }
    }

//
//  Main to do the DK in the 15 and 6 SU(3) flavor representations
//
int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);
  START_CODE();

  if (argc < 4){
    QDPIO::cerr << "Usage: "<<argv[0]<<" Nspace Ntime block1 [block2 ...] "<<endl;
    QDP_abort(1);
  }

  int Nspace = atoi(argv[1]);
  int Ntime  = atoi(argv[2]);
  string in_file(argv[3]);
  string out_file(basename(argv[3]));

  if(out_file.size() > 3)
    if(out_file.substr(out_file.size() - 3) == "dat")
      out_file.replace(out_file.end() - 3, out_file.end(), "lst");
    else
      out_file.append(".lst");
  else
    out_file.append(".lst");

  string D_file     = "D_"+out_file ;
  string kaon_file  = "kaon_"+out_file ;
  string DK_15_file = "DK_15_"+out_file ;
  string DK_6_file  = "DK_6_"+out_file ;

  multi1d <int> nrow(4) ;
  nrow[0] = nrow[1] = nrow[2] = Nspace; //assume Lx=Ly=Lz
  nrow[3] = Ntime;
  Layout::setLattSize(nrow);
  Layout::create();

  QDPIO::cout << " D(p1)-K(p2): contractions" << endl;

  XMLReader file_xml ;
  QDPFileReader from(file_xml,in_file,QDPIO_SERIAL) ;
  cout<<"OK\n";
  cout.flush();

  PropSourceConst_t src_head ;
  PropSinkSmear_t snk_head ;
  if(file_xml.count("/qqqNucNuc_w/Propagator_info/Propagator[1]/PropSource")!=0){
    read(file_xml,"/qqqNucNuc_w/Propagator_info/Propagator[1]/PropSource",src_head);
    read(file_xml,"/qqqNucNuc_w/Propagator_info/Propagator[1]/PropSink",snk_head);
  }
  else if(
      file_xml.count(
     "/qqqNucNuc_w/Propagator_info/Propagator/elem[1]/ForwardProp/PropSource"
     )!=0){
    read(file_xml,
     "/qqqNucNuc_w/Propagator_info/Propagator/elem[1]/ForwardProp/PropSource",
     src_head);
    read(file_xml,
     "/qqqNucNuc_w/Propagator_info/Propagator/elem[1]/ForwardProp/PropSink",
     snk_head);
  }
  else{
    QDPIO::cerr<<"OOPS! can't get the propagator source!\n";
    exit(1);
  }

  string snk;
  std::istringstream xml_s(snk_head.sink.xml);
  XMLReader snk_top(xml_s);
  XMLReader top(snk_top, snk_head.sink.path);
  read(top,"SinkType",snk);
  snk.replace(snk.end()-5,snk.end(),"");
  cout<<"Found a "<<snk<<" sink"<<endl;

  int j_decay = src_head.j_decay;

  // get number of momentum states
  int Nmom;
  read(file_xml,"/qqqNucNuc_w/MomNum",Nmom);

  //read all data
  read(file_xml,"/qqqNucNuc_w/ProgramInfo/Setgeom/latt_size",nrow) ;
  int Nt(nrow[j_decay]);
  int Nx(nrow[0]);
  cout<<"  NL="<<Nx<<" NT="<<Nt<<"\n"<<endl;
  //cout.flush();
  // Get Kaon block = "pion" from code
  ContainerQQbar QQbarK("pion",snk) ;
  cout<<"made ContainerQQbar QQBarPion"<<endl;
  QQbarK.setSize(Nmom,Nt,Nd-1);
  cout<<"set size of QQbar"<<endl;
  read_meson_block(from, QQbarK) ;
  cout<<"got kaon"<<endl;
  // Get D block = "kaon" from code
  ContainerQQbar QQbarD("kaon",snk) ;
  QQbarD.setSize(Nmom,Nt,Nd-1);
  read_meson_block(from, QQbarD) ;
  cout<<"got D-meson"<<endl;
  cout<<"Nmom = " << Nmom << endl;

  // close file reader
  close(from);

  // Now loop over momentum
  // We need two momentum loops since we are doing boosted systems
  multi1d <int> p(3);
  multi1d <int> q(3);
  multi1d <int> Ptot(3);
  int PtotSq=0;

  // loop over Kaon momentum
  for( int mom_p(0); mom_p < Nmom; mom_p++){
    p = QQbarK.Mom(mom_p);
    stringstream ps;
    ps <<"px"<<p[0]<<"py"<<p[1]<<"pz"<<p[2];
    stringstream psR;
    psR <<"qx"<<p[0]<<"qy"<<p[1]<<"qz"<<p[2];

    // loop over D-meson momentum
    for( int mom_q(0); mom_q < Nmom; mom_q++){
      q = QQbarD.Mom(mom_q);
      stringstream qs;
      qs <<"qx"<<q[0]<<"qy"<<q[1]<<"qz"<<q[2];
      stringstream qsR;
      qsR <<"px"<<q[0]<<"py"<<q[1]<<"pz"<<q[2];

      PtotSq = 0;
      for( int i(0); i < 3; i++){
    Ptot[i] = p[i] + q[i];
    PtotSq += Ptot[i] * Ptot[i];
      }
      stringstream pTs;
      pTs <<"Px"<<Ptot[0]<<"Py"<<Ptot[1]<<"Pz"<<Ptot[2];

      string ms     = pTs.str() + "_" + ps.str() + "_" + qs.str() + "_";
      string msR    = pTs.str() + "_" + qsR.str() + "_" + psR.str() + "_";

      string DK_15  = ms  + DK_15_file ;
      string DK_15R = msR + DK_15_file ;
      string DK_6   = ms  + DK_6_file ;
      string DK_6R  = ms  + DK_6_file ;

      if(Ptot[0] == 0 && Ptot[1] == 0 && Ptot[2] == 0){
    cout << "Ptot: "<< Ptot[0] <<" "<< Ptot[1] <<" "<< Ptot[2]<<": ";
    cout << "   p:"<< p[0]<<" "<< p[1] <<" "<< p[2]<<": ";
    cout << "   q:"<< q[0]<<" "<< q[1] <<" "<< q[2]<<": " << endl;
      }
      // DK_15
      multi1d<DComplex> C15(Nt) ;
      DoDK_15(QQbarK, QQbarD, p, q, Nx, C15) ; // pi(p) k(q)
      {
    ofstream jlist(DK_15.c_str());
    jlist <<"1 "<<Nt<<" 1 "<<Nx<<" 1 "<<endl ;
    for(int t(0); t<Nt;t++)
      {
        jlist <<t<<" "<< std::setprecision(16) <<real(C15[t])<<" "<<imag(C15[t])<<endl;
      }
    jlist.close();
      }
      // DK_6
      multi1d<DComplex> C6(Nt) ;
      DoDK_6(QQbarK, QQbarD, p, q, Nx, C6) ; // pi(p) k(q)
      {
    ofstream jlist(DK_6.c_str());
    jlist <<"1 "<<Nt<<" 1 "<<Nx<<" 1 "<<endl ;
    for(int t(0); t<Nt;t++)
      {
        jlist <<t<<" "<< std::setprecision(16) <<real(C6[t])<<" "<<imag(C6[t])<<endl;
      }
    jlist.close();
      }
    }// q-mom loop
  }// p-mom loop
  END_CODE();
  QDP_finalize();
  exit(0);
}
