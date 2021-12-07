// Contraction code for DK with qqbar blocks
// 2021-12-06: Andre Walker-Loud

#ifndef __DK_H__
#define __DK_H__

using namespace Chroma ;

//
//  Spin Contraction routine
//
void DK_15(multi1d<DComplex>& corr,
           multi1d<DPropagator>& meson_block_p1,
           multi1d<DPropagator>& meson_block_p2) ;
void DK_6(multi1d<DComplex>& corr,
          multi1d<DPropagator>& meson_block_p1,
          multi1d<DPropagator>& meson_block_p2) ;

//
//  make D-K correlators
//
void DoDK_15(ContainerQQbar& QQbarP, ContainerQQbar& QQbarQ,
    multi1d<int> p, multi1d<int> q, int Nx, multi1d<DComplex>& C)
{

  int Nt(QQbarP.Nt()) ;
  if(QQbarQ.Nt()!=Nt){
    QDPIO::cerr<<"OOPS missmatch in Nt!\n" ;
    QDP_abort(1) ;
  }
  multi1d<DPropagator> qqbar_p(Nt) ; //momentum p
  multi1d<DPropagator> qqbar_q(Nt) ; //momentum q

  qqbar_p = QQbarP[p] ;
  qqbar_q = QQbarQ[q] ;

  if( C.size() != Nt){
        C.resize(Nt);
  }

  DK_15(C, qqbar_p, qqbar_q);

}
void DoDK_6(ContainerQQbar& QQbarP, ContainerQQbar& QQbarQ,
    multi1d<int> p, multi1d<int> q, int Nx, multi1d<DComplex>& C)
{

  int Nt(QQbarP.Nt()) ;
  if(QQbarQ.Nt()!=Nt){
    QDPIO::cerr<<"OOPS missmatch in Nt!\n" ;
    QDP_abort(1) ;
  }
  multi1d<DPropagator> qqbar_p(Nt) ; //momentum p
  multi1d<DPropagator> qqbar_q(Nt) ; //momentum q

  qqbar_p = QQbarP[p] ;
  qqbar_q = QQbarQ[q] ;

  if( C.size() != Nt){
        C.resize(Nt);
  }

  DK_6(C, qqbar_p, qqbar_q);

}

#endif
