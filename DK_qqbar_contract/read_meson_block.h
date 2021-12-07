#include "chromabase.h"
#include "meas/hadron/qqbar_w.h"

#ifndef __read_meson_block_h__
#define __read_meson_block_h__

using namespace std;

namespace Chroma {
    class ContainerQQbar{
        private:
            multi2d <DPropagator> meson_block ;
            multi2d <int> mom ;
            string type ;
            string sink ;
        public:
            ContainerQQbar(const string& t, const string& s):type(t),sink(s){}
            ContainerQQbar(){}

        void setDescr(const string& t, const string& s){
            type = t ;
            sink = s ;
            }

        void setSize(int Nmom, int tlen, int dims){
            meson_block.resize(Nmom,tlen);
            mom.resize(Nmom,dims);
            }

        string Type() const  {return type;}
        string Sink() const  {return sink;}
    
        multi1d<int> Mom(int i){ return mom[i];}
    
        void setMom(const multi1d<int>& m,int i){ mom[i] = m ;}

        int NumMom() const { return meson_block.size2() ;}
        int Nt() const {return meson_block[0].size();}

        multi1d<DPropagator> operator[](int i){
            return meson_block[i] ;
            }

        int IndxOfMom(const multi1d<int>& p){
            int i(0);
            bool flag(true);

            while(flag&&(i<meson_block.size2())){
                bool eq(true);
                for(int k(0); (k<mom[i].size1())&&eq;k++)
                eq=eq&&(mom[i][k]==p[k]) ;
                if(eq)
                    flag = false ;
                else
                    i++ ;
                }
            if(i>=meson_block.size2()){
                QDPIO::cerr<<"IndxOfMom: Index error!\n" ;
                QDP_abort(1234);
                }
            return i ;
            }
        multi1d <DPropagator> operator[](const multi1d <int> & p){
            return meson_block[IndxOfMom(p)] ;
            }
    } ;


  void read_meson_block(QDPFileReader& from, ContainerQQbar& m_b) ;

}; // namespace chroma
#endif
