//  $Id$
//  $Log$
//  reads  qqbar propagators contracted at the sink with gamma5
//

#include "chromabase.h"
#include "qdp.h"
#include "util/ft/sftmom.h"
#include "util/ferm/transf.h"
#include "read_meson_block.h"

using namespace std;
using namespace Chroma ;
using namespace QDP;

namespace Chroma {
    void read_meson_block(QDPFileReader& from, ContainerQQbar& m_b){
        XMLReader record_xml;
        string sink, type;
        for(int p(0) ; p<m_b.NumMom();p++){
            multi1d<DPropagator> data(m_b.Nt());
            //find the sink and type we want
            read(from,record_xml,data);
            read(record_xml,"/qqbar_desc/sink",sink);
            read(record_xml,"/qqbar_desc/type",type);
            while(sink != m_b.Sink() || type != m_b.Type()){
                read(from,record_xml,data);
                read(record_xml,"/qqbar_desc/sink",sink);
                read(record_xml,"/qqbar_desc/type",type);
                }
            //Now read the thing
            multi1d <int> p_read ;
            read(record_xml,"/qqbar_desc/mom",p_read);
            m_b.setMom(p_read,p);
            m_b[p] = data ;
            } // loop over momenta
        }
    }; //Chroma namespace
