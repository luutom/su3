#CHROMA=${EBROOTCHROMA}
#CHROMA=/p/project/cslnpp/slnpp032/QCD/USQCD/JURECA/CHROMA/chroma/CHROMA_QPHIX_KNL_INSTALL/
CHROMA=/p/project/cslnpp/slnpp032/QCD/USQCD/JURECA/2020-REBUILD/INSTALL/CHROMA

# ./myChroma -i su3.ini.xml -by 2 -bz 2 -pxy 1 -pxyz1 -minct 1 -c 256 -sy 1 -sz 1

#CHROMA=/p/project/cslnpp/slnpp032/QCD/USQCD/JURECA/2020-REBUILD/INSTALL/JUR-DCCPU_GCC/CHROMA
CONFIG=$(CHROMA)/bin/chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I. 
LDFLAGS=$(shell $(CONFIG) --ldflags)
LIBS=$(shell $(CONFIG) --libs)

HDRS=su3_meas.h Cqbar_qqbar.h spin_matrix.h permutations.h mom_project.h readwrite.h

#HDRS=su3_meas.h Cqbar_qqbar.h kn_meas.h kn_sys.h meson_meas.h nucleon_meas.h spin_matrix.h permutations.h udd.h mom_project.h readwrite.h nucleon.h meson.h myMultiD.h

OBJS=chroma.o \
     su3_meas.o \
     spin_matrix.o \
     permutations.o \
     mom_project.o \
     readwrite.o  \
     Cqbar_qqbar.o

myChroma: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

gammaDriver: gammaDriver.o spin_matrix.o
	$(CXX) -o $@ $(CXXFLAGS) gammaDriver.o spin_matrix.o $(LDFLAGS) $(LIBS)

archive:
	tar -zcvf kn.tar.gz *.cc Makefile kn.ini.xml  *.h cfgs/ props/

%.o: %.cc $(HDRS)
	$(CXX) -w $(CXXFLAGS) -c $< 

clean:
	rm -rf $(OBJS) *~ *.optrpt

distclean:
	rm -rf *.exe *.o *~ myChroma
