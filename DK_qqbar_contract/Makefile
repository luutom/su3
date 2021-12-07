CHROMA_DIR = /p/lscratchd/walkloud/c51/code/install/aztec/chroma-double
CHROMA_DIR = /p/project/cslnpp/slnpp032/QCD/USQCD/JURECA/2020-REBUILD/INSTALL/JUR-DCCPU_GCC/CHROMA

################# DO NOT MODIFY BELOW THIS ####################

CPP       = `${CHROMA_DIR}/bin/chroma-config --cxx`
CPP_FLAGS = `${CHROMA_DIR}/bin/chroma-config --cxxflags` -I..
#CPP_FLAGS = -I/dist/scidac/chroma/HEAD/scalar/include -fargument-noalias-global -O2 -finline-limit=50000 -msse2 -march=opteron -I/dist/scidac/qdp++/HEAD/scalar/include -I/usr/include/libxml2 -I..
LD_FLAGS  = `${CHROMA_DIR}/bin/chroma-config --ldflags`
LIBS      = `${CHROMA_DIR}/bin/chroma-config --libs` -lgomp


.SUFFIXES:  .x .d

OBJ = read_meson_block.o
INC = read_meson_block.h DK_qqbar.h

CSRC = $(OBJ:.o=.cc)

all:: DK_qqbar.x

DK_qqbar.x::  $(OBJ) $(INC)


.o.x:
    ${CPP}  ${LD_FLAGS} $*.o  $(OBJ)  ${LIBS} -o $*


# this includes the dependecy files
-include $(OBJ:.o=.d)

# The next line gets the dependences
.cc.d:
    set -e; $(CPP) -MM $(CPP_FLAGS) $< \
    | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
    [ -s $@ ] || rm -f $@

.cc.o:
    ${CPP} ${CPP_FLAGS}  -c $(CPP_FLAGS)  $*.cc

SOURCE = $(CSRC) $(INC) Makefile DK_qqbar.cc

clean:
    rm *.o *.d DK_qqbar
