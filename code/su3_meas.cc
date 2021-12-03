#include "su3_meas.h"
#include "chroma.h"
#include "spin_matrix.h"
#include <iostream>
#include <string>

using namespace std;
using namespace QDP;

namespace Chroma 
{ 
  // Fill in the blanks
  namespace SU3Env 
  { 

    // Function to register with a factory
    // This is boilerplate stuff
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      // Create a Params from the XML
      // Use it to create a SU3 system
      return new SU3( Params(xml_in, path) );
    }

    // The name of my measurement for the XML file
    // Change this for each measurement
    const std::string name = "SU3_SYSTEM";

    // Register the measurement in the measurement factory
    namespace { 
      bool registered = false;
    }

    bool registerAll()
    {
      bool success = true;
      if (! registered) { 
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }

    // Param stuff
    Params::Params() { su3_system="";frequency = 0; gauge_id=""; light_prop_id=""; heavy_prop_id=""; xml_file=""; filedir=""; consistencyTests = false; Gamma1=15; Gamma2=15; r15 = true; r6 = true; r3bar = false; }
    
    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
	{
	  XMLReader paramtop(xml_in, path);

	  if (paramtop.count("Frequency") == 1)
	    read(paramtop, "Frequency", frequency);
	  else
	    frequency = 1;
      
	  /* Name of system to do contractions */
	  read(paramtop, "Name", su3_system);

	  // Read in the file directory
	  read(paramtop,"Param/file_dir",filedir);

	  // Read in the output propagator/source configuration info for light quark
	  read(paramtop, "NamedObject/light_prop_id", light_prop_id);

          // Read in the output propagator/source configuration info for heavy quark
          read(paramtop, "NamedObject/heavy_prop_id", heavy_prop_id);

	  // Get either the gauge_id in the XML or the ID of the default
	  // field of no explicit ID exists.
	  read(paramtop, "NamedObject/gauge_id", gauge_id);

	  // get the indexes for the relevan gamma matrices (default is 15 (=gamma_5) for both cases)
	  read(paramtop,"Param/Gamma1",Gamma1);
	  read(paramtop,"Param/Gamma2",Gamma2);

	  if( paramtop.count("xml_file") != 0 ) { 
	    read(paramtop, "xml_file", xml_file);
	  }
	  else { 
	    xml_file == "";
	  }
	
	}
      catch(const std::string& e) 
	{
	  QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	  QDP_abort(1);
	}

      try {
	XMLReader paramtop(xml_in, path);
      
	// Read in the file directory
	read(paramtop,"Param/ConsistencyTests",consistencyTests);
	if(consistencyTests){
	  //	  QDPIO::cout << "Performing consistency tests" << endl;
	} else {
	  //	  QDPIO::cout << "Not performing consistency tests" << endl;
	};
      }catch(const std::string& e) {
	QDPIO::cout << "No consistency-tests option stated in " << e << " : will NOT perform any consistency tests" << endl;
      }
      
      try {
	XMLReader paramtop(xml_in, path);
      
	// Read in the file directory
	read(paramtop,"Param/Rep15",r15);
	if(r15){
	  QDPIO::cout << "Calculating [15] system" << endl;
	} else {
	  //	  QDPIO::cout << "Not performing consistency tests" << endl;
	};
      }catch(const std::string& e) {
	QDPIO::cout << "No [15] rep option stated in " << e << " : will NOT perform calculations of this system" << endl;
      }

      try {
	XMLReader paramtop(xml_in, path);
      
	// Read in the file directory
	read(paramtop,"Param/Rep6",r6);
	if(r6){
	  QDPIO::cout << "Calculating [6] system" << endl;
	} else {
	  //	  QDPIO::cout << "Not performing consistency tests" << endl;
	};
      }catch(const std::string& e) {
	QDPIO::cout << "No [6] rep option stated in " << e << " : will NOT perform calculations of this system" << endl;
      }

      try {
	XMLReader paramtop(xml_in, path);
      
	// Read in the file directory
	read(paramtop,"Param/Rep3bar",r3bar);
	if(r3bar){
	  QDPIO::cout << "Calculating [3b] system" << endl;
	} else {
	  //	  QDPIO::cout << "Not performing consistency tests" << endl;
	};
      }catch(const std::string& e) {
	QDPIO::cout << "No [3b] rep option stated in " << e << " : will NOT perform calculations of this system" << endl;
      }
      // Read your own parameters  from Param/ here
    }


    void
    Params::write(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Write out the file directory
      QDP::write(xml_out,"file_dir",filedir);

      // Write out the output propagator/source configuration info
      QDP::write(xml_out, "light_prop_id", light_prop_id);

      // Write out the output propagator/source configuration info
      QDP::write(xml_out, "heavy_prop_id", heavy_prop_id);

      QDP::write(xml_out, "gauge_id", gauge_id);

      if( xml_file != "" ){ 
	QDP::write(xml_out, "xml_file", xml_file);
      }


      pop(xml_out);
    }
  }

  void 
  SU3::operator()(long unsigned int update_no,
			  XMLWriter& xml_out) 
  {

    // This bit merely supports providing an external xml_file 
    // for this measurement
    if ( params.xml_file == "" ) { 
      
      func( update_no, xml_out );
    }
    else { 

      // Hash an XML file name from the user specified string
      // and the update number
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      // IN the global output, make a note that the output went
      // to this separate XML file
      push(xml_out, "my_meas");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
  }

  void 
  SU3::func(unsigned long int update_no, XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch measure_time;
    measure_time.reset();
    measure_time.start();


    // Test that the gauge configuration and the propagator we need
    // exist in the map.
    XMLBufferWriter gauge_xml;
    XMLBufferWriter light_prop_xml;
    XMLBufferWriter heavy_prop_xml;
    try
      {
	// Try and get at the gauge field if it doesn't exist 
	// an exception will be thrown.
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);
	TheNamedObjMap::Instance().get(params.gauge_id).getRecordXML(gauge_xml);

	// Do the same with the light propagator 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.light_prop_id);

	// Do the same with the heavy propagator
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.heavy_prop_id);

	// Get its record XML
	TheNamedObjMap::Instance().get(params.light_prop_id).getRecordXML(light_prop_xml);

	// Get its record XML
	TheNamedObjMap::Instance().get(params.heavy_prop_id).getRecordXML(heavy_prop_xml);

      }
    catch( std::bad_cast ) 
      {

	// If the object exists and is of the wrong type we will 
	// land in this catch.
	QDPIO::cerr << params.su3_system << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
    catch (const string& e) 
      {
	// If the object is not in the map we will land in this 
	// catch
	QDPIO::cerr << params.su3_system << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }

    // If we got here, then the gauge field and prop are both in 
    // the map. Their XML will have been captured.
    // Let us bind the references to a local name so 
    // we don't have to type the long lookup string every time.
    //
    // Note const here means we can't change the field
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);

    const LatticePropagator& light_prop = 
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.light_prop_id);

    const LatticePropagator& heavy_prop =
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.heavy_prop_id);

    XMLReader input_light_prop_xml;
    TheNamedObjMap::Instance().get(params.light_prop_id).getRecordXML(input_light_prop_xml);

    XMLReader input_heavy_prop_xml;
    TheNamedObjMap::Instance().get(params.heavy_prop_id).getRecordXML(input_heavy_prop_xml);

    // I assume that the heavy and light propagators have the same source locations  (!! might want to change this later !!)

    multi1d<int> srce_pt;
    read(input_light_prop_xml,"/SinkSmear/PropSource/Source/t_srce",srce_pt);

    multi1d<int> latsize;
    read(input_light_prop_xml,"/SinkSmear/Config_info/Params/HMCTrj/nrow",latsize);

    std::string traj;
    read(input_light_prop_xml,"/SinkSmear/Config_info/Params/MCControl/StartUpdateNum",traj);

    std::string srcesmear;
    read(input_light_prop_xml,"/SinkSmear/PropSource/Source/SourceType",srcesmear);

    std::string sinksmear;
    read(input_light_prop_xml,"/SinkSmear/PropSink/Sink/SinkType",sinksmear);

    std::string basename="_"+srcesmear+"_"+sinksmear+"_"+traj;
    std::string dirname=params.filedir;
    int Nt=latsize[3];

    // Boilerplate stuff to the output XML
    push(xml_out, "my_measurement");
    write(xml_out, "update_no", update_no);

    // Write info about the program
    proginfo(xml_out);

    // Write out the input
    params.write(xml_out, "Input");

    multi1d<int> pN(4);
  
    pN[0]=0;pN[1]=0;pN[2]=0;pN[3]=0;
    
    
    QDPIO::cout << "Light propagator ID is " << params.light_prop_id << endl;
    QDPIO::cout << "Heavy propagator ID is " << params.heavy_prop_id << endl;

    Cqbar_qqbar_system Dpi;
    int Gamma1,Gamma2,n2,nstates;
    std::string nshell;
    std::vector<std::vector<int>> pp;

    Gamma1 = params.Gamma1;
    Gamma2 = params.Gamma2;

    QDPIO::cout << params.su3_system <<" Beginning " << endl;
    Dpi.doChecks = params.consistencyTests;
    Dpi.initialize(srce_pt,pN,latsize,dirname,basename);

    // ok, I hardwire in the different n2 shells
    for(n2=0;n2<=5;n2++) {
      switch(n2) {
      case 0:
	nstates = 1;
	pp.resize(nstates);
	for (int nj = 0;nj < nstates; nj++) pp[nj].resize(4);
	pp[0][0]=0;pp[0][1]=0;pp[0][2]=0;pp[0][3]=0;
	nshell = "_N0";
	break;
      case 1:
	nstates = 3;
	pp.resize(nstates);
	for (int nj = 0;nj < nstates; nj++) pp[nj].resize(4);
	pp[0][0]=1;pp[0][1]=0;pp[0][2]=0;pp[0][3]=0;
	pp[1][0]=0;pp[1][1]=1;pp[1][2]=0;pp[1][3]=0;
	pp[2][0]=0;pp[2][1]=0;pp[2][2]=1;pp[2][3]=0;
	nshell = "_N1";
	break;
      case 2:
	nstates = 6;
	pp.resize(nstates);
	for (int nj = 0;nj < nstates; nj++) pp[nj].resize(4);
	pp[0][0]= 1;pp[0][1]= 1;pp[0][2]= 0;pp[0][3]=0;
	pp[1][0]= 1;pp[1][1]=-1;pp[1][2]= 0;pp[1][3]=0;
	pp[2][0]= 1;pp[2][1]= 0;pp[2][2]= 1;pp[2][3]=0;
	pp[3][0]= 1;pp[3][1]= 0;pp[3][2]=-1;pp[3][3]=0;
	pp[4][0]= 0;pp[4][1]= 1;pp[4][2]= 1;pp[4][3]=0;
	pp[5][0]= 0;pp[5][1]= 1;pp[5][2]=-1;pp[5][3]=0;	
	nshell = "_N2";
	break;
      case 3:
	nstates = 4;
	pp.resize(nstates);
	for (int nj = 0;nj < nstates; nj++) pp[nj].resize(4);
	pp[0][0]= 1;pp[0][1]= 1;pp[0][2]= 1;pp[0][3]=0;
	pp[1][0]= 1;pp[1][1]=-1;pp[1][2]= 1;pp[1][3]=0;
	pp[2][0]= 1;pp[2][1]= 1;pp[2][2]=-1;pp[2][3]=0;
	pp[3][0]= 1;pp[3][1]=-1;pp[3][2]=-1;pp[3][3]=0;
	nshell = "_N3";
	break;
      case 4:
	nstates = 3;
	pp.resize(nstates);
	for (int nj = 0;nj < nstates; nj++) pp[nj].resize(4);
	pp[0][0]= 2;pp[0][1]= 0;pp[0][2]= 0;pp[0][3]=0;
	pp[1][0]= 0;pp[1][1]= 2;pp[1][2]= 0;pp[1][3]=0;
	pp[2][0]= 0;pp[2][1]= 0;pp[2][2]= 2;pp[2][3]=0;
	nshell = "_N4";
	break;
      case 5:
	nstates = 12;
	pp.resize(nstates);
	for (int nj = 0;nj < nstates; nj++) pp[nj].resize(4);
	pp[0][0]=  2;pp[0][1]=  1;pp[0][2]=  0;pp[0][3]=0;
	pp[1][0]=  2;pp[1][1]= -1;pp[1][2]=  0;pp[1][3]=0;
	pp[2][0]=  2;pp[2][1]=  0;pp[2][2]=  1;pp[2][3]=0;
	pp[3][0]=  2;pp[3][1]=  0;pp[3][2]= -1;pp[3][3]=0;
	pp[4][0]=  0;pp[4][1]=  2;pp[4][2]=  1;pp[4][3]=0;
	pp[5][0]=  0;pp[5][1]=  2;pp[5][2]= -1;pp[5][3]=0;
	pp[6][0]=  1;pp[6][1]=  2;pp[6][2]=  0;pp[6][3]=0;
	pp[7][0]=  1;pp[7][1]= -2;pp[7][2]=  0;pp[7][3]=0;
	pp[8][0]=  1;pp[8][1]=  0;pp[8][2]=  2;pp[8][3]=0;
	pp[9][0]=  1;pp[9][1]=  0;pp[9][2]= -2;pp[9][3]=0;
	pp[10][0]= 0;pp[10][1]= 1;pp[10][2]= 2;pp[10][3]=0;
	pp[11][0]= 0;pp[11][1]= 1;pp[11][2]=-2;pp[11][3]=0;	
	nshell = "_N5";
	break;
      default: // default is always zero momentum
	nstates = 1;
	pp.resize(nstates);
	for (int nj = 0;nj < nstates; nj++) pp[nj].resize(4);
	pp[0][0]=0;pp[0][1]=0;pp[0][2]=0;pp[0][3]=0;
	nshell = "_N0";
      }
      for (int ns = 0; ns < nstates; ns++){
	pN[0]=pp[ns][0];pN[1]=pp[ns][1];pN[2]=pp[ns][2];pN[3]=pp[ns][3];
	Dpi.reset_momentum(srce_pt,pN);
    
	// do calculations
	Dpi.calcDirect(light_prop,heavy_prop,Gamma1,Gamma2);
	if(params.r6 || params.r15) Dpi.calcExchange(light_prop,heavy_prop,Gamma1,Gamma2);
	if(params.r15) Dpi.calc_15_rep();
	if(params.r6) Dpi.calc_6_rep();
	if(params.r3bar) {
	  Dpi.calcExchange3bar(light_prop,heavy_prop,Gamma1,Gamma2);
	  Dpi.calc_3bar_rep();
	}

	// write out the correlators
	if (Gamma1 == Gamma2) Dpi.write_qq_and_Qq(Gamma1,Gamma2,nshell);// only write out individual correlators if Gamma1 == Gamma2
	if(params.r6 && !params.r15) Dpi.write_6(Gamma1,Gamma2,nshell);
	if(params.r15 && !params.r6) Dpi.write_15(Gamma1,Gamma2,nshell);
	if(params.r15 && params.r6) Dpi.write_15_and_6(Gamma1,Gamma2,nshell);
	if(params.r3bar) {
	  Dpi.write_direct(Gamma1,Gamma2,nshell);
	  Dpi.write_3bar(Gamma1,Gamma2,nshell);
	}
      }
    }
    
    //    Dpi.calcDirect_P0(light_prop,heavy_prop,Gamma1,Gamma2);
    //    Dpi.calcExchange_P0(light_prop,heavy_prop,Gamma1,Gamma2);
    //    Dpi.calc_15_rep();
    //    Dpi.calc_6_rep();
    //    Dpi.write_15_and_6(Gamma1,Gamma2,"_P0_");

    // free everything up
    Dpi.free();


    pop(xml_out);
    measure_time.stop();
    QDPIO::cout << params.su3_system << ": total time = "
		<< measure_time.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << params.su3_system << ": ran successfully" << endl;
    END_CODE();
  } 

};

