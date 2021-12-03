/* 
 * File:   su3_meas.h
 * Author: tomluu
 *
 * Created on December 6, 2013, 10:07 AM
 */

#ifndef SU3_MEAS_H
#define	SU3_MEAS_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "Cqbar_qqbar.h"

namespace Chroma 
{ 
  // A namespace for this particular   measurement
  namespace SU3Env 
  {
    extern const std::string name;
    bool registerAll();
  

    //! Parameter structure
    struct Params 
    {
      // Default constructor -- forward declaration
      Params();

      // Construct from XML -- forward declaration
      Params(XMLReader& xml_in, const std::string& path);

      // Write myself out
      void write(XMLWriter& xml_out, const std::string& path);

      // How often should I measure me in an HMC run
      unsigned long frequency;

      // Various parameters taken from XML file
      std::string su3_system;
      std::string gauge_id;
      std::string light_prop_id;
      std::string heavy_prop_id;
      std::string filedir;
      int Gamma1, Gamma2;
      bool consistencyTests; /* not mandatory */
      bool r6,r15,r3bar;

      std::string xml_file; // Support output to own private XML File
    }; // struct
  }; // namespace SU3Env

  class SU3 : public AbsInlineMeasurement 
  {
  public:
    // Constructor -- default -- do nothing
    ~SU3() {}

    // Constructor -- from param struct -- copy param struct inside me
    SU3(const SU3Env::Params& p) : params(p) {}

    // Constructor -- copy constructor -- copy param struct from argument
    SU3(const SU3& p) : params(p.params) {}

    // Boiler plate
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 


  private:

    void func(const unsigned long update_no,XMLWriter& xml_out);
    
    SU3Env::Params params;
  };

};

#endif	/* SU3_MEAS_H */

