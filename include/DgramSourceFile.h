#ifndef PSXTCINPUT_DGRAMSOURCEFILE_H
#define PSXTCINPUT_DGRAMSOURCEFILE_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class DgramSourceFile.
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <string>
#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/scoped_ptr.hpp>

//----------------------
// Base Class Headers --
//----------------------
#include "PSXtcInput/IDatagramSource.h"
#include "psana/Configurable.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
namespace XtcInput {
  class DgramQueue;
}

//		---------------------
// 		-- Class Interface --
//		---------------------

namespace PSXtcInput {

/// @addtogroup PSXtcInput

/**
 *  @ingroup PSXtcInput
 *
 *  @brief Implementation of IDatagramSource interface which reads data from input files.
 *
 *  @note This software was developed for the LCLS project.  If you use all or 
 *  part of it, please give an appropriate acknowledgment.
 *
 *  @version $Id$
 *
 *  @author Andy Salnikov
 */

class DgramSourceFile : public IDatagramSource, public psana::Configurable {
public:

  // Constructor takes the name of the input module, used for accessing module 
  // configuration parameters. 
  DgramSourceFile(const std::string& name);

  // Destructor
  virtual ~DgramSourceFile();

  /**
   *   Initialization method for datagram source, this is typically called
   *   in beginJob() method and it may contain initialization code which 
   *   cannot be executed during construction of an instance.
   */ 
  virtual void init();

  /**
   *  @brief Return next datagram(s) from the source.
   *
   *  This method returns two sets of datagrams - eventDg is the set of histograms
   *  belonging to the next event, nonEventDg is the set of datagrams which has some
   *  other data (like EPICS) which is needed for correct interpretation of current
   *  event. Currently eventDg should contain one datagram but potentially in the
   *  future we may start event building in offline and that list can grow longer.
   *  It nonEventDg is non-empty then it has to be processed first as those datagram
   *  should come from earlier time than eventDg and eventDg may contain data that
   *  overrides data in nonEventDg (e.g. some EPICS PV data may be contained in both
   *  nonEventDg and eventDg).
   *
   *  This method will called repeatedly until it returns false.
   *
   *  @param[out] eventDg    returned set of datagrams for current event
   *  @param[out] nonEventDg returned set of datagrams containing other information.
   *  @return false if there are no more events, both eventDg and nonEventDg will be empty in this case.
   */
  virtual bool next(std::vector<XtcInput::Dgram>& eventDg, std::vector<XtcInput::Dgram>& nonEventDg);

protected:

  bool sameEvent(const XtcInput::Dgram &eventDg, const XtcInput::Dgram &otherDg) const;

private:

  boost::scoped_ptr<XtcInput::DgramQueue> m_dgQueue;  ///< Input datagram queue
  boost::scoped_ptr<boost::thread> m_readerThread;    ///< Thread which does datagram reading
  std::vector<std::string> m_fileNames;               ///< List of file names/datasets to read data from
  int m_firstControlStream;                           ///< Starting index of control streams

};

} // namespace PSXtcInput

#endif // PSXTCINPUT_DGRAMSOURCEFILE_H
