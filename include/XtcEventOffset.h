#ifndef PSXTCINPUT_XTCEVENTOFFSET_H
#define PSXTCINPUT_XTCEVENTOFFSET_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class XtcEventOffset.
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------

//----------------------
// Base Class Headers --
//----------------------
#include "PSEvt/EventOffset.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

namespace PSXtcInput {

/**
 *  @ingroup PSXtcInput
 *
 *  @brief Implementation of the EventOffset interface for XTC events.
 *
 *  This software was developed for the LCLS project.  If you use all or
 *  part of it, please give an appropriate acknowledgment.
 *
 *  @version $Id$
 *
 *  @author Elliott Slaughter
 */

class XtcEventOffset : public PSEvt::EventOffset {
public:

  // Default constructor
  XtcEventOffset (const std::vector<std::string> &filenames,
                  const std::vector<int64_t> &offsets,
                  int64_t configureOffset,
                  int64_t beginRunOffset,
                  const std::string &lastBeginCalibCycleFilename,
                  int64_t lastBeginCalibCycleOffset) ;

  // Destructor
  ~XtcEventOffset () ;

  virtual std::vector<std::string> filenames() const;
  virtual std::vector<int64_t> offsets() const;
  virtual int64_t configureOffset() const;
  virtual int64_t beginRunOffset() const;
  virtual std::string lastBeginCalibCycleFilename() const;
  virtual int64_t lastBeginCalibCycleOffset() const;

  /// Dump object in human-readable format
  virtual void print(std::ostream& os) const;

protected:

private:

  // Data members
  std::vector<std::string> m_filenames;
  std::vector<int64_t> m_offsets;
  int64_t m_configureOffset;
  int64_t m_beginRunOffset;
  std::string m_lastBeginCalibCycleFilename;
  int64_t m_lastBeginCalibCycleOffset;
};

} // namespace PSXtcInput

#endif // PSXTCINPUT_XTCEVENTOFFSET_H