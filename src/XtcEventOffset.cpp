//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class XtcEventOffset...
//
// Author List:
//      Elliott Slaughter
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "PSXtcInput/XtcEventOffset.h"

//-----------------
// C/C++ Headers --
//-----------------
#include <ostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

namespace PSXtcInput {

//----------------
// Constructors --
//----------------
XtcEventOffset::XtcEventOffset (const std::vector<std::string> &filenames,
                                const std::vector<int64_t> &offsets,
                                const boost::shared_ptr<std::string> &lastBeginCalibCycleDgram)
  : PSEvt::EventOffset()
  , m_filenames(filenames)
  , m_offsets(offsets)
  , m_lastBeginCalibCycleDgram(lastBeginCalibCycleDgram)
{
}

//--------------
// Destructor --
//--------------
XtcEventOffset::~XtcEventOffset ()
{
}

std::vector<std::string>
XtcEventOffset::filenames() const
{
  return m_filenames;
}

std::vector<int64_t>
XtcEventOffset::offsets() const
{
  return m_offsets;
}

const boost::shared_ptr<std::string>
XtcEventOffset::lastBeginCalibCycleDgram() const
{
  return m_lastBeginCalibCycleDgram;
}

/// Dump object in human-readable format
void
XtcEventOffset::print(std::ostream& os) const
{
  os << "XtcEventOffset(offsets=[";
  for (size_t i = 0; i < m_offsets.size(); i++) {
    os << m_offsets[i];
    if (i + 1 < m_offsets.size()) {
      os << ", ";
    }
  }
  os << "], filenames=[";
  for (size_t i = 0; i < m_filenames.size(); i++) {
    os << m_filenames[i];
    if (i + 1 < m_filenames.size()) {
      os << ", ";
    }
  }
  os << "])";
}

} // namespace PSXtcInput
