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
                                int64_t configureOffset,
                                int64_t beginRunOffset,
                                const std::string &lastBeginCalibCycleFilename,
                                int64_t lastBeginCalibCycleOffset)
  : PSEvt::EventOffset()
  , m_filenames(filenames)
  , m_offsets(offsets)
  , m_configureOffset(configureOffset)
  , m_beginRunOffset(beginRunOffset)
  , m_lastBeginCalibCycleFilename(lastBeginCalibCycleFilename)
  , m_lastBeginCalibCycleOffset(lastBeginCalibCycleOffset)
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

int64_t
XtcEventOffset::configureOffset() const
{
  return m_configureOffset;
}

int64_t
XtcEventOffset::beginRunOffset() const
{
  return m_beginRunOffset;
}

std::string
XtcEventOffset::lastBeginCalibCycleFilename() const
{
  return m_lastBeginCalibCycleFilename;
}

int64_t
XtcEventOffset::lastBeginCalibCycleOffset() const
{
  return m_lastBeginCalibCycleOffset;
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
