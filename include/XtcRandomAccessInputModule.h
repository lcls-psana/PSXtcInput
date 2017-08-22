#ifndef PSXTCINPUT_XTCRANDOMACCESSINPUTMODULE_H
#define PSXTCINPUT_XTCRANDOMACCESSINPUTMODULE_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: XtcRandomAccessInputModule.h 7696 2014-02-27 00:40:59Z salnikov@SLAC.STANFORD.EDU $
//
// Description:
//	Class XtcRandomAccessInputModule.
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <string>

//----------------------
// Base Class Headers --
//----------------------
#include "PSXtcInput/XtcInputModuleBase.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "PSXtcInput/RandomAccess.h"
#include "PSXtcInput/DgramPieces.h"

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
 *  @brief Psana input module for random access to XTC files.
 *
 *  This software was developed for the LCLS project.  If you use all or
 *  part of it, please give an appropriate acknowledgment.
 *
 *  @version $Id: XtcRandomAccessInputModule.h 7696 2017-08-18 00:40:59Z eslaught@SLAC.STANFORD.EDU $
 *
 *  @author Elliott Slaughter
 */

  class XtcRandomAccessInputModule : public XtcInputModuleBase {
public:

  /// Constructor takes the name of the module.
  XtcRandomAccessInputModule (const std::string& name) ;

  /// Method which is called once at the beginning of the job
  virtual void beginJob(Event& evt, Env& env);

  psana::RandomAccess& randomAccess() {return _rax;}

  // Destructor
  virtual ~XtcRandomAccessInputModule () ;

protected:

private:
  std::queue<DgramPieces> _queue;
  RandomAccess _rax;
};

} // namespace PSXtcInput

#endif // PSXTCINPUT_XTCRANDOMACCESSINPUTMODULE_H
