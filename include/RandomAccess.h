#ifndef PSXTCINPUT_RANDOMACCESS_H
#define PSXTCINPUT_RANDOMACCESS_H

#include <boost/cstdint.hpp>
#include <vector>
#include <queue>
#include <string>

#include "psana/RandomAccess.h"
#include "psana/Configurable.h"
#include "PSXtcInput/DgramPieces.h"

namespace PSXtcInput {

/// @addtogroup PSXtcInput

/**
 *  @ingroup PSXtcInput
 *
 *  @brief Interface to allow XTC file random access.
 *
 *  @version $Id: RandomAccess.h 7696 2017-08-18 00:40:59Z eslaught@SLAC.STANFORD.EDU $
 *
 *  @author Elliott Slaughter
 */

class RandomAccessRun;
class RunMap;

class RandomAccess : public psana::RandomAccess, public psana::Configurable {
public:
  RandomAccess(const std::string& name, std::queue<DgramPieces>& queue);
  ~RandomAccess();
  int jump(const std::vector<std::string>& filenames, const std::vector<int64_t> &offsets, const std::string &lastBeginCalibCycleDgram, uintptr_t runtime, uintptr_t ctx);
  void setrun(int run);
  const std::vector<unsigned>& runs();
private:
  std::vector<std::string> _fileNames;
  std::queue<DgramPieces>& _queue;
  RandomAccessRun*         _raxrun;
  RunMap*                  _rmap;
  int                      _run;
};

}

#endif // PSXTCINPUT_RANDOMACCESS_H
