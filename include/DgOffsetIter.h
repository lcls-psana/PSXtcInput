#ifndef PSXTCINPUT_DGOFFSETITER_H
#define PSXTCINPUT_DGOFFSETITER_H

#include <stdint.h>

#include "pdsdata/xtc/XtcIterator.hh"

namespace Pds {
  class Xtc;
}

struct OffsetPayload {
  int64_t fileOffset;
  uint32_t extent;
  OffsetPayload(int64_t f, uint32_t e) : fileOffset(f), extent(e) {}
};

class DgOffsetIter : public Pds::XtcIterator {
public:
  enum {Stop, Continue};
  DgOffsetIter(Pds::Xtc* xtc);

  int process(Pds::Xtc* xtc);
// private:
  OffsetPayload payload;
};

#endif // PSXTCINPUT_DGOFFSETITER_H
