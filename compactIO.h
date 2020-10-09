#ifndef COMPACTIO_H
#define COMPACTIO_H

#include <deque>
#include "Event.h"
#include "GolemMCSet.h"

namespace golemfit {

namespace jason {
  void splatData(const std::string& filename, const std::deque<Event>& exp, const std::deque<Event>& sim);
  void unsplatData(const std::string& filename, std::deque<Event>& exp, std::deque<Event>& sim);
}

namespace dump {
  void splatData(const std::string& filename, const uint32_t progChecksum, const std::deque<Event>& exp, const std::deque<Event>& sim);
  void unsplatData(const std::string& filename, const uint32_t progChecksum, std::deque<Event>& exp, std::deque<Event>& sim);
}

}

#endif //COMPACTIO_H
