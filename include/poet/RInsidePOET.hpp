#ifndef RPOET_H_
#define RPOET_H_

#include <RInside.h>

class RInsidePOET : public RInside {
public:
  static RInsidePOET &getInstance() {
    static RInsidePOET instance;

    return instance;
  }

  RInsidePOET(RInsidePOET const &) = delete;
  void operator=(RInsidePOET const &) = delete;

private:
  RInsidePOET() : RInside(){};
};

#endif // RPOET_H_
