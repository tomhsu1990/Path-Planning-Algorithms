#ifndef BFINTERVAL_H
#define BFINTERVAL_H

/* BFInterval is a light weight version of BigFloat2
 * It is mostly used in root isolation routines (Sturm, etc)
 * It has a special interval [1,0] as error indicator.
 */

#include <vector>

CORE_BEGIN_NAMESPACE

typedef long extLong;
typedef std::pair<BigFloat, BigFloat> BFInterval;
typedef std::vector<BFInterval> BFVecInterval;

// global constant:
const BFInterval INVALID_BFInterval(1,0);

CORE_END_NAMESPACE

#endif
