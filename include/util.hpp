#ifndef UTIL
#define UTIL

/***********************************************
 * Constant configuration:
 ***********************************************/
/* Constant value of PI */
#define PI 3.14159265358979323846
/* Constant value of sqrt(2) */
#define S2I 0.70710678118654752440
/* Constant value of 0.5 */
#define HALF 0.5

#define DEFAULT_REPETITIONS 1024

namespace NWQSim
{
    /* Basic data type for indices */
    using IdxType = long long int;
    /* Basic data type for value */
    using ValType = double;
} // namespace NWQSim

#endif // UTIL