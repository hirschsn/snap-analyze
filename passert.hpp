// For license details see LICENSE.

#pragma once

#include <cstdio>
#include <cstdlib>

[[noreturn]] void __passert_fail(const char *expr, const char *file, int line, const char *function)
{
    std::fprintf(stderr, "p_assert assertion failed: `%s' in %s:%i (%s)", expr, file, line, function);
    std::abort();
}

/** Assert-equivalent that is *not* a no-op if NDEBUG is set.
 */
#define p_assert(cond)                                                          \
    ((cond)? (void) 0: __passert_fail(#cond, __FILE__, __LINE__, __FUNCTION__))

