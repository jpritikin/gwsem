#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdarg.h>

#include "openmx.h"

void string_vsnprintf(const char *fmt, va_list orig_ap, std::string &dest)
{
	// Cannot stack allocate std::string due to compiler bug that
	// causes clash with va_list when optimization (-O2) is enabled.
    int size = 100;
    while (1) {
        dest.resize(size);
	va_list ap;
	va_copy(ap, orig_ap);
        int n = vsnprintf((char *)dest.c_str(), size, fmt, ap);
	va_end(ap);
        if (n > -1 && n < size) {
            dest.resize(n);
            return;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
}

std::string string_snprintf(const char *fmt, ...)
{
	std::string str;
	va_list ap;
        va_start(ap, fmt);
	string_vsnprintf(fmt, ap, str);
        va_end(ap);
	return str;
}

void mxThrow(const char* msg, ...)
{
	std::string str;
	va_list ap;
	va_start(ap, msg);
	string_vsnprintf(msg, ap, str);
	va_end(ap);

	throw std::runtime_error(str);
}
