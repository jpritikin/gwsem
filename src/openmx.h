#ifndef _OPENMX_DATA_H_
#define _OPENMX_DATA_H_

#include <vector>
#include <map>
#include <string>
#include <cstring>

#include <Rcpp.h>
using namespace Rcpp;

enum ColumnDataType {
	COLUMNDATA_INVALID,
	COLUMNDATA_ORDERED_FACTOR,
	COLUMNDATA_UNORDERED_FACTOR,
	COLUMNDATA_INTEGER,
	COLUMNDATA_NUMERIC
};

union dataPtr {
	double *realData;
	int *intData;
	dataPtr(double *_p) : realData(_p) {};
	dataPtr(int *_p) : intData(_p) {};
	void clear() { realData=0; intData=0; };
};

struct ColumnData {
	const char *name;
	ColumnDataType type;
	dataPtr ptr;
	std::vector<std::string> levels;       // factors only

	const char *typeName();
};

struct cstrCmp {
	bool operator() (const char *s1, const char *s2) const
	{ return strcmp(s1,s2) < 0; }
};

typedef std::map< const char *, int, cstrCmp > ColMapType;

std::string string_snprintf(const char *fmt, ...) __attribute__((format (printf, 1, 2)));

#endif
