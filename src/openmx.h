#ifndef _OPENMX_DATA_H_
#define _OPENMX_DATA_H_

#include <vector>
#include <map>
#include <string>
#include <cstring>

#include <Rcpp.h>
using namespace Rcpp;

void mxThrow(const char* msg, ...) __attribute__((format (printf, 1, 2))) __attribute__((noreturn));
#define OOPS mxThrow("%s at %d: oops", __FILE__, __LINE__)

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
	dataPtr() : intData(0) {};
	dataPtr(double *_p) : realData(_p) {};
	dataPtr(int *_p) : intData(_p) {};
	void clear() { realData=0; intData=0; };
};

class ColumnData {
private:
	dataPtr ptr;
  bool owner;
  int minValue; // for count/ordinal only
  int maxValue; // for count/ordinal only
public:
	const char *name;
	ColumnDataType type;
	std::vector<std::string> levelNames;       // factors only

	const char *typeName();
  ColumnData(const char *_name) : owner(false), minValue(1), maxValue(NA_INTEGER),
                                  name(_name), type(COLUMNDATA_INVALID) {}
  ColumnData(const char *_name, ColumnDataType _type, int *col) :
    ptr(col), owner(true), minValue(1), maxValue(NA_INTEGER),
    name(_name), type(_type) {}
  ~ColumnData() { clear(); }
  void clear();
  ColumnData clone() const;
  void setMinValue(int mv) { minValue = mv; }
  void verifyMinValue(int nrows);
  void setZeroMinValue(int rows);
  int getMinValue() const { return minValue; }
  int getMaxValue() const { if (maxValue==NA_INTEGER) OOPS; return maxValue; }
  void setMaxValueFromLevels() { maxValue = minValue + levelNames.size() - 1; }
  void setMaxValueFromData(int nrows);
  int getNumThresholds() const { if (maxValue==NA_INTEGER) OOPS; return maxValue - minValue; }
  int getNumOutcomes() const { return 1 + getNumThresholds(); }
  int *i() { return ptr.intData; }
  double *d() { return ptr.realData; }
  void setOwn(double *_p) { clear(); ptr.realData = _p; owner=true; }
  void setOwn(int *_p) { clear(); ptr.intData = _p; owner=true; }
  void setOwn(dataPtr _p) { clear(); ptr = _p; owner=true; }
  void setBorrow(double *_p) { clear(); ptr.realData = _p; owner=false; }
  void setBorrow(int *_p) { clear(); ptr.intData = _p; owner=false; }
  void setBorrow(dataPtr _p) { clear(); ptr = _p; owner=false; }
  dataPtr steal() { dataPtr ret = ptr; ptr.clear(); return ret; }
};

inline void ColumnData::clear()
{
  if (ptr.intData && owner) {
    if (type == COLUMNDATA_NUMERIC) {
      delete [] ptr.realData;
    } else {
      delete [] ptr.intData;
    }
  }
  ptr.intData = 0;
}
struct cstrCmp {
	bool operator() (const char *s1, const char *s2) const
	{ return strcmp(s1,s2) < 0; }
};

typedef std::map< const char *, int, cstrCmp > ColMapType;

std::string string_snprintf(const char *fmt, ...) __attribute__((format (printf, 1, 2)));

#endif
