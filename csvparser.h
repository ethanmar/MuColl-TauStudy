#ifndef __CSVPARSE_H_2001_06_07__
#define __CSVPARSE_H_2001_06_07__

#include <string>
using namespace std;


class CSVParser {
 private:
  string m_sData;
  string::size_type m_nPos;
  void SkipSpaces(void);
 public:
  CSVParser();
  const CSVParser & operator << (const string &sIn);
  const CSVParser & operator << (const char *sIn);
  CSVParser & operator >> (int &nOut);
  CSVParser & operator >> (double &nOut);
  CSVParser & operator >> (string &sOut);
};

#endif