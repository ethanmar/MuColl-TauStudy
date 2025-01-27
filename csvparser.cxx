#include "MyTauFinder/csvparser.h"
#include <cstdlib>
#include <iostream>
using namespace std;

CSVParser::CSVParser() : m_sData(""), m_nPos(0) {}

void CSVParser::SkipSpaces(void) {
  while (m_nPos < m_sData.length() && m_sData[m_nPos] == ' ')
    m_nPos++;
}

const CSVParser &CSVParser::operator<<(const string &sIn) {
  this->m_sData = sIn;
  this->m_nPos = 0;
  return *this;
}

const CSVParser &CSVParser::operator<<(const char *sIn) {
  this->m_sData = sIn;
  this->m_nPos = 0;
  return *this;
}

CSVParser &CSVParser::operator>>(int &nOut) {
  string sTmp = "";
  SkipSpaces();
  while (m_nPos < m_sData.length() && m_sData[m_nPos] != ',')
    sTmp += m_sData[m_nPos++];

  m_nPos++; // skip past comma
  nOut = atoi(sTmp.c_str());
  return *this;
}

CSVParser &CSVParser::operator>>(double &nOut) {
  string sTmp = "";
  SkipSpaces();
  while (m_nPos < m_sData.length() && m_sData[m_nPos] != ',')
    sTmp += m_sData[m_nPos++];

  m_nPos++; // skip past comma
  nOut = atof(sTmp.c_str());
  return *this;
}

CSVParser &CSVParser::operator>>(string &sOut) {
  bool bQuotes = false;
  sOut = "";
  SkipSpaces();

  // Jump past first " if necessary
  if (m_nPos < m_sData.length() && m_sData[m_nPos] == '"') {
    bQuotes = true;
    m_nPos++;
  }

  while (m_nPos < m_sData.length()) {
    if (!bQuotes && m_sData[m_nPos] == ',')
      break;
    if (bQuotes && m_sData[m_nPos] == '"') {
      if (m_nPos + 1 >= m_sData.length() - 1)
        break;
      if (m_sData[m_nPos + 1] == ',')
        break;
    }
    sOut += m_sData[m_nPos++];
  }

  // Jump past last " if necessary
  if (bQuotes && m_nPos < m_sData.length() && m_sData[m_nPos] == '"')
    m_nPos++;

  // Jump past , if necessary
  if (m_nPos < m_sData.length() && m_sData[m_nPos] == ',')
    m_nPos++;

  return *this;
}