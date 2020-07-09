
//##########################################################################
//# This software is part of the Treescaper i
//# -- Version 0.1   
//# Copyright (C) 2010 Wen Huang
//# 
//# This program is free software; you can redistribute it and/or
//# modify it under the terms of the GNU General Public License
//# as published by the Free Software Foundation; either version 2
//# of the License, or (at your option) any later version.
//#
//# This program is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details. 
//# http://www.gnu.org/copyleft/gpl.html 
//##########################################################################

// wstring.h
// Definition of a String class
//           by whuang Nov/03/2008
#ifndef Q_DEBUG_STREAM_H
#define Q_DEBUG_STREAM_H

#include <iostream>
#include <streambuf>
#include <string>
#include "wfile.h"
#include "QTextCursor"
#include "QScrollBar"
#include "QMutexLocker"

#include "qtextedit.h"

class QDebugStream : public std::basic_streambuf<char>
{
public:
 QDebugStream(std::ostream &stream, QTextEdit* text_edit1,
              QTextEdit* text_edit2, QTextEdit* text_edit3,
              QTextEdit* text_edit4, QTextEdit* text_edit5,
              QTextEdit* text_edit6, char *logfilename) : m_stream(stream)
 {
     log_window1 = text_edit1;
     log_window2 = text_edit2;
     log_window3 = text_edit3;
     log_window4 = text_edit4;
     log_window5 = text_edit5;
     log_window6 = text_edit6;
     logfname = logfilename;
     m_old_buf = stream.rdbuf();
     stream.rdbuf(this);
     if(logfname != NULL)
     {
         logfile.open(logfname);
         logfile.clean();
     }
 }
 ~QDebugStream()
 {
  // output anything that is left
  if (!m_string.empty())
  {
      log_window1->append(m_string.c_str());
      log_window2->append(m_string.c_str());
      log_window3->append(m_string.c_str());
      log_window4->append(m_string.c_str());
      log_window5->append(m_string.c_str());
      log_window6->append(m_string.c_str());
//      log_window1->moveCursor(QTextCursor::End);
//      log_window2->moveCursor(QTextCursor::End);
//      log_window3->moveCursor(QTextCursor::End);
//      log_window4->moveCursor(QTextCursor::End);
//      log_window5->moveCursor(QTextCursor::End);
//      log_window6->moveCursor(QTextCursor::End);
      if(logfile.is_open())
          logfile << m_string.c_str() << std::endl;
  }

  m_stream.rdbuf(m_old_buf);
  if(logfile.is_open())
      logfile.close();
 }

protected:
 virtual int_type overflow(int_type v)
 {
     QMutexLocker locker(&mutex1);
  if (v == '\n')
  {
      log_window1->append(m_string.c_str());
      log_window2->append(m_string.c_str());
      log_window3->append(m_string.c_str());
      log_window4->append(m_string.c_str());
      log_window5->append(m_string.c_str());
      log_window6->append(m_string.c_str());
      if(logfile.is_open())
          logfile << m_string.c_str() << std::endl;
      m_string.erase(m_string.begin(), m_string.end());
  }
  else
  {
   m_string += v;
//   log_window1->moveCursor(QTextCursor::End);
//   log_window2->moveCursor(QTextCursor::End);
//   log_window3->moveCursor(QTextCursor::End);
//   log_window4->moveCursor(QTextCursor::End);
//   log_window5->moveCursor(QTextCursor::End);
//   log_window6->moveCursor(QTextCursor::End);
  }

  return v;
 }

 virtual std::streamsize xsputn(const char *p, std::streamsize n)
 {
     QMutexLocker locker(&mutex2);
  m_string.append(p, p + n);

  int pos = 0;
  while (pos != std::string::npos)
  {
   pos = m_string.find('\n');
   if (pos != std::string::npos)
   {
    std::string tmp(m_string.begin(), m_string.begin() + pos);
    log_window1->append(tmp.c_str());
    log_window2->append(tmp.c_str());
    log_window3->append(tmp.c_str());
    log_window4->append(tmp.c_str());
    log_window5->append(tmp.c_str());
    log_window6->append(tmp.c_str());
    if(logfile.is_open())
        logfile << m_string.c_str() << std::endl;
    m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
   } else
   {
//       log_window1->moveCursor(QTextCursor::End);
//       log_window2->moveCursor(QTextCursor::End);
//       log_window3->moveCursor(QTextCursor::End);
//       log_window4->moveCursor(QTextCursor::End);
//       log_window5->moveCursor(QTextCursor::End);
//       log_window6->moveCursor(QTextCursor::End);
   }
  }

  return n;
 }

private:
 std::ostream &m_stream;
 std::streambuf *m_old_buf;
 std::string m_string;
 QTextEdit* log_window1;
 QTextEdit* log_window2;
 QTextEdit* log_window3;
 QTextEdit* log_window4;
 QTextEdit* log_window5;
 QTextEdit* log_window6;
 char* logfname;
 File logfile;
 mutable QMutex mutex1;
 mutable QMutex mutex2;
};

#endif
