
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <list>

#include "timer.h"

split::split(std::string l) {
  label=l;
  gettimeofday(&time, NULL);
}

float timer::timesub(struct timeval t1, struct timeval t2) {
  float sec, usec;
  if (t2.tv_usec > t1.tv_usec) {
    sec = t2.tv_sec-t1.tv_sec;
    usec = t2.tv_usec-t1.tv_usec;
  } else {
    sec = t2.tv_sec-t1.tv_sec-1;
    usec = (1000000+t2.tv_usec-t1.tv_usec);
  }
  return (float) sec + ((float) usec / 1000000.0);
}

timer::timer(std::string t) {
  title=t;
  timer::addsplit("initialized");
}

void timer::addsplit(std::string l) {
  split *s=new split(l);
  splits.push_back(*s);
}
  
void timer::print() {
  std::list<split>::iterator s;
  float delta, fromstart;
  s=splits.begin();
  split &first = *s;
  split *last = &(*s);
  s++;
  std::cout << title << "\n";
  std::cout << std::setw(10) << "split" << std::setw(10) << "last" << std::setw(10) << "total\n";
  for (; s!=splits.end(); ++s) {
    fromstart = timesub(first.time, s->time);
    delta = timesub(last->time, s->time);
    last=&(*s);
    std::cout << std::setw(10) << s->label << "  " <<
      std::setw(10) << std::setprecision(6) << delta <<
      std::setw(10) << std::setprecision(6) << fromstart <<
      "\n";
  }
}

