
#include <unistd.h>
#include <cstdlib>
#include <sys/time.h>
#include <list>

class split {

public:
  std::string label;
  struct timeval time;

  split(std::string l);
};

class timer {
private:
  std::string title;
  std::list<split> splits;
  float timesub(struct timeval t1, struct timeval t2);
public:
  timer(std::string title);
  void addsplit(std::string l);
  void print();
};
