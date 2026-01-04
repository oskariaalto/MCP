#ifndef PROGRESS_BAR_HPP
#define PROGRESS_BAR_HPP

#include <string>
#include <iostream>
#include <sys/ioctl.h>
#include <unistd.h>

class ProgressBar
{
private:
  std::string name;     
  int total_iterations; 
  int bar_width;        
  int num_iterations = 0;

  
  int get_terminal_width() const;

  
  void draw_progress_bar() const;

public:
  
  ProgressBar(const std::string &name, int total_iterations);

  
  void update();
};

#endif 
