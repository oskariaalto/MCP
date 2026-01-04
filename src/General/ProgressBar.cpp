#include "ProgressBar.hpp"


ProgressBar::ProgressBar(const std::string &name, int total_iterations)
    : name(name), total_iterations(total_iterations)
{
  bar_width = get_terminal_width() - 30; 
  if (bar_width < 10)
  {
    bar_width = 10; 
  }
}


int ProgressBar::get_terminal_width() const
{
  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  return w.ws_col;
}


void ProgressBar::draw_progress_bar() const
{
  double progress = (double)num_iterations / total_iterations; 
  int pos = bar_width * progress;                              

  std::cout << "\r" << name << " ["; 
  for (int i = 0; i < bar_width; ++i)
  {
    if (i < pos)
    {
      std::cout << "#"; 
    }
    else
    {
      std::cout << " "; 
    }
  }
  std::cout << "] " << int(progress * 100.0) << " %"; 
  std::cout.flush();                                  
}


void ProgressBar::update()
{
  num_iterations++;
  draw_progress_bar();
  if (num_iterations == total_iterations)
  {
    std::cout << "\n"; 
  }
}
