#ifndef INPUT_READER_HPP
#define INPUT_READER_HPP

#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include <stdexcept>

class InputReader
{
public:
  
  InputReader(const std::string &filename) : filename_(filename) {}

  
  bool readInputs(std::unordered_map<std::string, double> &values)
  {
    std::ifstream input_file(filename_);
    if (!input_file)
    {
      std::cerr << "Error: Could not open file " << filename_ << std::endl;
      return false;
    }

    std::string name;
    std::string value;

    while (input_file >> name >> value)
    {
      try
      {
        values[name] = std::stod(value);
      }
      catch (const std::invalid_argument &e)
      {
        std::cerr << "Error: Invalid value format for " << name << " in file." << std::endl;
        return false;
      }
    }

    input_file.close();
    return true;
  }

private:
  std::string filename_;
};

#endif 
