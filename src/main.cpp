#include <iostream>
#include "Exercise_1/exercise1.hpp"
#include "Exercise_2/exercise2.hpp"
#include "Exercise_3/exercise3.hpp"
#include "Exercise_4/exercise4.hpp"
#include "Exercise_5/exercise5.hpp"

int main()
{
  int choice;

  
  std::cout << "Choose a task to run:" << std::endl;
  std::cout << "1. Run Exercise 1" << std::endl;
  std::cout << "2. Run Exercise 2" << std::endl;
  std::cout << "3. Run Exercise 3" << std::endl;
  std::cout << "4. Run Exercise 4" << std::endl;
  std::cout << "5. Run Exercise 5" << std::endl;
  std::cout << "Enter your choice: ";
  std::cin >> choice;

  
  switch (choice)
  {
  case 1:
    runExercise1();
    break;
  case 2:
    runExercise2();
    break;
  case 3:
    runExercise3();
    break;
  case 4:
    runExercise4();
    break;
  case 5:
    runExercise5();
    break;
  default:
    std::cout << "Invalid choice." << std::endl;
    break;
  }

  return 0;
}
