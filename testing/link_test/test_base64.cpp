
#include <base64.h>
#include <iostream>
#include <string>

int main()
{
  std::string ip = "base64";
  std::string op;

  bool result = Base64::Encode(ip, &op);
  std::cout << "Result = " << result << std::endl;

  return EXIT_SUCCESS;
}