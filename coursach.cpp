
#include <iostream>

#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "utils.h"

int main()
{
   std::ios::sync_with_stdio(false);
   size_t test_number;
   std::cout << "Menu: How much tests do you wont to run?\n";
   std::cin >> test_number;

   std::ifstream input;
   for (size_t i = 1; i <= test_number; i++)
   {
      std::string input_filename = "input/input" + std::to_string(i) + ".txt";
      input.open(input_filename);
      try
      {
         input.exceptions(std::ifstream::failbit);
         std::cout << "File input " + std::to_string(i) + " was opened\n";

         size_t n_in;
         double gamma_in;
         input >> n_in >> gamma_in;

         size_t n = n_in;
         double gamma = gamma_in;
         const size_t m = 2 * n + 1;

         std::vector<double> di(m);
         std::vector<size_t> ig(m + 1);
         std::vector<double> gg;
         std::vector<double> G(m);
         std::vector<double> F(9);

         std::vector<std::vector<double>> B(3, std::vector<double>(3));
         std::vector<std::vector<double>> C(3, std::vector<double>(3));

         std::vector<double> r;

         r.reserve(3 * n - 1);

         load_data_to(n, r, input, i);
         assembly(n, m, gamma, ig, di, gg, B, C, F, G, r, i);
         region(n, m, di, G, r, i);

         SLAU_LLT(n, m, ig, di, gg, G);

         std::ofstream output;
         std::string output_filename = "output/output" + std::to_string(i) + ".txt";
         output.open(output_filename);

         if (output.is_open())
         {
            std::cout << "File output " + std::to_string(i) + " was opened\n";
            export_data_to(n, m, G, r, output, i);
         }
         else
         {
            throw "File output " + std::to_string(i) + " was NOT opened\n";
         }

         input.close();
         output.close();
      }
      catch (const std::ios_base::failure& fail)
      {
         std::cout << fail.what() << '\n';
      }
   }
   return 0;
}
