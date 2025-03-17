#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./program_name <csv_file_path>" << std::endl;
        return 1;
    }

    std::string filePath = argv[1];
    std::ifstream file(filePath);

    double weights[4] = { 0.11, 0.04, 0.0001, 0.0005 };
    double predicition = 0.;
    const int GOP_SIZE = 4;

    if (file.is_open()) {
        std::string line;
        int lineCount = 0;
        while (getline(file, line)) {
            std::vector<double> numbers;
            std::stringstream ss(line);
            std::string cell;

            if (lineCount++ == 0)
                continue;

            // 读取每一行的前三列数字
            for (int i = 0; i < 3; i++) {
                if (getline(ss, cell, ',')) {
                    double number = std::stod(cell);
                    numbers.push_back(number);
                }
            }

            // 输出前三列数字
            int poc = static_cast<int>(numbers[0]) % 250;
            double E = numbers[1], h = numbers[2];

            if (poc == 0)
                predicition += E * weights[0];  
            else if(poc % GOP_SIZE == 0)
            {
                predicition += h * weights[1];
            }
            else if (poc % GOP_SIZE == 2)
            {
                predicition += h * weights[2];
            }
            else
            {
                predicition += h * weights[3];
            }


            //std::cout << poc << ", " << E << ", " << h << std::endl;

        }
        predicition /= lineCount;
        std::cout << predicition << std::endl;
        file.close();
    }

    else {
        std::cout << "Unable to open file: " << filePath << std::endl;
    }

    return 0;
}