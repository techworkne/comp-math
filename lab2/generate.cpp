#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>

using namespace std;

void generateRandomDiagonallyDominant(int n, vector<vector<double>>& A, vector<double>& b) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(-50.0, 50.0);

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                A[i][j] = dist(gen);
                sum += std::abs(A[i][j]);
            }
        }
        // Обеспечиваем диагональное преобладание
        A[i][i] = sum + std::abs(dist(gen)) + 1.0; 
        if (dist(gen) < 0) A[i][i] = -A[i][i];
        
        b[i] = dist(gen);
    }
}

void generateHilbert(int n, vector<vector<double>>& A, vector<double>& b) {
    for (int i = 0; i < n; ++i) {
        b[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            A[i][j] = 1.0 / (i + j + 1.0);
            b[i] += A[i][j]; // Если сумма строк, то точное решение x = [1, 1, ..., 1]
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Использование: " << argv[0] << " <размер_n> <тип: 1-Случайная, 2-Гильберт> <имя_файла>\n";
        return 1;
    }

    int n = stoi(argv[1]);
    int type = stoi(argv[2]);
    string filename = argv[3];

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    if (type == 1) {
        generateRandomDiagonallyDominant(n, A, b);
        cout << "Сгенерирована случайная матрица с диагональным преобладанием.\n";
    } else if (type == 2) {
        generateHilbert(n, A, b);
        cout << "Сгенерирована матрица Гильберта.\n";
    } else {
        cout << "Неизвестный тип матрицы!\n";
        return 1;
    }

    ofstream out(filename);
    if (!out.is_open()) {
        cerr << "Ошибка создания файла!\n";
        return 1;
    }

    out << n << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            out << fixed << setprecision(10) << A[i][j] << " ";
        }
        out << b[i] << "\n";
    }
    
    out.close();
    cout << "Данные записаны в " << filename << "\n";
    return 0;
}