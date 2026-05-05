#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return exp(cos(x));
}

double rect(double a, double b, int n) {
    double h = (b - a) / n;
    double S = 0.0;
    for (int i = 0; i < n; i++) {
        double x_mid = a + (i + 0.5) * h;
        S += f(x_mid);
    }
    return S * h;
}

double trap(double a, double b, int n) {
    double h = (b - a) / n;
    double S = (f(a) + f(b)) / 2.0;
    for (int i = 1; i < n; i++) {
        S += f(a + i * h);
    }
    return S * h;
}

double simpson(double a, double b, int n) {
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += (i % 2 == 0) ? 2 * f(x) : 4 * f(x);
    }
    return (h / 3.0) * sum;
}

double runge_error_rect(double a, double b, int n) {
    double I_n = rect(a, b, n);
    double I_2n = rect(a, b, 2 * n);
    return fabs(I_2n - I_n) / 3.0;
}

double runge_error_trap(double a, double b, int n) {
    double I_n = trap(a, b, n);
    double I_2n = trap(a, b, 2 * n);
    return fabs(I_2n - I_n) / 3.0;
}

double runge_error_simpson(double a, double b, int n) {
    double I_n = simpson(a, b, n);
    double I_2n = simpson(a, b, 2 * n);
    return fabs(I_2n - I_n) / 15.0;
}

int find_n_rect(double a, double b, double eps) {
    int n = 2;
    while (runge_error_rect(a, b, n) > eps) {
        n *= 2;
    }
    return n;
}

int find_n_trap(double a, double b, double eps) {
    int n = 2;
    while (runge_error_trap(a, b, n) > eps) {
        n *= 2;
    }
    return n;
}

int find_n_simpson(double a, double b, double eps) {
    int n = 2;
    while (runge_error_simpson(a, b, n) > eps) {
        n *= 2;
    }
    return n;
}

int main() {
    double a = 0.0, b = 1.0;
    double eps;

    cout << "вычисление интеграла от exp(cos(x)) на [0, 1]" << endl;
    cout << "введите точность: ";
    cin >> eps;

    cout << endl;

    int n_rect = find_n_rect(a, b, eps);
    int n_trap = find_n_trap(a, b, eps);
    int n_simp = find_n_simpson(a, b, eps);

    double rec = rect(a, b, n_rect);
    double trapp = trap(a, b, n_trap);
    double simp = simpson(a, b, n_simp);

    cout << "метод прямоугольников:" << endl;
    cout << "| требуется отрезков: " << n_rect << endl;
    cout << "| значение: " << rec << endl;
    cout << "| оценка погрешности: " << runge_error_rect(a, b, n_rect) << endl;
    cout << endl;

    cout << "метод трапеций:" << endl;
    cout << "| требуется отрезков: " << n_trap << endl;
    cout << "| значение: " << trapp << endl;
    cout << "| оценка погрешности: " << runge_error_trap(a, b, n_trap) << endl;
    cout << endl;

    cout << "метод Симпсона:" << endl;
    cout << "| требуется отрезков: " << n_simp << endl;
    cout << "| значение: " << simp << endl;
    cout << "| оценка погрешности: " << runge_error_simpson(a, b, n_simp) << endl;

    return 0;
}
