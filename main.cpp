#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

void mostrarSistema(const vector<vector<double>>& A, const vector<double>& b) {
    cout << "Sistema de ecuaciones:" << endl;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            cout << A[i][j] << "*x" << j + 1;
            if (j != A[i].size() - 1) cout << " + ";
        }
        cout << " = " << b[i] << endl;
    }
    cout << endl;
}

// Función para resolver un sistema de ecuaciones lineales usando el método de Jacobi
vector<double> jacobi(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x, double tolerancia, int max_iteraciones) {
    int n = A.size();
    vector<double> x_nuevo(n);
    vector<double> errores(n);
    
    cout << fixed << setprecision(6);
    
    for (int iter = 0; iter < max_iteraciones; ++iter) {
        cout << "Iteracion " << iter + 1 << ":" << endl;
        bool cumple_tolerancia = true;

        for (int i = 0; i < n; ++i) {
            x_nuevo[i] = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    x_nuevo[i] -= A[i][j] * x[j];
                }
            }
            x_nuevo[i] /= A[i][i];

            // Calcular el error relativo si el nuevo valor es distinto de 0
            errores[i] = fabs((x_nuevo[i] - x[i]) / (x_nuevo[i] != 0 ? x_nuevo[i] : 1)) * 100;
            if (errores[i] > tolerancia * 100) cumple_tolerancia = false;
        }

        // Mostrar valores aproximados y errores
        for (int i = 0; i < n; ++i) {
            cout << "x[" << i + 1 << "] = " << x_nuevo[i] << " \tError = " << errores[i] << "%" << endl;
        }
        cout << endl;

        if (cumple_tolerancia) {
            cout << "Convergencia alcanzada en " << iter + 1 << " iteraciones.\n" << endl;
            return x_nuevo;
        }

        x = x_nuevo;
    }

    cout << " No se alcanzo la convergencia en " << max_iteraciones << " iteraciones.\n" << endl;
    return x;
}

int main() {
    // Define la matriz de coeficientes A
    vector<vector<double>> A = {{10, 2, 1}, {1, 5, 1}, {2, 3, 8}};

    // Define el vector de términos independientes b
    vector<double> b = {7, -8, 15};

    // Valores iniciales para la solución x
    vector<double> x_inicial = {0, 0, 0};

    // Tolerancia y número máximo de iteraciones
    double tolerancia = 0.0001;
    int max_iteraciones = 100;

    // Mostrar sistema y tolerancia
    mostrarSistema(A, b);
    cout << "Tolerancia de error relativa maxima permitida: " << tolerancia * 100 << "%\n" << endl;

    // Llamar a la función jacobi
    vector<double> solucion = jacobi(A, b, x_inicial, tolerancia, max_iteraciones);

    // Imprimir la solución final
    cout << "Solucion final aproximada del sistema de ecuaciones:" << endl;
    for (int i = 0; i < solucion.size(); ++i) {
        cout << "x[" << i + 1 << "] = " << solucion[i] << endl;
    }

    return 0;
}
