import java.text.DecimalFormat;

// agregacja lokalnych macierzy [C] w dużą lokalną macierz
public class MatrixCGlobal
{
    double [][] globalMatrixC;
    Dane dane;

    public MatrixCGlobal(Dane dane)
    {
        globalMatrixC = new double[dane.nh][dane.nh];       // macierz o wymiarach: liczba węzłów x liczba węzłów
        this.dane = dane;

        for (int i = 0; i<dane.nh; i++) {                   // wyzerowanie, żeby pozbyć się zaśmiecenia
            for (int j = 0; j<dane.nh; j++) {
                globalMatrixC[i][j] = 0;
            }
        }
    }


    // w celu połączenia wszystkich lokalnych macierzy C w jedną globalną -------------------------------> AGREGACJA [C]
    void ageracjaLokalnychC(MatrixCLocal matrixCLocal) {                // przesyłanie jako argument macierzy C lokalnej
        for (int i = 0; i<4; i++) {                                     // 4 bo wymiary macierzy lokalnych są właśnie 4x4
            for (int j = 0; j<4; j++) {
                globalMatrixC[matrixCLocal.element.node[i].ID][matrixCLocal.element.node[j].ID] += matrixCLocal.matrixC[i][j];
                // ID pomaga w odpowiednim ustawianiu elementów w macierzy globalnej
            }
        }
    }


    // --------------------------------------------------------------------------------------------------------wypisanie
    DecimalFormat decimalFormat = new DecimalFormat("#0.000");
    void printMartixCGlobal() {
        System.out.println("\n\n==\t==\t==\t==\t==\t==\t==\t==\t==\tMACIERZ [C] " +
                "GLOBALNIE\t==\t==\t==\t==\t==\t==\t==\t==\t==\t");
        for (int i = 0; i<dane.nh; i++) {
            for (int j = 0; j<dane.nh; j++) {
                System.out.print(decimalFormat.format(globalMatrixC[i][j]) + "   \t"); }
            System.out.println();
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
    }

}
