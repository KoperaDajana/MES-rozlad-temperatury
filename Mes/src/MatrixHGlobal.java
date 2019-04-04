import java.text.DecimalFormat;

public class MatrixHGlobal
{
    double [][] globalMatrixH;
    Dane dane;

    public MatrixHGlobal(Dane dane)
    {
        globalMatrixH = new double[dane.nh][dane.nh];           // wielkość: ilość_węzłów x ilość_węzłów
        this.dane = dane;

        for (int i = 0; i<dane.nh; i++) {                       // wypełnienie zerami -> uniknięcie zaśmiecania
            for (int j = 0; j<dane.nh; j++) {
                globalMatrixH[i][j] = 0;
            }
        }
    }


    // agregacja [H]
    void ageracjaLokalnychH(MatrixHLocal matrixHLocal) {
        for (int i = 0; i<4; i++) {                             // 4 bo wymiary macierzy lokalnych są właśnie 4x4
            for (int j = 0; j < 4; j++) {
                globalMatrixH[matrixHLocal.element.node[i].ID][matrixHLocal.element.node[j].ID] += matrixHLocal.matrixH[i][j];
            }
        }
    }

    // agregacja [H_BC]
    void ageracjaLokalnychHBC(MatrixHLocal matrixHLocal) {
        for (int i = 0; i<4; i++) {                             // 4 bo wymiary macierzy lokalnych są właśnie 4x4
            for (int j = 0; j < 4; j++) {
                globalMatrixH[matrixHLocal.element.node[i].ID][matrixHLocal.element.node[j].ID] += matrixHLocal.matrixHBC[i][j];
            }
        }
    }



    // agregacja [H] = [H] + [C]/dT (gdzie dT to krok czasowy symulacji = 50)
    void ageracjaHplusCprzezKrokCzasowy(MatrixCGlobal matrixCGlobal) {
        for (int i = 0; i<dane.nh; i++) {                       //4 bo wymiary macierzy lokalnych są właśnie 4x4
            for (int j = 0; j < dane.nh; j++) {
                globalMatrixH[i][j] += matrixCGlobal.globalMatrixC[i][j]/dane.simulationStepTime;
            }
        }
    }



    //------------------------------------------------------------------------------------------------------------------
    DecimalFormat decimalFormat = new DecimalFormat("#0.000");
    void printMartixHGlobal() {
        System.out.println("==\t==\t==\t==\t==\t==\t==\t==\t==\tMACIERZ [H] " +
                "GLOBALNIE\t==\t==\t==\t==\t==\t==\t==\t==\t==\t");
        for (int i = 0; i<dane.nh; i++) {
            for (int j = 0; j<dane.nh; j++) {
                System.out.print(decimalFormat.format(globalMatrixH[i][j]) + "   \t");
            }
            System.out.println();
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
    }

}
