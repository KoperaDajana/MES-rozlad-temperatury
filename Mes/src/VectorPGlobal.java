import java.text.DecimalFormat;

public class VectorPGlobal {
    double[] vectorP;
    Dane dane;

    public VectorPGlobal(Dane dane) {
        this.dane = dane;
        vectorP = new double[dane.nh]; //ilość węzłów - nh

        for (int i = 0; i < dane.nh; i++) { vectorP[i] = 0; }                   // wyzerowanie wektora {P}
    }


    //zerowanie globalnego wektora {P}
    void zerowaniePGlobalnego() {
        for (int i = 0; i < dane.nh; i++) { vectorP[i] = 0; }
    }

    // funkcja, która agreguje lokalny wektor {P} do globalnego wektora {P}
    void agregacjaVektoraPLokalnegoDoGlobalnego(MatrixHLocal matrixHLocal) {
        for (int i = 0; i < 4; i++) {
            vectorP[matrixHLocal.element.node[i].ID] += matrixHLocal.vectorPLocal[i];
        }
    }

    // agregacja wektora {P} z uwzględnionym krokiem czasowym i temperaturą początkową --------> ({P} + {[C]/dT} * {T0})
    void agregacjaWekoraPFinal(MatrixCGlobal matrixCGlobal) {
        for (int i = 0; i < dane.nh; i++) {
            for (int j = 0; j < dane.nh; j++) {
                vectorP[i] += (matrixCGlobal.globalMatrixC[i][j] / dane.simulationStepTime) * dane.initialTemp;
            }
        }
    }



    // agregacja [C] globalnego/dT do Gaussa P globalnego {P} = {P} + {[C] / dT} * {Tn}
    // Tn zmienia się z każdym krokiem symulacji ---> gdzie Tn są to temperatury pobierane z vectorTempInNodes (Proces)
    void metodaGaussa(MatrixCGlobal matrixCGlobal, double[] temperatureOfNodes) {
        for (int i = 0; i<dane.nh; i++) {
            for (int j = 0; j<dane.nh; j++) {
                vectorP[i] += (matrixCGlobal.globalMatrixC[i][j] / dane.simulationStepTime) * temperatureOfNodes[j];
            }
        }
    }


    DecimalFormat decimalFormat = new DecimalFormat("#0.0");
    void printVectorP() {
        System.out.println("\n\n==\t==\t==\t==\t==\t==\t==\t==\t==\tWEKTOR {P} " +
                "GLOBALNIE\t==\t==\t==\t==\t==\t==\t==\t==\t==\t");
        for (int i = 0; i < dane.nh; i++) {
            System.out.print(decimalFormat.format(vectorP[i]) + "\t");
        }
        System.out.println("\n--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
    }

}
