import java.text.DecimalFormat;

public class MatrixCLocal
{
    double [][] matrixC_1pc;            // macierze przed lokalną agregacją do macierzy C, osobno dla każdego pc
    double [][] matrixC_2pc;
    double [][] matrixC_3pc;
    double [][] matrixC_4pc;

    double [][] matrixC;                // lokalna macierz [C]
    UniversalElement universalElement;
    Element element;

    public MatrixCLocal(Element element, Jakobian2D jakobian2D, Dane dane) {
        matrixC_1pc = new double[4][4];
        matrixC_2pc = new double[4][4];
        matrixC_3pc = new double[4][4];
        matrixC_4pc = new double[4][4];
        universalElement = new UniversalElement();
        this.element = element;

        // macierz c dla pierwszego punktu całkowania--------------------------------------------------------------C_1pc
        for (int i = 0; i < 4; i++) {
            matrixC_1pc[0][i] = universalElement.dN[0][i] * universalElement.dN[0][0] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[0];
            matrixC_1pc[1][i] = universalElement.dN[0][i] * universalElement.dN[1][0] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[0];
            matrixC_1pc[2][i] = universalElement.dN[0][i] * universalElement.dN[2][0] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[0];
            matrixC_1pc[3][i] = universalElement.dN[0][i] * universalElement.dN[3][0] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[0];
        }

        // macierz c dla drugiego punktu całkowania----------------------------------------------------------------C_2pc
        for (int i = 0; i<4; i++) {
            matrixC_2pc[0][i] = universalElement.dN[i][1] * universalElement.dN[0][1] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[1];
            matrixC_2pc[1][i] = universalElement.dN[i][1] * universalElement.dN[1][1] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[1];
            matrixC_2pc[2][i] = universalElement.dN[i][1] * universalElement.dN[2][1] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[1];
            matrixC_2pc[3][i] = universalElement.dN[i][1] * universalElement.dN[3][1] *
                               dane.specificHeat * dane.density * jakobian2D.detJ[1];
        }

        // macierz c dla trzeciego punktu całkowania---------------------------------------------------------------C_3pc
        for (int i = 0; i<4; i++) {
            matrixC_3pc[0][i] = universalElement.dN[i][2] * universalElement.dN[0][2] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[2];
            matrixC_3pc[1][i] = universalElement.dN[i][2] * universalElement.dN[1][2] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[2];
            matrixC_3pc[2][i] = universalElement.dN[i][2] * universalElement.dN[2][2] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[2];
            matrixC_3pc[3][i] = universalElement.dN[i][2] * universalElement.dN[3][2] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[2];
        }


        // macierz c dla czwartego punktu całkowania---------------------------------------------------------------C_4pc
        for (int i = 0; i<4; i++) {
            matrixC_4pc[0][i] = universalElement.dN[i][3] * universalElement.dN[0][3] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[3];
            matrixC_4pc[1][i] = universalElement.dN[i][3] * universalElement.dN[1][3] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[3];
            matrixC_4pc[2][i] = universalElement.dN[i][3] * universalElement.dN[2][3] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[3];
            matrixC_4pc[3][i] = universalElement.dN[i][3] * universalElement.dN[3][3] *
                                dane.specificHeat * dane.density * jakobian2D.detJ[3];
        }




        matrixC = new double[4][4]; //-----------------------------------------------------------------------lokalna [C]

        for (int i = 0; i < 4; i++) {
            matrixC[0][i] = matrixC_1pc[0][i] + matrixC_2pc[0][i] + matrixC_3pc[0][i] + matrixC_4pc[0][i];
            matrixC[1][i] = matrixC_1pc[1][i] + matrixC_2pc[1][i] + matrixC_3pc[1][i] + matrixC_4pc[1][i];
            matrixC[2][i] = matrixC_1pc[2][i] + matrixC_2pc[2][i] + matrixC_3pc[2][i] + matrixC_4pc[2][i];
            matrixC[3][i] = matrixC_1pc[3][i] + matrixC_2pc[3][i] + matrixC_3pc[3][i] + matrixC_4pc[3][i];
        }
    }










    DecimalFormat decimalFormat = new DecimalFormat("#0.000");
    void printMatrixCforPC() {
        System.out.println("\n\n----\t----\t----\t--\t MATRIX C LOKALNIE \t--\t----\t----\t----");
        System.out.println("Macierze składające się na macierz lokalną [C]: ");
        System.out.println("Dla 1 punktu całkowania: ");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(matrixC_1pc[i][0]) + "\t\t" +
                                decimalFormat.format(matrixC_1pc[i][1]) + "\t\t" +
                                decimalFormat.format(matrixC_1pc[i][2]) + "\t\t" +
                                decimalFormat.format(matrixC_1pc[i][3]) + "\t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 2 punktu całkowania: ");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(matrixC_2pc[i][0]) + "\t\t" +
                                decimalFormat.format(matrixC_2pc[i][1]) + "\t\t" +
                                decimalFormat.format(matrixC_2pc[i][2]) + "\t\t" +
                                decimalFormat.format(matrixC_2pc[i][3]) + "\t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 3 punktu całkowania: ");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(matrixC_3pc[i][0]) + "\t\t" +
                                decimalFormat.format(matrixC_3pc[i][1]) + "\t\t" +
                                decimalFormat.format(matrixC_3pc[i][2]) + "\t\t" +
                                decimalFormat.format(matrixC_3pc[i][3]) + "\t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 4 punktu całkowania: ");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(matrixC_4pc[i][0]) + "\t\t" +
                                decimalFormat.format(matrixC_4pc[i][1]) + "\t\t" +
                                decimalFormat.format(matrixC_4pc[i][2]) + "\t\t" +
                                decimalFormat.format(matrixC_4pc[i][3]) + "\t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        }

        void printLocalMatrixC () {
            System.out.println("\n\n=\t==\t==\t==\t==\t==\t==\tMACIERZ C LOKALNIE\t==\t==\t==\t==\t==\t==\t=");
            for (int i = 0; i < 4; i++) {
                System.out.println( decimalFormat.format(matrixC[i][0]) + "\t\t" +
                                    decimalFormat.format(matrixC[i][1]) + "\t\t" +
                                    decimalFormat.format(matrixC[i][2]) + "\t\t" +
                                    decimalFormat.format(matrixC[i][3]) + "\t\t"); }
            System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        }
}
