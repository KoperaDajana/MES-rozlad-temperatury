import java.text.DecimalFormat;

// klasa zawierająca jakobian oraz jego wyznacznik w celu przejścia z układu globalnego na lokalny
// arkusze excel: Jakobian2D oraz MatrixH
public class Jakobian2D
{
    UniversalElement universalElement = new UniversalElement();
    // macierz zawierająca pochodne dx/dEta oraz dy/dKsi
    double [][] Jakobian = new double[4][4];            // tabela kolorowa 1 --> Jakobian

    // wektor wyznaczników dla poszczególnych punktów całkowania
    double [] detJ = new double[4];                     // tabela kolorowa 2 --> det J (wektor)
    // macierz przekształceń
    double [][] jakobianPrzezDet = new double[4][4];    // tabela kolorowa 3 --> Jakobian/det J

    // dla przejścia na globalny układ
    double [][] dNdx = new double[4][4];                // excel MatrixH
    double [][] dNdy = new double[4][4];

    public Jakobian2D(Element element)
    {
        // ------------------------------------------------------------------------------------------JAKOBIAN (bez detJ)
        for (int j = 0; j<4; j++) {
            Jakobian[0][j] = universalElement.dNdKsi[0][j] * element.node[0].x +
                             universalElement.dNdKsi[1][j] * element.node[1].x +
                             universalElement.dNdKsi[2][j] * element.node[2].x +
                             universalElement.dNdKsi[3][j] * element.node[3].x;

            Jakobian[1][j] = universalElement.dNdKsi[0][j] * element.node[0].y +
                             universalElement.dNdKsi[1][j] * element.node[1].y +
                             universalElement.dNdKsi[2][j] * element.node[2].y +
                             universalElement.dNdKsi[3][j] * element.node[3].y;

            Jakobian[2][j] = universalElement.dNdEta[0][j] * element.node[0].x +
                             universalElement.dNdEta[1][j] * element.node[1].x +
                             universalElement.dNdEta[2][j] * element.node[2].x +
                             universalElement.dNdEta[3][j] * element.node[3].x;

            Jakobian[3][j] = universalElement.dNdEta[0][j] * element.node[0].y +
                             universalElement.dNdEta[1][j] * element.node[1].y +
                             universalElement.dNdEta[2][j] * element.node[2].y +
                             universalElement.dNdEta[3][j] * element.node[3].y;
        }


        // ---------------------------------------------------------------------------------------------------------detJ
        for (int i = 0; i<4; i++) {
            detJ[i] = Jakobian[0][i] * Jakobian[3][i] - Jakobian[1][i] * Jakobian [2][i];
        }

        // -------------------------------------------------------------------------------------------------JAKOBIAN/DET
        for (int i = 0; i<4; i++) {
            jakobianPrzezDet[0][i] = Jakobian[3][i]/detJ[i];
            jakobianPrzezDet[1][i] = Jakobian[1][i]/detJ[i];
            jakobianPrzezDet[2][i] = Jakobian[2][i]/detJ[i];
            jakobianPrzezDet[3][i] = Jakobian[0][i]/detJ[i];
        }


        // przejście z układu lokalnego dzieki stworzeniu Jakobianu, na układ globalny --->
        // stworzenie kolejnych dwóch macierzy, które będą odnosić się do układu globalnego ---> x, y

        // jakobian * dN/dKsi + jakobian * dN/dEta
        // dN / dx-------------------------------------------------------------------------------------------------dN/dx
        for (int i = 0; i<4; i++) {
            dNdx[0][i] = jakobianPrzezDet[0][0] * universalElement.dNdKsi[i][0]
                       + jakobianPrzezDet[1][0] * universalElement.dNdEta[i][0];
            dNdx[1][i] = jakobianPrzezDet[0][1] * universalElement.dNdKsi[i][1]
                       + jakobianPrzezDet[1][1] * universalElement.dNdEta[i][1];
            dNdx[2][i] = jakobianPrzezDet[0][2] * universalElement.dNdKsi[i][2]
                       + jakobianPrzezDet[1][2] * universalElement.dNdEta[i][2];
            dNdx[3][i] = jakobianPrzezDet[0][3] * universalElement.dNdKsi[i][3]
                       + jakobianPrzezDet[1][3] * universalElement.dNdEta[i][3];
        }


        // dN / dy-------------------------------------------------------------------------------------------------dN/dy
        for (int i = 0; i<4; i++) {
            dNdy[0][i] = jakobianPrzezDet[2][0] * universalElement.dNdKsi[i][0]
                       + jakobianPrzezDet[3][0] * universalElement.dNdEta[i][0];
            dNdy[1][i] = jakobianPrzezDet[2][1] * universalElement.dNdKsi[i][1]
                       + jakobianPrzezDet[3][1] * universalElement.dNdEta[i][1];
            dNdy[2][i] = jakobianPrzezDet[2][2] * universalElement.dNdKsi[i][2]
                       + jakobianPrzezDet[3][2] * universalElement.dNdEta[i][2];
            dNdy[3][i] = jakobianPrzezDet[2][3] * universalElement.dNdKsi[i][3]
                       + jakobianPrzezDet[3][3] * universalElement.dNdEta[i][3];
        }
        // dN/dx oraz dN/dy będą mogły być użyte do macierzy [H] dla 2D
    }








    // ------------------------------------------------------------------------------------------------------wypisywanie
    DecimalFormat decimalFormat = new DecimalFormat("#0.00000");
    void printJakobian() {
        // ------------------------------------------------------------------------------------------wypisanie Jakobianu
        System.out.println("\n\n----\t----\t----\t----\t JAKOBIAN 2D \t----\t----\t----\t----");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(Jakobian[i][0]) + " \t\t" +
                                decimalFormat.format(Jakobian[i][1]) + " \t\t" +
                                decimalFormat.format(Jakobian[i][2]) + " \t\t" +
                                decimalFormat.format(Jakobian[i][3]) + " \t\t"); }
         System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        // -----------------------------------------------------------------------------------------------wypisanie detJ
        System.out.println("\n\n--\t--\t--\t--\t--\t--\t--\t--\t Wyznacznik J\t--\t--\t--\t--\t--\t--\t");
        for (int i = 0; i<4; i++) { System.out.println("detJ" + (i+1) + ": " + decimalFormat.format(detJ[i])); }


        // --------------------------------------------------------------------------------------wypisanie Jakobian/detJ
        System.out.println("\n\n--\t--\t--\t--\t--\t--\t--\t--\t Jakobian/det J\t--\t--\t--\t--\t--\t--\t");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(jakobianPrzezDet[i][0]) + " \t\t" +
                                decimalFormat.format(jakobianPrzezDet[i][1]) + " \t\t" +
                                decimalFormat.format(jakobianPrzezDet[i][2]) + " \t\t" +
                                decimalFormat.format(jakobianPrzezDet[i][3]) + " \t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
    }



    void printJakobianGlobal() {
        // ------------------------------------------------------------------------------------wypisywanie dN/dx & dN/dy
        System.out.println("\n\n--\t--\t--\t--\t--\t--\t PRZEJŚCIE NA GLOBALNY \t--\t--\t--\t--\t--\t--\t");
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t dN/dx\t--\t--\t--\t--\t--\t--\t--\t--\t");
        System.out.println("dN1/dx" + "\t\t\t" + "dN2/dx" + "\t\t\t" + "dN3/dx" + "\t\t\t" + "dN4/dx");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(dNdx[i][0]) + " \t\t" +
                                decimalFormat.format(dNdx[i][1]) + " \t\t" +
                                decimalFormat.format(dNdx[i][2]) + " \t\t" +
                                decimalFormat.format(dNdx[i][3]) + " \t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("\n--\t--\t--\t--\t--\t--\t--\t--\t dN/dy\t--\t--\t--\t--\t--\t--\t--\t--\t");
        System.out.println("dN1/dy" + "\t\t\t" + "dN2/dy" + "\t\t\t" + "dN3/dy" + "\t\t\t" + "dN4/dy");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(dNdy[i][0]) + " \t\t" +
                                decimalFormat.format(dNdy[i][1]) + " \t\t" +
                                decimalFormat.format(dNdy[i][2]) + " \t\t" +
                                decimalFormat.format(dNdy[i][3]) + " \t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
    }
}
