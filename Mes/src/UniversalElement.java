import java.text.DecimalFormat;
import static java.lang.Math.sqrt;

// klasa reprezentująca uniwersalny element, który może być wykorzystywany do obliczania jakobianu
// zawiera 3 macierze (funkcje kształtu, funkcje kształtu po dKsi oraz po dEta)
public class UniversalElement
{
    // wypełnienie dwóch tablic dla ksi oraz eta wartościami
    double[] ksi = new double [] {((-1)/sqrt(3)), (1/sqrt(3)), (1/sqrt(3)), ((-1)/sqrt(3))};
    double[] eta = new double [] {((-1)/sqrt(3)), ((-1)/sqrt(3)), (1/sqrt(3)), (1/sqrt(3))};

    // macierze
    double [][] dN = new double [4][4];         // funkcje kształtu (N = 0.25 * (1 +/- E)*(1 +/- n))
    double [][] dNdKsi = new double [4][4];     // funkcje kształtu po dKsi (N = +/- 0.25 * (1 +/- n))
    double [][] dNdEta = new double [4][4];     // funkcje kształtu po dEta (N = +/- 0.25 * (1 +/- E))

    public UniversalElement() {
        // MACIERZ 1 -> Wartości funkcji kształtu---------------------------------------------------------------------dN
        // dla N1 (pierwsza kolumna)
        dN[0][0] = 0.25 * ((1 - ksi[0]) * (1 - eta[0]));    // N1 dla węzła 1
        dN[1][0] = 0.25 * ((1 - ksi[1]) * (1 - eta[1]));    // N1 dla węzła 2
        dN[2][0] = 0.25 * ((1 - ksi[2]) * (1 - eta[2]));    // N1 dla węzła 3
        dN[3][0] = 0.25 * ((1 - ksi[3]) * (1 - eta[3]));    // N1 dla węzła 4

        // dla N2 (druga kolumna)
        dN[0][1] = 0.25 * ((1 + ksi[0]) * (1 - eta[0]));    // N2 dla węzła 1
        dN[1][1] = 0.25 * ((1 + ksi[1]) * (1 - eta[1]));    // N2 dla węzła 2
        dN[2][1] = 0.25 * ((1 + ksi[2]) * (1 - eta[2]));    // N2 dla węzła 3
        dN[3][1] = 0.25 * ((1 + ksi[3]) * (1 - eta[3]));    // N2 dla węzła 4

        // dla N3 (trzecia kolumna)
        dN[0][2] = 0.25 * ((1 + ksi[0]) * (1 + eta[0]));    // N3 dla węzła 1
        dN[1][2] = 0.25 * ((1 + ksi[1]) * (1 + eta[1]));    // N3 dla węzła 2
        dN[2][2] = 0.25 * ((1 + ksi[2]) * (1 + eta[2]));    // N3 dla węzła 3
        dN[3][2] = 0.25 * ((1 + ksi[3]) * (1 + eta[3]));    // N3 dla węzła 4

        // dla N4 (czwarta kolumna)
        dN[0][3] = 0.25 * ((1 - ksi[0]) * (1 + eta[0]));    // N4 dla węzła 1
        dN[1][3] = 0.25 * ((1 - ksi[1]) * (1 + eta[1]));    // N4 dla węzła 2
        dN[2][3] = 0.25 * ((1 - ksi[2]) * (1 + eta[2]));    // N4 dla węzła 3
        dN[3][3] = 0.25 * ((1 - ksi[3]) * (1 + eta[3]));    // N4 dla węzła 4


        // MACIERZ 2 -> Funkcje kształtu po KSI-------------------------------------------------------------------dNdKSI
        // wypełniane "kolumnowo"
        dNdKsi[0][0] = -0.25 * (1 - eta[0]);                // N1 dla eta[0] = -0,577350269
        dNdKsi[1][0] =  0.25 * (1 - eta[0]);                // N2 dla eta[0]
        dNdKsi[2][0] =  0.25 * (1 + eta[0]);                // N3 dla eta[0]
        dNdKsi[3][0] = -0.25 * (1 + eta[0]);                // N4 dla eta[0]


        dNdKsi[0][1] = -0.25 * (1 - eta[1]);                // N1 dla eta[1] = -0,577350269
        dNdKsi[1][1] =  0.25 * (1 - eta[1]);                // N2 dla eta[1]
        dNdKsi[2][1] =  0.25 * (1 + eta[1]);                // N3 dla eta[1]
        dNdKsi[3][1] = -0.25 * (1 + eta[1]);                // N4 dla eta[1]


        dNdKsi[0][2] = -0.25 * (1 - eta[2]);                // N1 dla eta[2] = 0,577350269
        dNdKsi[1][2] =  0.25 * (1 - eta[2]);                // N2 dla eta[2]
        dNdKsi[2][2] =  0.25 * (1 + eta[2]);                // N3 dla eta[2]
        dNdKsi[3][2] = -0.25 * (1 + eta[2]);                // N4 dla eta[2]


        dNdKsi[0][3] = -0.25 * (1 - eta[3]);                // N1 dla eta[3] = 0,577350269
        dNdKsi[1][3] =  0.25 * (1 - eta[3]);                // N2 dla eta[3]
        dNdKsi[2][3] =  0.25 * (1 + eta[3]);                // N3 dla eta[3]
        dNdKsi[3][3] = -0.25 * (1 + eta[3]);                // N4 dla eta[3]


        // MACIERZ 3 -> Funkcje kształtu po ETA-------------------------------------------------------------------dNdETA
        // wypełniane "kolumnowo"
        dNdEta[0][0] = -0.25 * (1 - ksi[0]);                // N1 dla ksi[0] = -0,577350269
        dNdEta[1][0] = -0.25 * (1 + ksi[0]);                // N2 dla ksi[0]
        dNdEta[2][0] =  0.25 * (1 + ksi[0]);                // N3 dla ksi[0]
        dNdEta[3][0] =  0.25 * (1 - ksi[0]);                // N4 dla ksi[0]


        dNdEta[0][1] = -0.25 * (1 - ksi[1]);                // N1 dla ksi[1] = 0,577350269
        dNdEta[1][1] = -0.25 * (1 + ksi[1]);                // N2 dla ksi[1]
        dNdEta[2][1] =  0.25 * (1 + ksi[1]);                // N3 dla ksi[1]
        dNdEta[3][1] =  0.25 * (1 - ksi[1]);                // N4 dla ksi[1]


        dNdEta[0][2] = -0.25 * (1 - ksi[2]);                // N1 dla ksi[2] = 0,577350269
        dNdEta[1][2] = -0.25 * (1 + ksi[2]);                // N2 dla ksi[2]
        dNdEta[2][2] =  0.25 * (1 + ksi[2]);                // N3 dla ksi[2]
        dNdEta[3][2] =  0.25 * (1 - ksi[2]);                // N4 dla ksi[2]


        dNdEta[0][3] = -0.25 * (1 - ksi[3]);                 // N1 dla ksi[3] = -0,577350269
        dNdEta[1][3] = -0.25 * (1 + ksi[3]);                 // N2 dla ksi[3]
        dNdEta[2][3] =  0.25 * (1 + ksi[3]);                 // N3 dla ksi[3]
        dNdEta[3][3] =  0.25 * (1 - ksi[3]);                 // N4 dla ksi[3]

    }

    // ------------------------------------------------------------------------------------------------------wypisywania
    DecimalFormat decimalFormat = new DecimalFormat("#0.0000");
    // wypisanie wartości KSI oraz ETA zgodnie z danymi z excela
    void printKsiEta() {
        System.out.println("\n\n----\t----\t----\t----\t ELEMENT UNIWERSALNY \t----\t----\t----\t----");
        System.out.println("Wypisanie wartości KSI oraz ETA (zgodnie z arkuszem: Jakobian2D)");

        System.out.println("\n\t\t\t KSI (globalnie x)");
        for (int ks = 0; ks < 4; ks++) {
            System.out.print(decimalFormat.format(ksi[ks]) + " \t\t"); }

        System.out.println("\n\n\t\t\t ETA (globalnie y)");
        for (int et = 0; et < 4; et++) {
            System.out.print(decimalFormat.format(eta[et]) + " \t\t"); }
        System.out.println("\n--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t\n");
    }

    // wypisanie macierzy
    void printMacierze_dN_dNdKsi_dNdEta() {
        //------------------------------------------------------------------------------------------------------------dN
        DecimalFormat decimalFormat = new DecimalFormat("#0.00000");
        System.out.println("\n\n----\t----\t----\t----\t ELEMENT UNIWERSALNY \t----\t----\t----\t----");
        System.out.println("--\t--\t--\t--\t--\t--\t--\t MACIERZ 1\t--\t--\t--\t--\t--\t--\t--\t--\t\n\t\t\t\t\t\t" +
                "(funkcje kształtu N)");
        System.out.println("Wartości funkcji kształtu N1: ");
        for (int i = 0; i<4; i++) { System.out.println("Węzeł " + (i+1) + " = " + decimalFormat.format(dN[i][0])); }
        System.out.println("Suma funkcji kształtu: " + decimalFormat.format(dN[0][0] + dN[1][0] + dN[2][0] + dN[3][0]));

        System.out.println("\nWartości funkcji kształtu N2: ");
        for (int i = 0; i<4; i++) { System.out.println("Węzeł " + (i+1) + " = " + decimalFormat.format(dN[i][1])); }
        System.out.println("Suma funkcji kształtu: " + decimalFormat.format(dN[0][1] + dN[1][1] + dN[2][1] + dN[3][1]));

        System.out.println("\nWartości funkcji kształtu N3: ");
        for (int i = 0; i<4; i++) { System.out.println("Węzeł " + (i+1) + " = " + decimalFormat.format(dN[i][2])); }
        System.out.println("Suma funkcji kształtu: " + decimalFormat.format(dN[0][2] + dN[1][2] + dN[2][2] + dN[3][2]));

        System.out.println("\nWartości funkcji kształtu N4: ");
        for (int i = 0; i<4; i++) { System.out.println("Węzeł " + (i+1) + " = " + decimalFormat.format(dN[i][3])); }
        System.out.println("Suma funkcji kształtu: " + decimalFormat.format(dN[0][3] + dN[1][3] + dN[2][3] + dN[3][3]));


        //--------------------------------------------------------------------------------------------------------dNdKSI
        System.out.println("\n\n--\t--\t--\t--\t--\t--\t--\t--\t MACIERZ 2\t--\t--\t--\t--\t--\t--\t--\t\t\n\t\t" +
                "\t\t\t\t(funkcje kształtu po KSI)");
        System.out.println("\t\t\t\t\t\t\t\t  dN/dKSI");
        System.out.println("Wartości funkcji kształtu dNdKsi dla ETA 1: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dKsi = " +
                decimalFormat.format(dNdKsi[i][0])); }

        System.out.println("\nWartości funkcji kształtu dNdKsi dla ETA 2: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dKsi = " +
                decimalFormat.format(dNdKsi[i][1])); }

        System.out.println("\nWartości funkcji kształtu dNdKsi dla ETA 3: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dKsi = " +
                decimalFormat.format(dNdKsi[i][2])); }

        System.out.println("\nWartości funkcji kształtu dNdKsi dla ETA 4: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dKsi = " +
                decimalFormat.format(dNdKsi[i][3])); }


        //--------------------------------------------------------------------------------------------------------dNdETA
        System.out.println("\n\n--\t--\t--\t--\t--\t--\t--\t--\t MACIERZ 3\t--\t--\t--\t--\t--\t--\t--\t\t\n\t\t" +
                "\t\t\t\t(funkcje kształtu po ETA)");
        System.out.println("\t\t\t\t\t\t\t\t  dNdETA");
        System.out.println("Wartości funkcji kształtu dNdETA dla KSI 1: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dEta = " +
                decimalFormat.format(dNdEta[i][0])); }

        System.out.println("\nWartości funkcji kształtu dNdETA dla KSI 2: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dEta = " +
                decimalFormat.format(dNdEta[i][1])); }

        System.out.println("\nWartości funkcji kształtu dNdETA dla KSI 3: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dEta = " +
                decimalFormat.format(dNdEta[i][2])); }

        System.out.println("\nWartości funkcji kształtu dNdETA dla KSI 4: ");
        for (int i = 0; i<4; i++) { System.out.println("dN" + (i+1) + "/dEta = " +
                decimalFormat.format(dNdEta[i][3])); }
    }
}
