import java.text.DecimalFormat;

// Klasa agregująca macierze [C] oraz [H] do globalnych, wypisywanie wektora {P}
// ustawianie temperatury w węzłach oraz obliczenie max i minimalnych temperatur metodą gaussa
public class Proces
{
    Grid grid;
    Dane dane;
    MatrixHGlobal matrixHGlobal; // = new MatrixHGlobal(dane); --> wywala błąd
    MatrixCGlobal matrixCGlobal;
    VectorPGlobal vectorPGlobal;
    Gauss gauss = new Gauss();
    double [] vectorTempInNodes;

    public Proces() {
        grid = new Grid();                              // siatka
        dane = new Dane();                              // pobranie danych
        matrixCGlobal = new MatrixCGlobal(dane);        // C globalne
        matrixHGlobal = new MatrixHGlobal(dane);        // H globalne

        vectorPGlobal = new VectorPGlobal(dane);        // P wektor globalny
        vectorTempInNodes = new double[dane.nh];


        for (int i = 0; i<dane.ne; i++) {
            // liczenie jakobianu dla każdego elementu siatki
            Jakobian2D jakobian2D = new Jakobian2D(grid.elements[i]);
            // liczenie macierzy lokalnej H dla każdego elementu siatki korzystając z utworzonego jakobianu i danych
            MatrixCLocal matrixCLocal = new MatrixCLocal(grid.elements[i], jakobian2D, dane);       // macierz lok [C]
            MatrixHLocal matrixHLocal = new MatrixHLocal(grid.elements[i], jakobian2D, dane);       // macierz lok [H]

            // agregacja lokalnych [C] do globalnej
            matrixCGlobal.ageracjaLokalnychC(matrixCLocal);

            // agregacja z lokalnych [H] i [H_BC] do globanych macierzy [H]
            matrixHGlobal.ageracjaLokalnychH(matrixHLocal);     // agregacja do globalnego [H] lokalnych [H]
            matrixHGlobal.ageracjaLokalnychHBC(matrixHLocal);   // agregacja do globalnego [H] lokalnych [H_BC]


            // agregacja lokalnego wektora {P} tworzonego w klasie MatrixHLocal do globalnego {P}
            vectorPGlobal.agregacjaVektoraPLokalnegoDoGlobalnego(matrixHLocal);                            // wektor {P}
        }

        // wypisanie globalnej macierzy [C]
        matrixCGlobal.printMartixCGlobal();

        // wypisanie globalnej macierzy [H] (samego [H] bez uwzględniania warunków brzegowych)
        // matrixHGlobal.printMartixHGlobal();             // należy zakomentować aglomeracjeLokalnych HBC

        // wydruk macierzy globalnej [H] po dodaniu [C]/krok czasowy dT
        matrixHGlobal.ageracjaHplusCprzezKrokCzasowy(matrixCGlobal);
        System.out.println("\n\n[H] = [H] + [C]/dT");
        matrixHGlobal.printMartixHGlobal();


        // agregacja globalnego {P}
        vectorPGlobal.agregacjaWekoraPFinal(matrixCGlobal);
        vectorPGlobal.printVectorP();




// =====================================================================================================================
// ------------------------------------------------------------------------------------------------------symulacja Gauss
        System.out.println("\n\nMAKSYMALNE I MINIMALNE TEMPERATURY W KOLEJNYCH KROKACH CZASOWYCH");
        System.out.println("Time[s] \t MinTemp[s] \t MaxTemp[s]");
        for (int i = 0; i<dane.simulationTime; i += dane.simulationStepTime) {
            vectorTempInNodes = gauss.gauss(matrixHGlobal, vectorPGlobal, dane);

//            for(double n:vectorTempInNodes) { System.out.println(n); }    // wypisanie całego vectorTempInNodes

            // ustawienie temperatury w węzłach siatki korzystając z danych w vectorTempInNodes
            grid.setTemperatureOfNodes(vectorTempInNodes);

            vectorPGlobal.zerowaniePGlobalnego();                           // wyzerowanie {P} globalnego

            // wypełnienie wektora {P} nowymi wartościami
            for (int j = 0; j<dane.ne; j++) {
                Jakobian2D jakobian2D = new Jakobian2D(grid.elements[j]);
                MatrixHLocal matrixHLocal = new MatrixHLocal(grid.elements[j], jakobian2D, dane);
                vectorPGlobal.agregacjaVektoraPLokalnegoDoGlobalnego(matrixHLocal);
            }

            // podstawienie do wzoru z metodaGaussa (VectorPGlobal)
            vectorPGlobal.metodaGaussa(matrixCGlobal, vectorTempInNodes);

            // w celu sprawdzenia jak to chula dla jednej iteracji, zmiana w forze i odpalenie:
//            matrixHGlobal.printMartixHGlobal();
//            vectorPGlobal.printVectorP();

            // sortowanie temperatur
            double min, max;
            min = vectorTempInNodes[0];
            max = vectorTempInNodes[0];
            for (int k = 0; k<vectorTempInNodes.length; k++) {
                if (min > vectorTempInNodes[k]) { min = vectorTempInNodes[k]; }
                if (max < vectorTempInNodes[k]) { max = vectorTempInNodes[k]; }
            }

            DecimalFormat decimalFormat = new DecimalFormat("#0.000");
            System.out.println((i+50)+ "  \t\t" + decimalFormat.format(min) + "  \t\t" + decimalFormat.format(max));


            // żeby wygenerować jak to się nagrzewa ogólnie ---> wygenerowanie mapy ciepła w excelu
//            for(int x = 0; x<dane.nh; x++) {
//                if (x % dane.N_H == 0) { System.out.println(); }
//                System.out.format("%.2f\t", grid.nodes[x].nodeTemp);
//            } System.out.println("\n");
        }
    }
}