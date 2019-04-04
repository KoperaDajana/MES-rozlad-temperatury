public class Main
{
    public static void main(String[] args)
    {
        // stworzenie obiektu danych
//        Dane dane = new Dane();
//        dane.dB = 0.025;
//        dane.dH = 0.025;
//        dane.printData();   //wypisanie parametrów procesu


        // próbne stworzenie zwykłego elementu
//        Element element = new Element(0,0, dane);
        // element.printElement();

        // wypisanie Ksi Eta oraz macierzy z Jakobian2D
//        UniversalElement universalElement = new UniversalElement();
        // universalElement.printKsiEta();
        // universalElement.printMacierze_dN_dNdKsi_dNdEta();

        // utworzenie jakobianu
//        Jakobian2D jakobian2D = new Jakobian2D(element);
        // jakobian2D.printJakobian();                     // wersja dla Jakobian2D
        // jakobian2D.printJakobianGlobal();               // wynik jakobianu dla współrzędnych globalnych

        // stworzenie siatki MES
//        Grid grid = new Grid();


        // lokalna macierz [H]
//        MatrixHLocal matrixHLocal = new MatrixHLocal(element, jakobian2D, dane);

        // wypisanie lokalnej macierzy H z arkusza excel: Matrix H
//        matrixHLocal.printMatrixHLokalnie();
//        matrixHLocal.printHBC_final();

        // lokalna macierz [C]
//        MatrixCLocal matrixCLocal = new MatrixCLocal(element, jakobian2D, dane);
        // matrixCLocal.printMatrixCforPC();                   // wypisanie 4 macierzy osobno dla każdego pc
//        matrixCLocal.printLocalMatrixC();                   // wypisanie lokalnej macierzy [C]


// -----------------------------------------------------------------koniec rzeczy tworzonych w celu sprawdzenia w excelu
// ---------------------------------------------------------------------------------------------------------------------
        // globalna macierz [C] zostaje zagregowana w Proces
        // globalna macierz [H] zostaje zagregowana w Proces
        // wektor globlany {P} wypisywany w Proces
        // metoda Gaussa również w Proces
        Proces proces = new Proces();
    }
}
