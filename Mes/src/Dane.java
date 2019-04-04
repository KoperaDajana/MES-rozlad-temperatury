import java.text.DecimalFormat;

// wczytywanie danych (initial data)
public class Dane
{
    int initialTemp;        // temperatura początkowa
    int simulationTime;     // czas trwania symulacji
    int simulationStepTime; // krok czasowy
    int ambientTemp;        // temperatura otoczenia
    int alfa;               // współczynnik przewodzenia ciepła k
    double H, B;            // H - wysokość przekroju, B - szerokość przekroju
    int N_H, N_B;           // N_H - ilość węzłów na wysokości, N_B -ilość węzłów na szerokości
    int specificHeat;       // pojemność cieplna
    int conductivity;       // przewodzenie
    int density;            // gęstość
    double dH, dB;          // odległość między węzłami (dH - na wysokości, dB - na szerokości)
    int nh;                 // ilość węzłów OGÓLNIE
    int ne;                 // ilość wszystkich elementów OGÓLNIE


    public Dane()
    {
        this.initialTemp = 100;
        this.simulationTime = 500;
        this.simulationStepTime = 50;
        this.ambientTemp = 1200;
        this.alfa = 300;
        H = 0.1;
        B = 0.1;
        N_H = 4;
        N_B = 4;
        this.specificHeat = 700;
        this.conductivity = 25;
        this.density = 7800;

        nh = N_H * N_B;                 // ilość węzłów (wysokość x szerokość)
        ne = (N_H - 1) * (N_B - 1);     // ilość wszystkich elementów
        dB = B / (N_B  - 1);            // odległości między węzłami po szerokości
        dH = H / (N_H - 1);             // odległości między węzłami po wysokości

    }

    DecimalFormat decimalFormat = new DecimalFormat("#0.0000");
    public void printData() {
        System.out.println("\n----\t----\t----\t\tPARAMETRY (initial data)\t\t----\t----\t----" +
                "\nTemperatura początkowa:  " + initialTemp + "\t[C]" +
                "\nCzas trwania symulacji:  " + simulationTime + "\t[s]" +
                "\nKrok czasowy symulacji:  " + simulationStepTime + "\t\t[s]" +
                "\nTemperatura otoczenia:   " + ambientTemp + "\t[C]" +
                "\nWspółczynnik Alfa:       " + alfa + "\t[W/m^2k]" +
                "\nH (wysokość przekroju):  " + H + "\t[m]" +
                "\nB (szerokość przekroju): " + B + "\t[m]" +
                "\nN_H (liczba węzłów na wysokości):  " + N_H +
                "\nN_B (liczba węzłów na szerokości): " + N_B +
                "\nPojemność cieplna: " + specificHeat + "\t[J/kg*C]" +
                "\nPrzewodzenie:      " + conductivity + "\t[W/m*C]" +
                "\nGęstość:           " + density + "\t[kg/m^3]" +
                "\n\nLiczba wszystkich węzłów:   " + nh +
                "\nLiczba wszstkich elementów: " + ne +
                "\nOdległość pomiędzy węzłami (szerokość): " + decimalFormat.format(dB) +
                "\nOdległość pomiędzy węzłami (wysokość):  " + decimalFormat.format(dH) +
                "\n----\t----\t----\t----\t----\t----\t----\t----\t----\t----\t----");
    }
}
