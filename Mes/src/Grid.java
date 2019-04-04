// siatka mesowska (ustawienie w siatce węzłów oraz elementów)

public class Grid
{
    Node [] nodes;             // tablica węzłów
    Element [] elements;       // tablica wszystkich elementów siatki mes
    Dane dane;

    public Grid() {
        this.dane = new Dane();             // wczytanie danych
        elements = new Element[dane.ne];    // tablica zawierająca ilość wszystkich elementów
        nodes = new Node[dane.nh];          // tablica  zawierająca ilość wszystkich węzłów

        int licznik = 0;    // dla numeracji podczas przechodzenia przez poszczególne elementy/węzły

        // pętla, żeby przechodzić pomiędzy węzłami w elementach (co jest w elemencie)
        for (int i = 0; i < (dane.N_H - 1); i++) {      //N_H - ilość węzłów na wysokości
            for(int j = 0; j<(dane.N_B - 1); j++) {
                elements[licznik++] = new Element(i * dane.dB, j * dane.dH, dane);
            }
        }

        licznik = 0;
        for(int i = 0; i < (dane.N_H); i++) {
            for(int j = 0; j < (dane.N_B); j++) {
                nodes[licznik++] = new Node(i * dane.dB, j * dane.dH, dane);
            }
        }
        idNodes();      // ustawienie ID dla poszczególnych węzłów
    }


    // funkcja przyporządkowująca id węzłów
    void idNodes() {
        for(int i = 0; i<dane.ne; i++) {
            for (int j = 0; j<4; j++) {
                for(int k = 0; k<dane.nh; k++) {
                    //jeśli współrzędna z elementu jest równa tej z tablicy węzłów ustawiam id
                    if (elements[i].node[j].x == nodes[k].x && elements[i].node[j].y == nodes[k].y) {
                            elements[i].node[j].ID = k;}
                }
            }
        }
    }

    // funkcja ustawiająca temperatury w węzłach siatki (ustawienie temperatury w procesie po Gausie)
    void setTemperatureOfNodes (double[] temperature) {
        for(int i = 0; i<dane.ne; i++) {
            for (int j = 0; j<4; j++) {
                for(int k = 0; k<dane.nh; k++) {
                    if (elements[i].node[j].x == nodes[k].x && elements[i].node[j].y == nodes[k].y) {
                            elements[i].node[j].nodeTemp = temperature[k];
                            nodes[k].nodeTemp = temperature[k]; }
                }
            }
        }
    }
}
