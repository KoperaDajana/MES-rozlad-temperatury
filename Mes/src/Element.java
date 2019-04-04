//klasa reprezentująca Element, stworzenie elementu, wypisanie go

public class Element
{
    public Node[] node;     // tablica węzłów do tworzenia elementu
    double x, y;            // współrzędne, żeby wiedzieć jak się poruszać po siatce mesowskiej

    public Element(double x, double y, Dane dane)
    {
        this.x = x;
        this.y = y;
        // utworzenie tablicy reprezentującej współrzędne elementu
        node = new Node[4];

        node[0] = new Node(x, y, dane);                             // x = 0,        y =0
        node[1] = new Node(x + dane.dB, y, dane);                // x + w prawo,  y
        node[2] = new Node(x + dane.dB, y + dane.dH, dane);   // x + w prawo,  y + do góry
        node[3] = new Node(x, y + dane.dH, dane);                // x,            y + do góry
    }

//        void printElement() {
//        System.out.println("\n\n----\t----\t----\tPRÓBNE WYPISANIE ELEMENTU\t----\t----\t----");
//        System.out.println("Początek siatki x: " + x + ", y: " + y);
//        System.out.println("Węzeł 1 : \t(x: " + node[0].x + ", y: " + node[0].y + ")" + "\tID: " + node[0].ID +
//                            "\n\t\t\tWarunek brzegowy:    " + node[0].status +
//                            "\n\t\t\tTemperatura w węźle: " + node[0].nodeTemp);
//        System.out.println("\nWęzeł 2: \t(x: " + node[1].x + ", y: " + node[1].y + ")" + "\tID: " + node[1].ID +
//                            "\n\t\t\tWarunek brzegowy:    " + node[1].status +
//                            "\n\t\t\tTemperatura w węźle: " + node[1].nodeTemp);
//        System.out.println("\nWęzeł 3: \t(x: " + node[2].x + ", y: " + node[2].y + ")" + "\tID: " + node[2].ID +
//                            "\n\t\t\tWarunek brzegowy:    " + node[2].status +
//                            "\n\t\t\tTemperatura w węźle: " + node[2].nodeTemp);
//        System.out.println("\nWęzeł 4: \t(x: " + node[3].x + ", y: " + node[3].y + ")" + "\tID: " + node[3].ID +
//                            "\n\t\t\tWarunek brzegowy:    " + node[3].status +
//                            "\n\t\t\tTemperatura w węźle: " + node[3].nodeTemp);
//        System.out.println("----\t----\t----\t----\t----\t----\t----\t----\t----\t----\t----");
//    }
}
