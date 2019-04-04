// klasa reprezentująca węzeł, informacje o nim - czy zawiera warunek brzegowy oraz jaką temperature

public class Node
{
    double x, y;        // współrzędne węzła
    int status;         // pokazuje, czy na danym węźle jest warunek brzegowy (0 - nie, jeśli 1 - tak)
    int ID;
    double nodeTemp;    // temperatura jaka jest w danym węźle


    public Node(double x, double y, Dane dane)
    {
        this.x = x;
        this.y = y;

        // te warunki sprawdzane w MatrixHLocal podczas tworzenia macierzy H z sum (proces końcowy lokalnego H)
        // jeżeli x = 0 lub y = 0 lub x>= szerokość lub y>= wysokość to istnieje warunek brzegowy w węźle
        if (this.x == 0. || this.y == 0. || this.x >= dane.B || this.y >= dane.H) {
            status = 1;
        // jeśli w "środku" to nie ma warunku brzegowego
        } else { status = 0;}
    }
}
