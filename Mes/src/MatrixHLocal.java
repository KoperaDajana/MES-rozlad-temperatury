import java.text.DecimalFormat;
import static java.lang.Math.sqrt;

// arkusz z Excela: Matrix_H oraz Matrix_H_BC_2d, posiada macierz H lokalną oraz wektor {P} lokalny
public class MatrixHLocal
{
    Element element;
    // wartości pobierane z dNdx, "wiersze", ponieważ ważniejsze punkty całkowania
    double [][] dN1dNT_dx_1pc; // dla 1 punktu całkowania N x N' po x (arkusz 2 excel MACIERZ H)
    double [][] dN2dNT_dx_2pc; // dla 2 punktu całkowania N x N' po x
    double [][] dN3dNT_dx_3pc; // dla 3 punktu całkowania N x N' po x
    double [][] dN4dNT_dx_4pc; // dla 4 punktu całkowania N x N' po x

    double [][] dN1dNT_dy_1pc; // dla 1 punktu całkowania N x N' po y (arkusz 2 excel MACIERZ H)
    double [][] dN2dNT_dy_2pc; // dla 2 punktu całkowania N x N' po y
    double [][] dN3dNT_dy_3pc; // dla 3 punktu całkowania N x N' po y
    double [][] dN4dNT_dy_4pc; // dla 4 punktu całkowania N x N' po y

    // ---------------------------------------------------------------------------------------------------------- * detJ
    double [][] dN1dNT_dx_1pc_detJ; // dla 1 punktu całkowania N x N' po x (arkusz 2 excel MACIERZ H)
    double [][] dN2dNT_dx_2pc_detJ; // dla 2 punktu całkowania N x N' po x
    double [][] dN3dNT_dx_3pc_detJ; // dla 3 punktu całkowania N x N' po x
    double [][] dN4dNT_dx_4pc_detJ; // dla 4 punktu całkowania N x N' po x

    double [][] dN1dNT_dy_1pc_detJ; // dla 1 punktu całkowania N x N' po y (arkusz 2 excel MACIERZ H)
    double [][] dN2dNT_dy_2pc_detJ; // dla 2 punktu całkowania N x N' po y
    double [][] dN3dNT_dy_3pc_detJ; // dla 3 punktu całkowania N x N' po y
    double [][] dN4dNT_dy_4pc_detJ; // dla 4 punktu całkowania N x N' po y

    // --------------------------------------------------------------------------------------------------- k * XY * detJ
    double [][] kXY_detJ_1pc;
    double [][] kXY_detJ_2pc;
    double [][] kXY_detJ_3pc;
    double [][] kXY_detJ_4pc;

    double [][] matrixH;    // pierwsza część macierzy H (ta niezależna, bez warunków brzegowych)
    // koniec na rzecz arkusza excel: MatrixH---------------------------------------------------------------------------






    // tworzone na rzecz excel:===========================================================================>MatrixH_BC_2D
    // druga część macierzy H, która uwzględnia warunek brzegowy (pierwszą część warunku brzegowego konwekcji)
    double [][] matrixHBC;

    double [][] NHBC;       // w celu uzupełnienia tabel "suma" w excelu arkusz Matrix_H_BC, funkcje kształtu
    double [][] NHBC_1pc;   // dla 1 punktu całkowania
    double [][] NHBC_2pc;   // dla 2 punktu całkowania

    double [] dbAndDetJ;    // wektor reprezentujący wartość długości boku oraz wyznacznika

    double [][] suma1pc;    // sumy z Matrix_HCB żeby stworzyć H z warunkami brzegowymi
    double [][] suma2pc;
    double [][] suma3pc;
    double [][] suma4pc;


    // =========================================================================================================WEKTOR P
    double [] vektorPLokalnie1pc;
    double [] vektorPLokalnie2pc;
    double [] vektorPLokalnie3pc;
    double [] vektorPLokalnie4pc;

    double [] vectorPLocal = new double[]{0, 0, 0, 0}; //wypełniam zerami wektor p na starcie

    public MatrixHLocal(Element element, Jakobian2D jakobian2D, Dane dane) {
        // wszystkie new w konstruktorze, żeby dopiero teraz mi przypisało pamięć, co by jej nie zapychało tak
        vektorPLokalnie1pc = new double[4];
        vektorPLokalnie2pc = new double[4];
        vektorPLokalnie3pc = new double[4];
        vektorPLokalnie4pc = new double[4];
        zerowanieWektoraP();

        this.element = element;

        // N x N'(Transponowane) dla 1 punktu całkowania (po x)----------------------------------------------------1pc/x
        dN1dNT_dx_1pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN1dNT_dx_1pc[0][i] = jakobian2D.dNdx[0][i] * jakobian2D.dNdx[0][0];
            dN1dNT_dx_1pc[1][i] = jakobian2D.dNdx[0][i] * jakobian2D.dNdx[0][1];
            dN1dNT_dx_1pc[2][i] = jakobian2D.dNdx[0][i] * jakobian2D.dNdx[0][2];
            dN1dNT_dx_1pc[3][i] = jakobian2D.dNdx[0][i] * jakobian2D.dNdx[0][3];
        }
        // N x N' dla 2 punktu całkowania (po x)-------------------------------------------------------------------2pc/x
        dN2dNT_dx_2pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN2dNT_dx_2pc[0][i] = jakobian2D.dNdx[1][i] * jakobian2D.dNdx[1][0];
            dN2dNT_dx_2pc[1][i] = jakobian2D.dNdx[1][i] * jakobian2D.dNdx[1][1];
            dN2dNT_dx_2pc[2][i] = jakobian2D.dNdx[1][i] * jakobian2D.dNdx[1][2];
            dN2dNT_dx_2pc[3][i] = jakobian2D.dNdx[1][i] * jakobian2D.dNdx[1][3];
        }
        // N x N' dla 3 punktu całkowania (po x)-------------------------------------------------------------------3pc/x
        dN3dNT_dx_3pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN3dNT_dx_3pc[0][i] = jakobian2D.dNdx[2][i] * jakobian2D.dNdx[2][0];
            dN3dNT_dx_3pc[1][i] = jakobian2D.dNdx[2][i] * jakobian2D.dNdx[2][1];
            dN3dNT_dx_3pc[2][i] = jakobian2D.dNdx[2][i] * jakobian2D.dNdx[2][2];
            dN3dNT_dx_3pc[3][i] = jakobian2D.dNdx[2][i] * jakobian2D.dNdx[2][3];
        }
        // N x N' dla 4 punktu całkowania (po x)-------------------------------------------------------------------4pc/x
        dN4dNT_dx_4pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN4dNT_dx_4pc[0][i] = jakobian2D.dNdx[3][i] * jakobian2D.dNdx[3][0];
            dN4dNT_dx_4pc[1][i] = jakobian2D.dNdx[3][i] * jakobian2D.dNdx[3][1];
            dN4dNT_dx_4pc[2][i] = jakobian2D.dNdx[3][i] * jakobian2D.dNdx[3][2];
            dN4dNT_dx_4pc[3][i] = jakobian2D.dNdx[3][i] * jakobian2D.dNdx[3][3];
        }
        // =============================================================================================================
        // N x N'(Transponowane) dla 1 punktu całkowania (po x)----------------------------------------------------1pc/y
        dN1dNT_dy_1pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN1dNT_dy_1pc[0][i] = jakobian2D.dNdy[0][i] * jakobian2D.dNdy[0][0];
            dN1dNT_dy_1pc[1][i] = jakobian2D.dNdy[0][i] * jakobian2D.dNdy[0][1];
            dN1dNT_dy_1pc[2][i] = jakobian2D.dNdy[0][i] * jakobian2D.dNdy[0][2];
            dN1dNT_dy_1pc[3][i] = jakobian2D.dNdy[0][i] * jakobian2D.dNdy[0][3];
        }
        // N x N' dla 2 punktu całkowania (po x)-------------------------------------------------------------------2pc/y
        dN2dNT_dy_2pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN2dNT_dy_2pc[0][i] = jakobian2D.dNdy[1][i] * jakobian2D.dNdy[1][0];
            dN2dNT_dy_2pc[1][i] = jakobian2D.dNdy[1][i] * jakobian2D.dNdy[1][1];
            dN2dNT_dy_2pc[2][i] = jakobian2D.dNdy[1][i] * jakobian2D.dNdy[1][2];
            dN2dNT_dy_2pc[3][i] = jakobian2D.dNdy[1][i] * jakobian2D.dNdy[1][3];
        }
        // N x N' dla 3 punktu całkowania (po x)-------------------------------------------------------------------3pc/y
        dN3dNT_dy_3pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN3dNT_dy_3pc[0][i] = jakobian2D.dNdy[2][i] * jakobian2D.dNdy[2][0];
            dN3dNT_dy_3pc[1][i] = jakobian2D.dNdy[2][i] * jakobian2D.dNdy[2][1];
            dN3dNT_dy_3pc[2][i] = jakobian2D.dNdy[2][i] * jakobian2D.dNdy[2][2];
            dN3dNT_dy_3pc[3][i] = jakobian2D.dNdy[2][i] * jakobian2D.dNdy[2][3];
        }
        // N x N' dla 4 punktu całkowania (po x)-------------------------------------------------------------------4pc/y
        dN4dNT_dy_4pc = new double[4][4];
        for (int i = 0; i < 4; i++) {
            // kolumnami
            dN4dNT_dy_4pc[0][i] = jakobian2D.dNdy[3][i] * jakobian2D.dNdy[3][0];
            dN4dNT_dy_4pc[1][i] = jakobian2D.dNdy[3][i] * jakobian2D.dNdy[3][1];
            dN4dNT_dy_4pc[2][i] = jakobian2D.dNdy[3][i] * jakobian2D.dNdy[3][2];
            dN4dNT_dy_4pc[3][i] = jakobian2D.dNdy[3][i] * jakobian2D.dNdy[3][3];
        }


        // --------------------------------------DETJ--------------------------------------------------------------DET J
        dN1dNT_dx_1pc_detJ = new double[4][4];
        dN2dNT_dx_2pc_detJ = new double[4][4];
        dN3dNT_dx_3pc_detJ = new double[4][4];
        dN4dNT_dx_4pc_detJ = new double[4][4];
        // dla po x
        for (int i = 0; i<4; i++) {
            dN1dNT_dx_1pc_detJ[0][i] = dN1dNT_dx_1pc[0][i] * jakobian2D.detJ[i];
            dN1dNT_dx_1pc_detJ[1][i] = dN1dNT_dx_1pc[1][i] * jakobian2D.detJ[i];
            dN1dNT_dx_1pc_detJ[2][i] = dN1dNT_dx_1pc[2][i] * jakobian2D.detJ[i];
            dN1dNT_dx_1pc_detJ[3][i] = dN1dNT_dx_1pc[3][i] * jakobian2D.detJ[i];
        }
        for (int i = 0; i<4; i++) {
            dN2dNT_dx_2pc_detJ[0][i] = dN2dNT_dx_2pc[0][i] * jakobian2D.detJ[i];
            dN2dNT_dx_2pc_detJ[1][i] = dN2dNT_dx_2pc[1][i] * jakobian2D.detJ[i];
            dN2dNT_dx_2pc_detJ[2][i] = dN2dNT_dx_2pc[2][i] * jakobian2D.detJ[i];
            dN2dNT_dx_2pc_detJ[3][i] = dN2dNT_dx_2pc[3][i] * jakobian2D.detJ[i];
        }
        for (int i = 0; i<4; i++) {
            dN3dNT_dx_3pc_detJ[0][i] = dN3dNT_dx_3pc[0][i] * jakobian2D.detJ[i];
            dN3dNT_dx_3pc_detJ[1][i] = dN3dNT_dx_3pc[1][i] * jakobian2D.detJ[i];
            dN3dNT_dx_3pc_detJ[2][i] = dN3dNT_dx_3pc[2][i] * jakobian2D.detJ[i];
            dN3dNT_dx_3pc_detJ[3][i] = dN3dNT_dx_3pc[3][i] * jakobian2D.detJ[i];
        }
        for (int i = 0; i<4; i++) {
            dN4dNT_dx_4pc_detJ[0][i] = dN4dNT_dx_4pc[0][i] * jakobian2D.detJ[i];
            dN4dNT_dx_4pc_detJ[1][i] = dN4dNT_dx_4pc[1][i] * jakobian2D.detJ[i];
            dN4dNT_dx_4pc_detJ[2][i] = dN4dNT_dx_4pc[2][i] * jakobian2D.detJ[i];
            dN4dNT_dx_4pc_detJ[3][i] = dN4dNT_dx_4pc[3][i] * jakobian2D.detJ[i];
        }

        //det J po y
        dN1dNT_dy_1pc_detJ = new double[4][4];
        dN2dNT_dy_2pc_detJ = new double[4][4];
        dN3dNT_dy_3pc_detJ = new double[4][4];
        dN4dNT_dy_4pc_detJ = new double[4][4];

        for (int i = 0; i<4; i++) {
            dN1dNT_dy_1pc_detJ[0][i] = dN1dNT_dy_1pc[0][i] * jakobian2D.detJ[i];
            dN1dNT_dy_1pc_detJ[1][i] = dN1dNT_dy_1pc[1][i] * jakobian2D.detJ[i];
            dN1dNT_dy_1pc_detJ[2][i] = dN1dNT_dy_1pc[2][i] * jakobian2D.detJ[i];
            dN1dNT_dy_1pc_detJ[3][i] = dN1dNT_dy_1pc[3][i] * jakobian2D.detJ[i];
        }
        for (int i = 0; i<4; i++) {
            dN2dNT_dy_2pc_detJ[0][i] = dN2dNT_dy_2pc[0][i] * jakobian2D.detJ[i];
            dN2dNT_dy_2pc_detJ[1][i] = dN2dNT_dy_2pc[1][i] * jakobian2D.detJ[i];
            dN2dNT_dy_2pc_detJ[2][i] = dN2dNT_dy_2pc[2][i] * jakobian2D.detJ[i];
            dN2dNT_dy_2pc_detJ[3][i] = dN2dNT_dy_2pc[3][i] * jakobian2D.detJ[i];
        }
        for (int i = 0; i<4; i++) {
            dN3dNT_dy_3pc_detJ[0][i] = dN3dNT_dy_3pc[0][i] * jakobian2D.detJ[i];
            dN3dNT_dy_3pc_detJ[1][i] = dN3dNT_dy_3pc[1][i] * jakobian2D.detJ[i];
            dN3dNT_dy_3pc_detJ[2][i] = dN3dNT_dy_3pc[2][i] * jakobian2D.detJ[i];
            dN3dNT_dy_3pc_detJ[3][i] = dN3dNT_dy_3pc[3][i] * jakobian2D.detJ[i];
        }
        for (int i = 0; i<4; i++) {
            dN4dNT_dy_4pc_detJ[0][i] = dN4dNT_dy_4pc[0][i] * jakobian2D.detJ[i];
            dN4dNT_dy_4pc_detJ[1][i] = dN4dNT_dy_4pc[1][i] * jakobian2D.detJ[i];
            dN4dNT_dy_4pc_detJ[2][i] = dN4dNT_dy_4pc[2][i] * jakobian2D.detJ[i];
            dN4dNT_dy_4pc_detJ[3][i] = dN4dNT_dy_4pc[3][i] * jakobian2D.detJ[i];
        }

        kXY_detJ_1pc = new double[4][4];
        kXY_detJ_2pc = new double[4][4];
        kXY_detJ_3pc = new double[4][4];
        kXY_detJ_4pc = new double[4][4];

        // dla pierwszego punktu całkowania
        for (int j = 0; j<4; j++) {
            kXY_detJ_1pc[0][j] = dane.conductivity * (dN1dNT_dx_1pc_detJ[0][j] + dN1dNT_dy_1pc_detJ[0][j]);
            kXY_detJ_1pc[1][j] = dane.conductivity * (dN1dNT_dx_1pc_detJ[1][j] + dN1dNT_dy_1pc_detJ[1][j]);
            kXY_detJ_1pc[2][j] = dane.conductivity * (dN1dNT_dx_1pc_detJ[2][j] + dN1dNT_dy_1pc_detJ[2][j]);
            kXY_detJ_1pc[3][j] = dane.conductivity * (dN1dNT_dx_1pc_detJ[3][j] + dN1dNT_dy_1pc_detJ[3][j]);
        }
         // dla drugiego punktu całkowania
        for (int j = 0; j<4; j++) {
            kXY_detJ_2pc[0][j] = dane.conductivity * (dN2dNT_dx_2pc_detJ[0][j] + dN2dNT_dy_2pc_detJ[0][j]);
            kXY_detJ_2pc[1][j] = dane.conductivity * (dN2dNT_dx_2pc_detJ[1][j] + dN2dNT_dy_2pc_detJ[1][j]);
            kXY_detJ_2pc[2][j] = dane.conductivity * (dN2dNT_dx_2pc_detJ[2][j] + dN2dNT_dy_2pc_detJ[2][j]);
            kXY_detJ_2pc[3][j] = dane.conductivity * (dN2dNT_dx_2pc_detJ[3][j] + dN2dNT_dy_2pc_detJ[3][j]);
        }
        // dla trzeciego punktu całkowania
        for (int j = 0; j<4; j++) {
            kXY_detJ_3pc[0][j] = dane.conductivity * (dN3dNT_dx_3pc_detJ[0][j] + dN3dNT_dy_3pc_detJ[0][j]);
            kXY_detJ_3pc[1][j] = dane.conductivity * (dN3dNT_dx_3pc_detJ[1][j] + dN3dNT_dy_3pc_detJ[1][j]);
            kXY_detJ_3pc[2][j] = dane.conductivity * (dN3dNT_dx_3pc_detJ[2][j] + dN3dNT_dy_3pc_detJ[2][j]);
            kXY_detJ_3pc[3][j] = dane.conductivity * (dN3dNT_dx_3pc_detJ[3][j] + dN3dNT_dy_3pc_detJ[3][j]);
        }
        // dla czwartego punktu całkowania
        for (int j = 0; j<4; j++) {
            kXY_detJ_4pc[0][j] = dane.conductivity * (dN4dNT_dx_4pc_detJ[0][j] + dN4dNT_dy_4pc_detJ[0][j]);
            kXY_detJ_4pc[1][j] = dane.conductivity * (dN4dNT_dx_4pc_detJ[1][j] + dN4dNT_dy_4pc_detJ[1][j]);
            kXY_detJ_4pc[2][j] = dane.conductivity * (dN4dNT_dx_4pc_detJ[2][j] + dN4dNT_dy_4pc_detJ[2][j]);
            kXY_detJ_4pc[3][j] = dane.conductivity * (dN4dNT_dx_4pc_detJ[3][j] + dN4dNT_dy_4pc_detJ[3][j]);
        }


        matrixH = new double[4][4];     // część macierzy [H], która nie uwzględnia warunków brzegowych
        for (int i = 0; i<4; i++) {
            matrixH[0][i] = kXY_detJ_1pc[0][i] + kXY_detJ_2pc[0][i] + kXY_detJ_3pc[0][i] + kXY_detJ_4pc[0][i];
            matrixH[1][i] = kXY_detJ_1pc[1][i] + kXY_detJ_2pc[1][i] + kXY_detJ_3pc[1][i] + kXY_detJ_4pc[1][i];
            matrixH[2][i] = kXY_detJ_1pc[2][i] + kXY_detJ_2pc[2][i] + kXY_detJ_3pc[2][i] + kXY_detJ_4pc[2][i];
            matrixH[3][i] = kXY_detJ_1pc[3][i] + kXY_detJ_2pc[3][i] + kXY_detJ_3pc[3][i] + kXY_detJ_4pc[3][i];
        }

// =====================================================================================================================
// =====================================================================================================================
        // MATRIX H BC (druga część warunku brzegowego uwzględniający warunki brzegowe)
        matrixHBC = new double[4][4];
        dbAndDetJ = new double[2];      // wektor zawierający długość boku oraz długość_boku/2
        dbAndDetJ[0] = dane.dB;
        dbAndDetJ[1] = dane.dB/2;

        double ksi1 = (-1)/sqrt(3);     double ksi2 =  (1/sqrt(3));
        double eta1 = -1;               double eta2 = -1;

        NHBC = new double[2][4];
        NHBC_1pc = new double [4][4];
        NHBC_2pc = new double [4][4];
        //--------------------------------------------------------------------------------------------------------SUMA 1
        // dla 1 powierzchni wypełnienie funkcji kształtu dla nowo zdefiniowanych ksi i eta
        zerowanieNHBC();
        suma1pc = new double[4][4];

        NHBC[0][0] = 0.25 * (1 - ksi1) * (1- eta1);//1
        NHBC[1][0] = 0.25 * (1 - ksi2) * (1 - eta2);
        NHBC[0][1] = 0.25 * (1 + ksi1) * (1 - eta1);//2
        NHBC[1][1] = 0.25 * (1 + ksi2) * (1 - eta2);
        NHBC[0][2] = 0.25 * (1 + ksi1) * (1 + eta1);//3
        NHBC[1][2] = 0.25 * (1 + ksi2) * (1 + eta2);
        NHBC[0][3] = 0.25 * (1 - ksi1) * (1 + eta1);//4
        NHBC[1][3] = 0.25 * (1 - ksi2) * (1 + eta2);
        // printHBC();
        for (int i = 0; i<4; i++) {
            NHBC_1pc[0][i] = NHBC[0][i] * NHBC[0][0] * dane.alfa;
            NHBC_1pc[1][i] = NHBC[0][i] * NHBC[0][1] * dane.alfa;
            NHBC_1pc[2][i] = NHBC[0][i] * NHBC[0][2] * dane.alfa;
            NHBC_1pc[3][i] = NHBC[0][i] * NHBC[0][3] * dane.alfa;
        }
        // printHBC_1pc();
        for (int i = 0; i<4; i++) {
            NHBC_2pc[0][i] = NHBC[1][i] * NHBC[1][0] * dane.alfa;
            NHBC_2pc[1][i] = NHBC[1][i] * NHBC[1][1] * dane.alfa;
            NHBC_2pc[2][i] = NHBC[1][i] * NHBC[1][2] * dane.alfa;
            NHBC_2pc[3][i] = NHBC[1][i] * NHBC[1][3] * dane.alfa;
        }
        // printHBC_2pc();

        // i - wiersz, j - kolumno, wypełniane osobno suma1pc
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                suma1pc[i][j] = (NHBC_1pc[i][j] + NHBC_2pc[i][j]) * dbAndDetJ[1];
            }
            // uzupełnienie wartości lokalnego wektora {P} dla pierwszego punktu całkowania -------------------> {P} 1pc
            vektorPLokalnie1pc[i] = (NHBC[0][i] + NHBC[1][i]) * dbAndDetJ[1];
        }
        // printSuma1();

        //--------------------------------------------------------------------------------------------------------SUMA 2
        zerowanieNHBC(); //pow2
        suma2pc = new double[4][4];

        double ksi11 = 1;                       double ksi22 = 1;
        double eta11 = (-1)/sqrt(3);            double eta22 = 1/sqrt(3);

        NHBC[0][0] = 0.25 * (1 - ksi11) * (1- eta11);//1
        NHBC[1][0] = 0.25 * (1 - ksi22) * (1 - eta22);
        NHBC[0][1] = 0.25 * (1 + ksi11) * (1 - eta11);//2
        NHBC[1][1] = 0.25 * (1 + ksi22) * (1 - eta22);
        NHBC[0][2] = 0.25 * (1 + ksi11) * (1 + eta11);//3
        NHBC[1][2] = 0.25 * (1 + ksi22) * (1 + eta22);
        NHBC[0][3] = 0.25 * (1 - ksi11) * (1 + eta11);//4
        NHBC[1][3] = 0.25 * (1 - ksi22) * (1 + eta22);
        // printHBC();
        for (int i = 0; i<4; i++) {
            NHBC_1pc[0][i] = NHBC[0][i] * NHBC[0][0] * dane.alfa;
            NHBC_1pc[1][i] = NHBC[0][i] * NHBC[0][1] * dane.alfa;
            NHBC_1pc[2][i] = NHBC[0][i] * NHBC[0][2] * dane.alfa;
            NHBC_1pc[3][i] = NHBC[0][i] * NHBC[0][3] * dane.alfa;
        }
        // printHBC_1pc();
        for (int i = 0; i<4; i++) {
            NHBC_2pc[0][i] = NHBC[1][i] * NHBC[1][0] * dane.alfa;
            NHBC_2pc[1][i] = NHBC[1][i] * NHBC[1][1] * dane.alfa;
            NHBC_2pc[2][i] = NHBC[1][i] * NHBC[1][2] * dane.alfa;
            NHBC_2pc[3][i] = NHBC[1][i] * NHBC[1][3] * dane.alfa;
        }
        // printHBC_2pc();

        // i - wiersz, j - kolumno, wypełniane osobno dla suma2pc
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                suma2pc[i][j] = (NHBC_1pc[i][j] + NHBC_2pc[i][j]) * dbAndDetJ[1];
            }
            // uzupełnienie wartości lokalnego wektora {P} dla drugiego punktu całkowania ---------------------> {P} 2pc
            vektorPLokalnie2pc[i] = (NHBC[0][i] + NHBC[1][i]) * dbAndDetJ[1];
        }
        // printSuma2();

        //--------------------------------------------------------------------------------------------------------SUMA 3
        zerowanieNHBC(); //pow3
        suma3pc = new double[4][4];

        double ksi111 = 1/sqrt(3);      double ksi222 = (-1)/sqrt(3);
        double eta111 = 1;              double eta222 = 1;

        NHBC[0][0] = 0.25 * (1 - ksi111) * (1- eta111);//1
        NHBC[1][0] = 0.25 * (1 - ksi222) * (1 - eta222);
        NHBC[0][1] = 0.25 * (1 + ksi111) * (1 - eta111);//2
        NHBC[1][1] = 0.25 * (1 + ksi222) * (1 - eta222);
        NHBC[0][2] = 0.25 * (1 + ksi111) * (1 + eta111);//3
        NHBC[1][2] = 0.25 * (1 + ksi222) * (1 + eta222);
        NHBC[0][3] = 0.25 * (1 - ksi111) * (1 + eta111);//4
        NHBC[1][3] = 0.25 * (1 - ksi222) * (1 + eta222);
        // printHBC();
        for (int i = 0; i<4; i++) {
            NHBC_1pc[0][i] = NHBC[0][i] * NHBC[0][0] * dane.alfa;
            NHBC_1pc[1][i] = NHBC[0][i] * NHBC[0][1] * dane.alfa;
            NHBC_1pc[2][i] = NHBC[0][i] * NHBC[0][2] * dane.alfa;
            NHBC_1pc[3][i] = NHBC[0][i] * NHBC[0][3] * dane.alfa;
        }
        // printHBC_1pc();
        for (int i = 0; i<4; i++) {
            NHBC_2pc[0][i] = NHBC[1][i] * NHBC[1][0] * dane.alfa;
            NHBC_2pc[1][i] = NHBC[1][i] * NHBC[1][1] * dane.alfa;
            NHBC_2pc[2][i] = NHBC[1][i] * NHBC[1][2] * dane.alfa;
            NHBC_2pc[3][i] = NHBC[1][i] * NHBC[1][3] * dane.alfa;
        }
        // printHBC_2pc();

        // i - wiersz, j - kolumno, wypełniane osobno dla suma3pc
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                suma3pc[i][j] = (NHBC_1pc[i][j] + NHBC_2pc[i][j]) * dbAndDetJ[1];
            }
            // uzupełnienie wartości lokalnego wektora {P} dla trzeciego punktu całkowania --------------------> {P} 3pc
            vektorPLokalnie3pc[i] = (NHBC[0][i] + NHBC[1][i]) * dbAndDetJ[1];
        }
        // printSuma3();

        //--------------------------------------------------------------------------------------------------------SUMA 4
        zerowanieNHBC(); //pow4
        suma4pc = new double[4][4];

        double ksi1111 = -1;                    double ksi2222 = -1;
        double eta1111 = 1/sqrt(3);             double eta2222 = (-1)/sqrt(3);
        NHBC[0][0] = 0.25 * (1 - ksi1111) * (1- eta1111);//1
        NHBC[1][0] = 0.25 * (1 - ksi2222) * (1 - eta2222);
        NHBC[0][1] = 0.25 * (1 + ksi1111) * (1 - eta1111);//2
        NHBC[1][1] = 0.25 * (1 + ksi2222) * (1 - eta2222);
        NHBC[0][2] = 0.25 * (1 + ksi1111) * (1 + eta1111);//3
        NHBC[1][2] = 0.25 * (1 + ksi2222) * (1 + eta2222);
        NHBC[0][3] = 0.25 * (1 - ksi1111) * (1 + eta1111);//4
        NHBC[1][3] = 0.25 * (1 - ksi2222) * (1 + eta2222);
        // printHBC();
        for (int i = 0; i<4; i++) {
            NHBC_1pc[0][i] = NHBC[0][i] * NHBC[0][0] * dane.alfa;
            NHBC_1pc[1][i] = NHBC[0][i] * NHBC[0][1] * dane.alfa;
            NHBC_1pc[2][i] = NHBC[0][i] * NHBC[0][2] * dane.alfa;
            NHBC_1pc[3][i] = NHBC[0][i] * NHBC[0][3] * dane.alfa;
        }
        // printHBC_1pc();
        for (int i = 0; i<4; i++) {
            NHBC_2pc[0][i] = NHBC[1][i] * NHBC[1][0] * dane.alfa;
            NHBC_2pc[1][i] = NHBC[1][i] * NHBC[1][1] * dane.alfa;
            NHBC_2pc[2][i] = NHBC[1][i] * NHBC[1][2] * dane.alfa;
            NHBC_2pc[3][i] = NHBC[1][i] * NHBC[1][3] * dane.alfa;
        }
        // printHBC_2pc();

        // i - wiersz, j - kolumno, wypełniane dla suma4pc
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                suma4pc[i][j] = (NHBC_1pc[i][j] + NHBC_2pc[i][j]) * dbAndDetJ[1];
            }
            // uzupełnienie wartości lokalnego wektora {P} dla czwartego punktu całkowania --------------------> {P} 4pc
            vektorPLokalnie4pc[i] = (NHBC[0][i] + NHBC[1][i]) * dbAndDetJ[1];
        }
        // printSuma4();



//======================================================================================================================
//======================================================================================================================
        // MACIERZ HBC (zerowanie, sprawdzanie warunków brzegowych oraz uzupełnianie jej
        // wartościami z macierz suma1/2/3/4pc), tabela po prawej stronie MatrixH w arkuszu: MatrixH_CB_2D
        zerowanieH();
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                if(element.node[0].status == 1 && element.node[1].status == 1) {
                    matrixHBC[i][j] = matrixHBC[i][j] + suma1pc[i][j]; }

                if(element.node[1].status == 1 && element.node[2].status == 1) {
                    matrixHBC[i][j] = matrixHBC[i][j] + suma2pc[i][j]; }

                if(element.node[2].status == 1 && element.node[3].status == 1) {
                    matrixHBC[i][j] = matrixHBC[i][j] + suma3pc[i][j]; }

                if(element.node[3].status == 1 && element.node[0].status == 1) {
                    matrixHBC[i][j] = matrixHBC[i][j] + suma4pc[i][j]; }
            }

            // sprawdzanie warunków dla wektora P (czy jest w danym węźle warunek brzegowy)
            // żeby na końcu móc przemnożyć go przez współczynniki dla danego warunku
            if(element.node[0].status == 1 && element.node[1].status == 1) {
                vectorPLocal[i] += vektorPLokalnie1pc[i]; }

            if(element.node[1].status == 1 && element.node[2].status == 1) {
                vectorPLocal[i] += vektorPLokalnie2pc[i]; }

            if(element.node[2].status == 1 && element.node[3].status == 1) {
                vectorPLocal[i] += vektorPLokalnie3pc[i]; }

            if(element.node[3].status == 1 && element.node[0].status == 1) {
                vectorPLocal[i] += vektorPLokalnie4pc[i]; }
        }

        //przemnożenia wektora P przez alfa i temperaturę otoczenia (warunek brzegowy konwekcji II część)
        for(int i=0; i<4; i++){
            // przez wartość powierzchni (S) zostało pomnożone wcześniej w pętlach dla sum
            vectorPLocal[i] = vectorPLocal[i] * dane.alfa * dane.ambientTemp;
        }
    }





    // -------------------------------------------------------------------------------------------------funkcje zerujące
    void zerowanieNHBC() {
        for (int i = 0; i<2; i++) {
            for (int j = 0; j<4; j++) { NHBC[i][j] = 0; }
        }
        for (int g = 0; g<4; g++) {
            for (int u = 0; u<4; u++) {
                NHBC_1pc[g][u] = 0;
                NHBC_2pc[g][u] = 0; }
        }
    }

    void zerowanieH() {
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) { matrixHBC[i][j] = 0; }
        }
    }

    void zerowanieWektoraP() {
        for (int i = 0; i<4; i++) {
            vektorPLokalnie1pc[i] = 0;
            vektorPLokalnie2pc[i] = 0;
            vektorPLokalnie3pc[i] = 0;
            vektorPLokalnie4pc[i] = 0;
        }
    }




// ==========================================================================================================WYPISYWANIE
    DecimalFormat decimalFormat = new DecimalFormat("#0.000");
    void printMatrixHLokalnie() {
        System.out.println("\n\n----\t----\t----\t--\t MATRIX H LOKALNIE\t--\t----\t----\t----");
        System.out.println("\n==\t==\t==\t==\t==\t MACIERZ H LOKALNA\t==\t==\t==\t==\t==\t==\t==");
        System.out.println("Z excel: arkusz MatrixH (bez uwzględnienia warunków brzegowych)");
        for (int i = 0; i < 4; i++) {
            System.out.println( decimalFormat.format(matrixH[i][0]) + "\t\t" +
                                decimalFormat.format(matrixH[i][1]) + "\t\t" +
                                decimalFormat.format(matrixH[i][2]) + "\t\t" +
                                decimalFormat.format(matrixH[i][3]) + "\t\t"); }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
    }

    // wypisanie ostatniej macierzy H z excel: arkusz MatrixH_BC_2d
    void printHBC_final() {
        System.out.println("\n--\t--\t--\t--\t--\t--\t--\t Macierz H_BC\t--\t--\t--\t--\t--\t--\t--\t");
        System.out.println("Z excel: arkusz Matrix H_BC (uwzględniająca warunek brzegowy)");
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                System.out.print(decimalFormat.format(matrixHBC[i][j]) + "\t\t"); }
            System.out.println();
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println();
    }


// =====================================================================================================================
// ==========================================================================wypisywanie dla macierzy H funkcji kształtu
    void printHBC() {
        System.out.println("Macierz funkcji kształtu N HBC");
        for (int i = 0; i<2; i++) {
            for (int j = 0; j<4; j++) {
                System.out.print(decimalFormat.format(NHBC[i][j]) + "\t\t"); }
            System.out.println();
        }
        System.out.println();
    }

    void printHBC_1pc() {
        System.out.println("Macierz funkcji kształtu N HBC dla 1 pc");
            for (int i = 0; i<4; i++) {
                for (int j = 0; j<4; j++) {
                    System.out.print(decimalFormat.format(NHBC_1pc[i][j]) + "\t\t"); }
                System.out.println();
            }
        System.out.println();
    }

    void printHBC_2pc() {
        System.out.println("Macierz funkcji kształtu N HBC dla 2 pc");
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                System.out.print(decimalFormat.format(NHBC_2pc[i][j]) + "\t\t"); }
            System.out.println();
        }
        System.out.println();
    }
    // wypisanie sum, które później zostają wpisane do macierzy H końcowej w MatrixH_BC_2d
    void printSuma1() {
        System.out.println("SUMA 1");
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                System.out.print(decimalFormat.format(suma1pc[i][j]) + "\t\t"); }
            System.out.println();
        }
        System.out.println();
    }
    void printSuma2() {
        System.out.println("SUMA 2");
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                System.out.print(decimalFormat.format(suma2pc[i][j]) + "\t\t"); }
            System.out.println();
        }
        System.out.println();
    }
    void printSuma3() {
        System.out.println("SUMA 3");
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                System.out.print(decimalFormat.format(suma3pc[i][j]) + "\t\t"); }
            System.out.println();
        }
        System.out.println();
    }
    void printSuma4() {
        System.out.println("SUMA 4");
        for (int i = 0; i<4; i++) {
            for (int j = 0; j<4; j++) {
                System.out.print(decimalFormat.format(suma4pc[i][j]) + "\t\t"); }
            System.out.println();
        }
        System.out.println();
    }


    // wypisanie rzeczy z arkusza Matrix H przed lokalną agregacją z elementów w lokalną macierz H
    void printMatrixHprzedLokalnaAgregacja() {
        // ----------------------------------------------------------------------------------------------WYPISANIE PO DX
        System.out.println("\n\nX --\t--\t--\t--\t--\t--\t--\t--\t [dN][dN]'/dx X \t--\t--\t--\t--\t--\t--\t--\t--\tX");
        System.out.println("Dla 1pc po x");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN1dNT_dx_1pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dx_1pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dx_1pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dx_1pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 2pc po x");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN2dNT_dx_2pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dx_2pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dx_2pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dx_2pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 3pc po x");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN3dNT_dx_3pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dx_3pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dx_3pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dx_3pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 4pc po x");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN4dNT_dx_4pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dx_4pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dx_4pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dx_4pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("\n\nY --\t--\t--\t--\t--\t--\t--\t--\t [dN][dN]'/dy Y \t--\t--\t--\t--\t--\t--\t--\t--\tY");
        System.out.println("Dla 1pc po y");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN1dNT_dy_1pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dy_1pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dy_1pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dy_1pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 2pc po y");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN2dNT_dy_2pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dy_2pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dy_2pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dy_2pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 3pc po y");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN3dNT_dy_3pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dy_3pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dy_3pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dy_3pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        System.out.println("Dla 4pc po y");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN4dNT_dy_4pc[i][0]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dy_4pc[i][1]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dy_4pc[i][2]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dy_4pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        //========================================================================================================DET Jx
        System.out.println("\n\n--\t--\t--\t--\t--\t--\t--\t [dN][dN]'/dx * det J\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 1pc po x * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN1dNT_dx_1pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dx_1pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dx_1pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dx_1pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 2pc po x * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN2dNT_dx_2pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dx_2pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dx_2pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dx_2pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 3pc po x * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN3dNT_dx_3pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dx_3pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dx_3pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dx_3pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 4pc po x * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN4dNT_dx_4pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dx_4pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dx_4pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dx_4pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        // -------------------------------------------------------------------------------------------------------DET Jy
        System.out.println("\n\n--\t--\t--\t--\t--\t--\t--\t [dN][dN]'/dy * det J\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 1pc po y * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN1dNT_dy_1pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dy_1pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dy_1pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN1dNT_dy_1pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 2pc po y * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN2dNT_dy_2pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dy_2pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dy_2pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN2dNT_dy_2pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 3pc po y * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN3dNT_dy_3pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dy_3pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dy_3pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN3dNT_dy_3pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 4pc po y * det J");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(dN4dNT_dy_4pc_detJ[i][0]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dy_4pc_detJ[i][1]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dy_4pc_detJ[i][2]) + "\t\t" +
                    decimalFormat.format(dN4dNT_dy_4pc_detJ[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");

        // -------------------------------------------------------------------------k * ([N][N]'/dx + [N][N]'/dy) * detJ
        System.out.println("\n\n--\t--\t--\t--\t--\t k * ([N][N]'/dx + [N][N]'/dy) * detJ\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 1pc: ");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(kXY_detJ_1pc[i][0]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_1pc[i][1]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_1pc[i][2]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_1pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 2pc: ");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(kXY_detJ_2pc[i][0]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_2pc[i][1]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_2pc[i][2]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_2pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 3pc: ");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(kXY_detJ_3pc[i][0]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_3pc[i][1]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_3pc[i][2]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_3pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
        System.out.println("Dla 4pc: ");
        for (int i = 0; i < 4; i++) {
            System.out.println(decimalFormat.format(kXY_detJ_4pc[i][0]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_4pc[i][1]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_4pc[i][2]) + "\t\t" +
                    decimalFormat.format(kXY_detJ_4pc[i][3]) + "\t\t");
        }
        System.out.println("--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--");
    }


}
