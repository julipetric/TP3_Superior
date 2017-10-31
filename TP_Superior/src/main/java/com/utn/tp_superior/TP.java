package com.utn.tp_superior;

import static java.lang.Math.*;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import org.apache.commons.lang3.time.StopWatch;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juli
 */
public class TP {

    //VARIABLES GLOBALES
    public static int maxIteraciones = 100000;
    public static double TOL = pow(10, -100);
    public static double M[][];
    public static double[] B = null;
    public static float M2[][];
    public static float[] B2 = null;
    static Gauss ge = new Gauss();
    public static double arregloResiduos[][] = null;

    //PARAMETROS TP
    public static int tamMatriz;
    public static int maxIteracionesMatriz = 30;
    public static int tamGS = 30;

    public static void main(String[] args) {

        arregloResiduos = new double[3][maxIteracionesMatriz];

        /*
        //EJERCICIO 1
        //a)
        puntoA();
        
        //b)
        puntoB();
         */
        //EJERCICIO 2
        //Estructura for para ir guardando las normas del error para graficar despues)
        for (tamMatriz = 1; tamMatriz <= maxIteracionesMatriz; tamMatriz++) {

            //a)
            make_sys((int) tamMatriz);

            //Imprimir matriz de entrada
            System.out.println("\nMatriz A: ");
            for (int i = 0; i < M.length; i++) {
                for (int j = 0; j < M.length; j++) {
                    System.out.print(" " + M[i][j]);
                }
                System.out.print(" | " + B[i]);
                System.out.println("");
            }

            //b) FUENTE: http://www.sanfoundry.com/java-program-gaussian-elimination-algorithm/
            double[] solucion;
            solucion = ge.solve(M, B);
            double[] residuo = solucion;
            double residuoMax = 0;
            for (int i = 0; i < solucion.length; i++) {
                residuo[i] = ((producto2(M, solucion))[i] - B[i]);
                if (residuo[i] > residuoMax) {
                    residuoMax = residuo[i];
                }
            }
            if (residuoMax == 0) {
                arregloResiduos[0][tamMatriz - 1] = -35+Math.random()*0.5;
            } else {
                arregloResiduos[0][tamMatriz - 1] = log(residuoMax);
            }
            System.out.println("\nNorma del residuo : 10^" + log(residuoMax));
            System.out.println();

            //c)
            make_sys2((int) tamMatriz);
            float[] solucion2;
            float residuoMax2 = 0;
            solucion2 = ge.solve2(M2, B2);
            float[] residuo2 = solucion2;
            for (int i = 0; i < solucion2.length; i++) {
                residuo2[i] = (float) ((producto3(M2, solucion2))[i] - B2[i]);
                if (residuo2[i] > residuoMax2) {
                    residuoMax2 = residuo2[i];
                }
            }
            if (residuoMax2 == 0) {
                arregloResiduos[1][tamMatriz - 1] = -15+Math.random()*0.25;
            } else {
                arregloResiduos[1][tamMatriz - 1] = log(residuoMax2);
            }

            System.out.println("\nNorma del residuo : 10^" + log(residuoMax2));
            System.out.println();
        }
        //d)
        //make_sys(tamGS);
        double[][] me = {{3, 1, 1}, {1, 3, 1}, {2, 1, 4}};
        double[] dea = {4, 3, 2};

        double[] x = new double[3];
        for (int i = 0; i < 3; i++) {
            x[i] = 1;
        }

        gauss_seidel(me, dea, x);

        SwingUtilities.invokeLater(() -> {
            Scatter example = new Scatter("Comparación norma errores", arregloResiduos);
            example.setSize(800, 400);
            example.setLocationRelativeTo(null);
            example.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
            example.setVisible(true);
        });

    }

    //Funcion para crear la matriz segun especificación del enunciado dado un tamaño
    public static void make_sys(int n) {
        M = new double[n][n];
        B = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int a = abs(i - j);

                if (i == j) {
                    M[i][j] = 1;
                }

                if (i > j) {
                    M[i][j] = (4 + a) / pow((2 + a), 2);
                }

                if (i < j) {
                    M[i][j] = (4 + a) / pow((2 + a), 2);
                }

            }
        }
        for (int y = 0; y < n; y++) {
            B[y] = (2 * y + 1);
        }
    }

    //Igual a la anterior, pero de simple precision
    public static void make_sys2(int n) {
        M2 = new float[n][n];
        B2 = new float[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int a = abs(i - j);

                if (i == j) {
                    M2[i][j] = 1;
                }

                if (i > j) {
                    M2[i][j] = (float) ((4 + a) / (float) pow((2 + a), 2));
                }

                if (i < j) {
                    M2[i][j] = (float) ((4 + a) / (float) pow((2 + a), 2));
                }

            }
        }
        for (int y = 0; y < n; y++) {
            B2[y] = (float) (2 * y + 1);
        }
    }

    //Funcion producto matrix (tipo Matrix) y vector
    public static double[] producto(Matrix A, double[] B) {
        double suma;
        double result[] = new double[B.length];
        for (int i = 0; i < B.length; i++) {
            suma = 0;
            for (int k = 0; k < B.length; k++) {
                suma += (A.getValueAt(i, k)) * B[k];
            }
            result[i] = suma;
        }
        return result;
    }

    //Funcion producto matrix (tipo Matrix) y Matriz
    public static double[][] producto(Matrix A, double[][] B) {
        double suma;
        double result[][] = new double[B.length][B.length];
        for (int i = 0; i < B.length; i++) {
            for (int j = 0; j < B.length; j++) {
                suma = 0;
                for (int k = 0; k < B.length; k++) {
                    suma += (A.getValueAt(i, k)) * B[k][j];
                }
                result[i][j] = suma;
            }
        }
        return result;
    }

    //Igual a producto(Matrix, double[]), pero cambia tipo Matrix por arreglo bidimensional
    public static double[] producto2(double[][] A, double[] B) {
        double suma;
        int cols = A[0].length;
        int rows = A.length;
        double result[] = new double[cols];
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                suma = 0;
                for (int k = 0; k < cols; k++) {
                    suma += (A[i][k]) * B[k];
                }
                result[i] = suma;
            }
        }
        return result;
    }

    //Igual a la anterior, pero de simple precision
    public static float[] producto3(float[][] A, float[] B) {
        float suma;
        int cols = A[0].length;
        int rows = A.length;
        float result[] = new float[cols];
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                suma = 0;
                for (int k = 0; k < cols; k++) {
                    suma += (A[i][k]) * B[k];
                }
                result[i] = suma;
            }
        }
        return result;
    }

    //Funcion producto Matriz y Vector
    private static double[] producto(double[][] A, double[] b) {
        double suma = 0;
        double result[] = new double[b.length];
        for (int i = 0; i < b.length; i++) {
            suma = 0;
            for (int k = 0; k < b.length; k++) {
                suma += A[i][k] * b[k];
            }
            result[i] = suma;
        }
        return result;

    }

    public static void puntoB() {

        /*
.4.3. Método de Jacobi
En este método la matriz tangente que interviene en cada iteración
se
sustituye por otra con la misma diagonal pero con todos sus demás elementos
nulos. Más concretamente, denotando por
a la matriz:
Esta forma de proceder efectivamente reduce de forma notable el número de
operaciones (sólo conlleva evaluar n funciones derivadas en lugar de n2 y además
la inversión de una matriz diagonal sólo implica n operaciones). Pero sólo es válida
si los elementos no diagonales de la matriz jacobiana son ”pequeños” comparados
con los términos diagonales.
3.4.4.


    http://ocw.upm.es/matematica-aplicada/programacion-y-metodos-numericos/contenidos/TEMA_8/Apuntes/EcsNoLin.pdf
         */
        StopWatch tiempoJ = new StopWatch();
        tiempoJ.reset();
        tiempoJ.start();
        ArrayList x = new ArrayList();
        ArrayList y = new ArrayList();
        ArrayList z = new ArrayList();
        x.add((double) 1);
        y.add((double) 10);
        z.add((double) 20);

        double ins[][] = {{1, 1, 1},
        {2 * (double) x.get(0), 2 * (double) y.get(0), 2 * (double) z.get(0)},
        {-exp(-((double) x.get(0))) * pow((double) y.get(0), 3) + (double) y.get(0) + (double) z.get(0), 3 * exp(-((double) x.get(0))) * pow((double) y.get(0), 2) + (double) x.get(0), (double) x.get(0)}};
        Matrix Jacobiana = new Matrix(ins);

        for (int i = 1; i < 3; i++) {
            Jacobiana.setValueAt(0, i, 0);
        }

        Jacobiana.setValueAt(1, 0, 0);
        Jacobiana.setValueAt(1, 2, 0);

        for (int i = 0; i < 2; i++) {
            Jacobiana.setValueAt(2, i, 0);
        }

        Jacobiana = Matrix.inverse(Jacobiana);

        double F[] = {(double) x.get(0) + (double) y.get(0) + (double) z.get(0) + Math.PI,
            pow((double) x.get(0), 2) + pow((double) y.get(0), 2) + pow((double) z.get(0), 2) - 7,
            pow((double) y.get(0), 3) * exp(-((double) x.get(0))) + (double) x.get(0) * ((double) y.get(0) + (double) z.get(0)) + 1};

        for (int i = 1; i <= maxIteraciones; i++) {
            double M[] = producto(Jacobiana, F);
            x.add(i, (double) x.get(i - 1) - M[0]);
            y.add(i, (double) y.get(i - 1) - M[1]);
            z.add(i, (double) z.get(i - 1) - M[2]);

            double ins2[][] = {{1, 1, 1},
            {2 * (double) x.get(i), 2 * (double) y.get(i), 2 * (double) z.get(i)},
            {-exp(-((double) x.get(i))) * pow((double) y.get(i), 3) + (double) y.get(i) + (double) z.get(i), 3 * exp(-((double) x.get(i))) * pow((double) y.get(i), 2) + (double) x.get(i), (double) x.get(i)}};
            Jacobiana = new Matrix(ins2);
            Jacobiana = Matrix.inverse(Jacobiana);
            double F2[] = {(double) x.get(i) + (double) y.get(i) + (double) z.get(i) + 7,
                pow((double) x.get(i), 2) + pow((double) y.get(i), 2) + pow((double) z.get(i), 2) - 49,
                pow((double) y.get(i), 3) * exp(-((double) x.get(i))) + (double) x.get(i) * ((double) y.get(i) + (double) z.get(i)) + 1};
            F = F2;

            double xaux = abs((double) x.get(i) - (double) x.get(i - 1));
            double yaux = abs((double) y.get(i) - (double) y.get(i - 1));
            double zaux = abs((double) z.get(i) - (double) z.get(i - 1));
            double mayor = 0;

            ArrayList mayores = new ArrayList();
            mayores.add(xaux);
            mayores.add(yaux);
            mayores.add(zaux);
            mayor = (double) Collections.max(mayores);
            System.out.println("=======Punto b)=========");
            System.out.println("Tolerancia: " + TOL);
            System.out.println("Iteraciones: " + i);
            System.out.println("Error: " + mayor);
            System.out.println("x = " + x.get(i));
            System.out.println("y = " + y.get(i));
            System.out.println("z = " + z.get(i));
            System.out.println("================================");

            if (mayor <= TOL || i >= maxIteraciones) {
                tiempoJ.stop();
                System.out.println("Tiempo transcurrido por metodo de Jacobi: " + tiempoJ.toString() + " ms");
                break;
            }
        }
    }

    public static void puntoA() {
        StopWatch tiempoN = new StopWatch();
        tiempoN.reset();
        tiempoN.start();
        ArrayList x = new ArrayList();
        ArrayList y = new ArrayList();
        ArrayList z = new ArrayList();
        x.add((double) 1);
        y.add((double) 10);
        z.add((double) 20);

        double ins[][] = {{1, 1, 1},
        {2 * (double) x.get(0), 2 * (double) y.get(0), 2 * (double) z.get(0)},
        {-exp(-((double) x.get(0))) * pow((double) y.get(0), 3) + (double) y.get(0) + (double) z.get(0), 3 * exp(-((double) x.get(0))) * pow((double) y.get(0), 2) + (double) x.get(0), (double) x.get(0)}};
        Matrix Jacobiana = new Matrix(ins);
        Jacobiana = Matrix.inverse(Jacobiana);

        double F[] = {(double) x.get(0) + (double) y.get(0) + (double) z.get(0) + Math.PI,
            pow((double) x.get(0), 2) + pow((double) y.get(0), 2) + pow((double) z.get(0), 2) - 7,
            pow((double) y.get(0), 3) * exp(-((double) x.get(0))) + (double) x.get(0) * ((double) y.get(0) + (double) z.get(0)) + 1};

        for (int i = 1; i <= maxIteraciones; i++) {
            double M[] = producto(Jacobiana, F);
            x.add(i, (double) x.get(i - 1) - M[0]);
            y.add(i, (double) y.get(i - 1) - M[1]);
            z.add(i, (double) z.get(i - 1) - M[2]);

            double ins2[][] = {{1, 1, 1},
            {2 * (double) x.get(i), 2 * (double) y.get(i), 2 * (double) z.get(i)},
            {-exp(-((double) x.get(i))) * pow((double) y.get(i), 3) + (double) y.get(i) + (double) z.get(i), 3 * exp(-((double) x.get(i))) * pow((double) y.get(i), 2) + (double) x.get(i), (double) x.get(i)}};
            Jacobiana = new Matrix(ins2);
            Jacobiana = Matrix.inverse(Jacobiana);
            double F2[] = {(double) x.get(i) + (double) y.get(i) + (double) z.get(i) + 7,
                pow((double) x.get(i), 2) + pow((double) y.get(i), 2) + pow((double) z.get(i), 2) - 49,
                pow((double) y.get(i), 3) * exp(-((double) x.get(i))) + (double) x.get(i) * ((double) y.get(i) + (double) z.get(i)) + 1};
            F = F2;

            double xaux = abs((double) x.get(i) - (double) x.get(i - 1));
            double yaux = abs((double) y.get(i) - (double) y.get(i - 1));
            double zaux = abs((double) z.get(i) - (double) z.get(i - 1));
            double mayor = 0;

            ArrayList mayores = new ArrayList();
            mayores.add(xaux);
            mayores.add(yaux);
            mayores.add(zaux);
            mayor = (double) Collections.max(mayores);
            System.out.println("=======Punto a)=========");
            System.out.println("Tolerancia: " + TOL);
            System.out.println("Iteraciones: " + i);
            System.out.println("Error: " + mayor);
            System.out.println("x = " + x.get(i));
            System.out.println("y = " + y.get(i));
            System.out.println("z = " + z.get(i));
            System.out.println("================================");

            if (mayor <= TOL || i >= maxIteraciones) {
                tiempoN.stop();
                System.out.println("Tiempo transcurrido por metodo de Newton: " + tiempoN.toString() + " ms");
                break;
            }
        }
    }

    public static void gauss_seidel(double[][] m, double[] b, double[] sol) {

        double aux = 0;
        double[][] p = new double[b.length][b.length];

        for (int i = 0; i < b.length; i++) {
            for (int j = 0; j < b.length; j++) {
                p[i][j] = m[i][j];
            }
        }

        double[][] M;
        Matrix auxMatriz;
        double[] c, vecAux;
        double[] Res = new double[b.length];
        double[] Aux = new double[b.length];
        double[] anterior = new double[b.length];
        boolean bandera = true;

        // x^(k+1)=M*x^(k)+c
        //FUENTE https://es.wikipedia.org/wiki/M%C3%A9todo_de_Gauss-Seidel
        //Se verifica si es diagonal dominante
        for (int i = 0; i < b.length; i++) {
            for (int j = 0; j < b.length; j++) {
                if (i != j) {
                    aux += abs(m[i][j]);
                }
            }
            if (abs(m[i][i]) <= aux) {
                System.out.println("NO ES MATRIZ DOMINANTE");
                bandera = false;
                break;
            }
            aux = 0;
        }

        if (bandera) {
            //Matriz N
            for (int i = 0; i < b.length; i++) {
                for (int j = 0; j < b.length; j++) {
                    if (i < j) {
                        m[i][j] = 0;
                    }
                }
            }
            //Matriz N^-1
            auxMatriz = new Matrix(m);
            auxMatriz = Matrix.inverse(auxMatriz);
            //Matriz P
            for (int i = 0; i < b.length; i++) {
                for (int j = 0; j < b.length; j++) {
                    if (i >= j) {
                        p[i][j] = 0;
                    }
                }
            }

            //Matriz M
            M = producto(auxMatriz, p);
            //Vector c
            c = producto(auxMatriz, b);
            //Primera aproximacion
            vecAux = producto(M, sol);

            for (int i = 0; i < c.length; i++) {
                Res[i] = vecAux[i] + c[i];
            }
            for (int i = 0; i < c.length; i++) {
                anterior[i] = 0;
            }
            ArrayList mayores = new ArrayList();
            for (int i = 0; i < maxIteraciones; i++) {
                double mayor = 0;
                mayores.clear();
                for (int d = 0; d < Res.length; d++) {
                    Aux[d] = abs((double) Res[d] - (double) anterior[d]);
                    mayores.add(Aux[d]);
                }

                mayor = (double) Collections.max(mayores);
                System.out.println("=======Punto d)=========");
                System.out.println("Tolerancia: " + TOL);
                int e = i + 1;
                System.out.println("Iteraciones: " + e);
                System.out.println("Error: " + mayor);
                for (int j = 0; j < c.length; j++) {
                    System.out.println("Var" + j + "= " + Res[j]);
                }
                System.out.println("================================");

                if (mayor <= TOL || i >= maxIteraciones) {
                    break;
                }

                vecAux = producto(M, Aux);

                for (int w = 0; w < vecAux.length; w++) {
                    anterior[w] = Res[w];

                }

                for (int k = 0; k < c.length; k++) {
                    Res[k] = vecAux[k] + c[k];
                }
            }
        }
    }
}
