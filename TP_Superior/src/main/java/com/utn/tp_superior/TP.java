package com.utn.tp_superior;

import static java.lang.Math.*;
import java.util.ArrayList;
import java.util.Collections;
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
    public static int maxIteraciones = 200;
    public static double TOL = pow(10, -30);
    public static double M[][];
    public static double[] B = null;
    public static float M2[][];
    public static float[] B2 = null;
    static Gauss ge = new Gauss();

    //PARAMETROS TP
    public static int tamMatriz = 100;
    public static int tamGS = 15;

    public static void main(String[] args) {

        /*
        //EJERCICIO 1
        //a)
        puntoA();
        
        //b)
        puntoB();
         */
        //EJERCICIO 2
        //a)
        make_sys((int) tamMatriz);

        //Imprimir matriz de entrada
        /*
        System.out.println("\nMatriz A: ");
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M.length; j++) {
                System.out.print(" " + M[i][j]);
            }
            System.out.print(" | " + B[i]);
            System.out.println("");
        }
         */
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
        System.out.println("\nNorma del residuo : 10^" + log(residuoMax2));
        System.out.println();
        
        //d)
        double[] x0 = new double [B.length];
        make_sys(15);
        gauss_seidel(M, B, x0);
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
        double result[] = new double[3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 1; j++) {
                suma = 0;
                for (int k = 0; k < 3; k++) {
                    suma += (A.getValueAt(i, k)) * B[k];
                }
                result[i] = suma;
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
    
    public static void gauss_seidel(double[][] m, double[] b, double[] sol){        
    }
}
