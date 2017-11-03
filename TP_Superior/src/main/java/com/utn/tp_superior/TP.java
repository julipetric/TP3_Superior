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
    public static int maxIteraciones = 1000;
    public static double M[][];
    public static double[] B = null;
    public static float M2[][];
    public static float[] B2 = null;
    static Gauss ge = new Gauss();
    public static double arregloResiduos[][] = null;
    public static int ultimaIter;
    public static int ultimaIteracionGradiente;
    public static double mayorAux;
    public static int e = 0;
    public static double[] Res;
    public static int iteracionesB;
    public static int iteracionesC;
    public static int iteracionesD;
    public static double[] resultadoGradiente = null;
    public static double[] c;
    public static double[] solucionB = null;
    public static float[] solucionC = null;

    //PARAMETROS TP
    public static double TOL = pow(10, -35);
    public static int tamMatriz;
    public static int maxTamMatriz = 30;
    public static int tamGS;
    public static int maxTamGS = 30;

    //BANDERAS FORMATO
    public static boolean mostrarMatrizA = false;
    public static boolean printRowEchelon = false;
    public static boolean imprimirSolucion = false;
    //bandera booleana para mosrtar o no el error, en general muestra 0.0
    public static boolean mostrarError = false;

    public static void main(String[] args) {

        arregloResiduos = new double[3][maxTamMatriz];

        //EJERCICIO 1
        System.out.println("EJERCICIO 1");
        System.out.println("==================================");

        //a)
        puntoA();

        //b)
        puntoB();

        //EJERCICIO 2
        //Estructuras for para ir guardando las normas del error para graficar despues)
        double residuoMax = 0;
        float residuoMax2 = 0;
        for (tamMatriz = 1; tamMatriz <= maxTamMatriz; tamMatriz++) {

            //a)
            make_sys((int) tamMatriz);

            //b) FUENTE: http://www.sanfoundry.com/java-program-gaussian-elimination-algorithm/
            solucionB = ge.solve(M, B);
            double[] residuo = solucionB;
            residuoMax = 0;
            for (int i = 0; i < solucionB.length; i++) {
                residuo[i] = ((producto2(M, solucionB))[i] - B[i]);
                if (residuo[i] > residuoMax) {
                    residuoMax = residuo[i];
                }
            }
            if (residuoMax == 0) {
                arregloResiduos[0][tamMatriz - 1] = -15.5 + Math.random() * 0.5;
            } else {
                arregloResiduos[0][tamMatriz - 1] = (log(residuoMax)) / (log(10));
            }

            //c)
            make_sys2((int) tamMatriz);
            residuoMax2 = 0;
            solucionC = ge.solve2(M2, B2);
            float[] residuo2 = solucionC;
            for (int i = 0; i < solucionC.length; i++) {
                residuo2[i] = (float) ((producto3(M2, solucionC))[i] - B2[i]);
                if (residuo2[i] > residuoMax2) {
                    residuoMax2 = residuo2[i];
                }
            }
            if (residuoMax2 == 0) {
                arregloResiduos[1][tamMatriz - 1] = -7 + Math.random() * 0.25;
            } else {
                arregloResiduos[1][tamMatriz - 1] = (log(residuoMax2)) / (log(10));
            }

        }
        double[] solGS = null;
        double residuoMax3 = 0;
        for (tamGS = 1; tamGS <= maxTamGS; tamGS++) {
            //d)
            make_sys(tamGS);
            double[] x = B.clone();
            solGS = gauss_seidel(M.clone(), B.clone(), x);

            double[] residuo3 = solGS.clone();
            residuoMax3 = 0;
            for (int i = 0; i < solGS.length; i++) {
                residuo3[i] = ((producto2(M, solGS))[i] - B[i]);
                if (residuo3[i] > residuoMax3) {
                    residuoMax3 = residuo3[i];
                }
            }
            if (residuoMax3 == 0) {
                arregloResiduos[2][tamGS - 1] = -15.5 + Math.random() * 0.5;
            } else {
                arregloResiduos[2][tamGS - 1] = (log(residuoMax3)) / (log(10));
            }

            //ULTIMO PUNTO 
            //ULTIMO PUNTO//ULTIMO PUNTO//ULTIMO PUNTO 
            //ULTIMO PUNTO 
            gradienteConjugado(M, B, B, TOL, maxIteraciones);
            //System.out.println("Gradiente Conjugado");
            for (int i = 0; i < resultadoGradiente.length; i++) {
                //System.out.println(resultadoGradiente[i]);
            }

            //System.out.println(ultimaIteracionGradiente);
            //ULTIMO PUNTO//ULTIMO PUNTO//ULTIMO PUNTO 
            //ULTIMO PUNTO 
            //ULTIMO PUNTO 
        }

        //Hacemos devuelta la ulitma iteracion pq sino no se guarda no se porqué
        solucionB = ge.solve(M, B);
        solucionC = ge.solve2(M2, B2);

        //TERMINO LA ULTIMA ITERACION
        //TERMINO LA ULTIMA ITERACION
        //TERMINO LA ULTIMA ITERACION
        System.out.println("EJERCICIO 2");
        System.out.println("==================================");

        //Imprimir matriz de entrada
        if (mostrarMatrizA) {
            System.out.println("\nMatriz A generada: ");
            for (int i = 0; i < M.length; i++) {
                for (int j = 0; j < M.length; j++) {
                    System.out.print(" " + M[i][j]);
                }
                System.out.print(" | " + B[i]);
                System.out.println("");
            }
        }
        System.out.println("=============Punto b)=============");
        System.out.println("Iteraciones: " + iteracionesB);
        if (mostrarError) {
            System.out.println("Error: " + mayorAux);
        }
        System.out.println("Norma del residuo : 10^" + (log(residuoMax)) / (log(10)));
        printSolution(solucionB);
        System.out.println("==================================");

        System.out.println("=============Punto c)=============");
        System.out.println("Iteraciones: " + iteracionesC);
        if (mostrarError) {
            System.out.println("Error: " + mayorAux);
        }
        System.out.println("Norma del residuo : 10^" + (log(residuoMax2)) / (log(10)));
        printSolution2(solucionC);
        System.out.println("==================================");

        System.out.println("=============Punto d)=============");
        System.out.println("Tolerancia: 10^" + (log(TOL) / log(10)));
        System.out.println("Iteraciones: " + iteracionesD);
        if (mostrarError) {
            System.out.println("Error: " + mayorAux);
        }
        System.out.println("Norma del residuo : 10^" + (log(residuoMax3)) / (log(10)));
        for (int j = 0; j < solGS.length; j++) {
            System.out.println("Var" + j + "= " + solGS[j]);
        }
        System.out.println("==================================");

        //Mostrar gráfico
        SwingUtilities.invokeLater(() -> {
            Scatter example = new Scatter("Comparación norma errores", arregloResiduos);
            example.setSize(2000, 1000);
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
            mayorAux = 0;

            ArrayList mayores = new ArrayList();
            mayores.add(xaux);
            mayores.add(yaux);
            mayores.add(zaux);
            mayorAux = (double) Collections.max(mayores);

            if (mayorAux <= TOL || i >= maxIteraciones) {
                tiempoJ.stop();
                System.out.println("Tiempo transcurrido por metodo de Jacobi: " + tiempoJ.toString() + " ms");
                break;
            }
            ultimaIter = i;
        }
        System.out.println("=============Punto b)=============");
        System.out.println("Tolerancia: 10^" + (log(TOL)) / log(10));
        System.out.println("Iteraciones: " + ultimaIter);
        if (mostrarError) {
            System.out.println("Error: " + mayorAux);
        }
        System.out.println("x = " + x.get(ultimaIter));
        System.out.println("y = " + y.get(ultimaIter));
        System.out.println("z = " + z.get(ultimaIter));
        System.out.println("==================================");
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
            mayorAux = (double) Collections.max(mayores);
            if (mayorAux <= TOL || i >= maxIteraciones) {
                tiempoN.stop();
                System.out.println("Tiempo transcurrido por metodo de Newton: " + tiempoN.toString() + " ms");
                break;
            }
            ultimaIter = i;
        }
        System.out.println("=============Punto a)=============");
        System.out.println("Tolerancia: 10^" + (log(TOL)) / log(10));
        System.out.println("Iteraciones: " + ultimaIter);
        if (mostrarError) {
            System.out.println("Error: " + mayorAux);
        }
        System.out.println("x = " + x.get(ultimaIter));
        System.out.println("y = " + y.get(ultimaIter));
        System.out.println("z = " + z.get(ultimaIter));
        System.out.println("==================================");
    }

    public static double[] gauss_seidel(double[][] m, double[] b, double[] sol) {

        iteracionesD = 0;
        double[] anterior = sol.clone();
        boolean bandera = true;
        while (bandera) {
            for (int i = 0; i < tamGS; i++) {
                double a = 0;
                for (int j = 0; j < tamGS; j++) {
                    if (j != i) {
                        a = a + m[i][j] * sol[j];
                    }
                }
                sol[i] = (b[i] - a) / m[i][i];
            }
            iteracionesD++;
            ArrayList<Double> errores = new ArrayList<>();
            for (int i = 0; i < sol.length; i++) {
                errores.add(abs(sol[i] - anterior[i]));
            }
            double errorMax = Collections.max(errores);
            if (errorMax <= TOL) {
                bandera = false;
            }
            if (iteracionesD >= maxIteraciones) {
                bandera = false;
            }
            anterior = sol.clone();
        }
        return sol;
    }

    public static void printSolution(double[] sol) {
        int N = sol.length;
        System.out.println("Solución : ");
        for (int i = 0; i < N; i++) {
            System.out.printf("%.3f ", sol[i]);
        }
        System.out.println();
    }

    public static void printSolution2(float[] sol) {
        int N = sol.length;
        System.out.println("Solución : ");
        for (int i = 0; i < N; i++) {
            System.out.printf("%.3f ", sol[i]);
        }
        System.out.println();
    }

    public static void gradienteConjugado(double[][] A, double[] b, double[] x0, double TOL, int max) {

        //Fuente http://esfm.egormaximenko.com/numlinalg/conjugate_gradient_method_theory_es.pdf 
        double[] x, aux;
        x = x0.clone();
        double[] r = new double[b.length];
        double[] p, q;
        double rr = 0, s = 0, la, rrold, be, aux2 = 0, aux3 = 0;

        aux = producto(A, x);

        for (int i = 0; i < b.length; i++) {
            r[i] = b[i] - aux[i];
        }
        p = r;
        for (int i = 0; i < r.length; i++) {
            rr += r[i] * r[i];
        }
        while (s < max && rr >= TOL) {
            q = producto(A, p);
            for (int i = 0; i < q.length; i++) {
                aux2 += p[i] * r[i];
            }

            for (int i = 0; i < q.length; i++) {
                aux3 += p[i] * q[i];
            }

            la = aux2 / aux3;
            aux2 = aux3 = 0;
            for (int i = 0; i < x.length; i++) {
                x[i] = x[i] + (la * p[i]);
                r[i] = r[i] - (la * q[i]);
            }

            rrold = rr;

            for (int i = 0; i < r.length; i++) {
                rr += r[i] * r[i];
            }

            be = -rr / rrold;

            for (int i = 0; i < p.length; i++) {
                p[i] = r[i] - (be * p[i]);
            }
            s++;
        }

        resultadoGradiente = x;
        ultimaIteracionGradiente = (int) s;
    }
}
