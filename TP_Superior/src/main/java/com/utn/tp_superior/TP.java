package com.utn.tp_superior;

import static java.lang.Math.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;
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

    public static int maxIteraciones = 200;
    public static double TOL = pow(10, -130);
    public static double M[][];
    public static double[] B = null;
    public static void main(String[] args) {

  //============================================EJERCICIO 1===========================================
        //============================================PUNTO a)===========================================
        //puntoA();
        //============================================PUNTO a)===========================================

        //============================================PUNTO b)===========================================
        //puntoB();
        //============================================PUNTO b)===========================================
 //============================================EJERCICIO 2===========================================
        //============================================PUNTO a)===========================================
        int tamMatriz = 1000;
        make_sys((int)tamMatriz);
        
        System.out.println("\nMatriz A: "); 
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M.length; j++) {
                 System.out.print(" "+M[i][j]);
            }
            System.out.print(" | "+B[i]);
            System.out.println("");
        }  
                
        //============================================PUNTO b)===========================================
        
        // FUENTE DE ESTO http://www.sanfoundry.com/java-program-gaussian-elimination-algorithm/
        Gauss ge = new Gauss();
        ge.solve(M,B);  
        
        //============================================PUNTO c)===========================================
        
    }      
 
public static void make_sys(int n){
    double a=0;
    M= new double[n][n];
    B= new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a = abs(i-j);
                
                if(i==j){
                    M[i][j]=1;
                }
                
                if(i>j){
                    M[i][j]=(4+a)/pow((2+a),2);
                   }
                
                if(i<j){
                    M[i][j]=(4+a)/pow((2+a),2);
                   }
                
            }
        }
        for (int y = 0; y<n; y++){
            B[y]=(2*y+1);
        }
    }
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
public static void puntoB() {

        /*
.4.3. M´etodo de Jacobi
En este m´etodo la matriz tangente que interviene en cada iteraci´on
se
sustituye por otra con la misma diagonal pero con todos sus dem´as elementos
nulos. M´as concretamente, denotando por
a la matriz:
Esta forma de proceder efectivamente reduce de forma notable el n´umero de
operaciones (s´olo conlleva evaluar n funciones derivadas en lugar de n2 y adem´as
la inversi´on de una matriz diagonal s´olo implica n operaciones). Pero s´olo es v´alida
si los elementos no diagonales de la matriz jacobiana son ”pequeños” comparados
con los t´erminos diagonales.
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

        double F[] = {(double) x.get(0) + (double) y.get(0) + (double) z.get(0) + 7,
            pow((double) x.get(0), 2) + pow((double) y.get(0), 2) + pow((double) z.get(0), 2) - 49,
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

        double F[] = {(double) x.get(0) + (double) y.get(0) + (double) z.get(0) + 7,
            pow((double) x.get(0), 2) + pow((double) y.get(0), 2) + pow((double) z.get(0), 2) - 49,
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
}
