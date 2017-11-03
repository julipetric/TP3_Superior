/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.utn.tp_superior;

/**
 *
 * @author Juli
 */
public class Gauss extends TP {

    public double[] solve(double[][] A, double[] B) {
        int N = B.length;

        for (int k = 0; k < N; k++) {
            /**
             * find pivot row *
             */
            int max = k;
            for (int i = k + 1; i < N; i++) {
                if (Math.abs(A[i][k]) > Math.abs(A[max][k])) {
                    max = i;
                }
            }

            /**
             * swap row in A matrix *
             */
            double[] temp = A[k];
            A[k] = A[max];
            A[max] = temp;

            /**
             * swap corresponding values in constants matrix *
             */
            double t = B[k];
            B[k] = B[max];
            B[max] = t;

            /**
             * pivot within A and B *
             */
            for (int i = k + 1; i < N; i++) {
                double factor = A[i][k] / A[k][k];
                B[i] -= factor * B[k];
                for (int j = k; j < N; j++) {
                    A[i][j] -= factor * A[k][j];
                }
            }
            iteracionesB = k;
        }

        /**
         * Print row echelon form *
         */
        if (printRowEchelon) {
            printRowEchelonForm(A, B);
        }

        /**
         * back substitution *
         */
        double[] solution = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * solution[j];
            }
            solution[i] = (B[i] - sum) / A[i][i];
        }
        /**
         * Print solution *
         */
        if(imprimirSolucion)printSolution(solution);

        return solution;
    }

    public float[] solve2(float[][] A, float[] B) {
        int N = B.length;

        for (int k = 0; k < N; k++) {
            /**
             * find pivot row *
             */
            int max = k;
            for (int i = k + 1; i < N; i++) {
                if (Math.abs(A[i][k]) > Math.abs(A[max][k])) {
                    max = i;
                }
            }

            /**
             * swap row in A matrix *
             */
            float[] temp = A[k];
            A[k] = A[max];
            A[max] = temp;

            /**
             * swap corresponding values in constants matrix *
             */
            float t = B[k];
            B[k] = B[max];
            B[max] = t;

            /**
             * pivot within A and B *
             */
            for (int i = k + 1; i < N; i++) {
                float factor = (float) (A[i][k] / A[k][k]);
                B[i] -= factor * B[k];
                for (int j = k; j < N; j++) {
                    A[i][j] -= factor * A[k][j];
                }
            }
            iteracionesC = k;
        }
        /**
         * Print row echelon form *
         */
        if (printRowEchelon) {
            printRowEchelonForm2(A, B);
        }

        /**
         * back substitution *
         */
        float[] solution = new float[N];
        for (int i = N - 1; i >= 0; i--) {
            float sum = 0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * solution[j];
            }
            solution[i] = (float) ((B[i] - sum) / A[i][i]);
        }
        /**
         * Print solution *
         */
        if (imprimirSolucion) printSolution2(solution);
 
        return solution;
    }

    public void printRowEchelonForm(double[][] A, double[] B) {
        int N = B.length;
        System.out.println("\nRow Echelon form : ");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%.3f ", A[i][j]);
            }
            System.out.printf("| %.3f\n", B[i]);
        }
        System.out.println();
    }

    public void printRowEchelonForm2(float[][] A, float[] B) {
        int N = B.length;
        System.out.println("\nRow Echelon form : ");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%.3f ", A[i][j]);
            }
            System.out.printf("| %.3f\n", B[i]);
        }
        System.out.println();
    }
}
