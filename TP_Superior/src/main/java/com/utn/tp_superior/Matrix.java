package com.utn.tp_superior;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Juli
 */
public class Matrix {

    private int nrows;
    private int ncols;
    private double[][] data;

    public Matrix(double[][] dat) {
        this.data = dat;
        this.nrows = dat.length;
        this.ncols = dat[0].length;
    }



    public int getNrows() {
        return nrows;
    }

    public int getNcols() {
        return ncols;
    }

    public double getValueAt(int i, int j) {
        return data[i][j];
    }

    public void setValueAt(int i, int j, double data) {
        this.data[i][j] = data;
    }

    public Matrix(int nrow, int ncol) {
        this.nrows = nrow;
        this.ncols = ncol;
        data = new double[nrow][ncol];
    }

    public static Matrix transpose(Matrix matrix) {
        Matrix transposedMatrix = new Matrix(matrix.getNcols(), matrix.getNrows());
        for (int i = 0; i < matrix.getNrows(); i++) {
            for (int j = 0; j < matrix.getNcols(); j++) {
                transposedMatrix.setValueAt(j, i, matrix.getValueAt(i, j));
            }
        }
        return transposedMatrix;
    }

    public boolean isSquare() {
        return (this.nrows == this.ncols);
    }

    public static double changeSign(int i) {
        if (i % 2 == 0) {
            return 1;
        } else {
            return -1;
        }
    }

    public static Matrix createSubMatrix(Matrix matrix, int excluding_row, int excluding_col) {
        Matrix mat = new Matrix(matrix.getNrows() - 1, matrix.getNcols() - 1);
        int r = -1;
        for (int i = 0; i < matrix.getNrows(); i++) {
            if (i == excluding_row) {
                continue;
            }
            r++;
            int c = -1;
            for (int j = 0; j < matrix.getNcols(); j++) {
                if (j == excluding_col) {
                    continue;
                }
                mat.setValueAt(r, ++c, matrix.getValueAt(i, j));
            }
        }
        return mat;
    }

    public static double determinant(Matrix matrix) {
        if (!matrix.isSquare()) {
            System.out.println("matrix need to be square.");
        }
        if (matrix.getNrows() == 1) {
            return matrix.getValueAt(0, 0);
        }
        if (matrix.getNrows() == 2) {
            return (matrix.getValueAt(0, 0) * matrix.getValueAt(1, 1))
                    - (matrix.getValueAt(0, 1) * matrix.getValueAt(1, 0));
        }
        double sum = 0.0;
        for (int i = 0; i < matrix.getNcols(); i++) {
            sum += changeSign(i) * matrix.getValueAt(0, i) * determinant(createSubMatrix(matrix, 0, i));
        }
        return sum;
    }

    public static Matrix cofactor(Matrix matrix) {
        Matrix mat = new Matrix(matrix.getNrows(), matrix.getNcols());
        for (int i = 0; i < matrix.getNrows(); i++) {
            for (int j = 0; j < matrix.getNcols(); j++) {
                mat.setValueAt(i, j, changeSign(i) * changeSign(j)
                        * determinant(createSubMatrix(matrix, i, j)));
            }
        }
        return mat;
    }
    
    public Matrix multiplyByConstant(double c){
        for(int i = 0; i<this.getNcols(); i++){
            for(int j = 0; j<this.getNrows(); j++){
                this.data[i][j]=this.data[i][j]*c;
            }
        }
        return this;
    }
    
    public static Matrix inverse(Matrix matrix) {
            return (transpose(cofactor(matrix)).multiplyByConstant(1.0 / determinant(matrix)));
    }
}
