/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.utn.tp_superior;

import java.awt.Color;
import static java.lang.Boolean.FALSE;
import static java.lang.Boolean.TRUE;
import javafx.scene.chart.NumberAxis;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author Julian
 */
public class Scatter extends JFrame {

    public Scatter(String title, double arreglo[][]) {
        super(title);

        // Create dataset
        XYDataset dataset = createDataset(arreglo);

        // Create chart
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Comparación norma errores",
                "Norma error (escala log10)",
                "Tamaño matriz (n)",
                dataset,
                PlotOrientation.HORIZONTAL,
                TRUE,
                FALSE,
                FALSE);

        //Changes background color
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(new Color(255, 228, 196));
        plot.setRangeAxisLocation(AxisLocation.BOTTOM_OR_LEFT);

        // Create Panel
        ChartPanel panel = new ChartPanel(chart);
        setContentPane(panel);
    }

    private XYDataset createDataset(double arreglo[][]) {
        XYSeriesCollection dataset = new XYSeriesCollection();

        XYSeries series1 = new XYSeries("Gauss");
        for (int i = 0; i < arreglo[0].length; i++) {
            series1.add(arreglo[0][i], i + 1);
        }

        XYSeries series2 = new XYSeries("Gauss (simple precisión");
        for (int i = 0; i < arreglo[1].length; i++) {
            series2.add(arreglo[1][i], i + 1);
        }
        
        XYSeries series3 = new XYSeries("Gauss-Seidel");
        for (int i = 0; i < arreglo[2].length; i++) {
            series3.add(arreglo[2][i], i + 1);
        }

        dataset.addSeries(series1);
        dataset.addSeries(series2);
        dataset.addSeries(series3);

        return dataset;
    }
    
    public Scatter(String title, double arregloR[][], double[][] arregloT) {
        super(title);

        // Create dataset
        XYDataset dataset = createDataset(arregloR, arregloT);

        // Create chart
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Comparación de tiempo en cada iteracion",
                "Tamaño de la matriz",
                "Tiempo de cálculo (ms)",
                dataset,
                PlotOrientation.HORIZONTAL,
                TRUE,
                FALSE,
                FALSE);

        //Changes background color
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(new Color(255, 228, 196));
        plot.setRangeAxisLocation(AxisLocation.BOTTOM_OR_LEFT);

        // Create Panel
        ChartPanel panel = new ChartPanel(chart);
        setContentPane(panel);
    }

    private XYDataset createDataset(double arregloR[][], double[][] arregloT) {
        XYSeriesCollection dataset = new XYSeriesCollection();

        XYSeries series2 = new XYSeries("Gauss-Seidel");
        for (int i = 0; i < arregloR[1].length; i++) {
            series2.add(i+1, arregloT[1][i]);
        }
        
        XYSeries series1 = new XYSeries("Gauss");
        for (int i = 0; i < arregloR[1].length; i++) {
            series1.add(i+1, arregloT[0][i]);
        }
        
        XYSeries series3 = new XYSeries("Gradientes conjugados");
        for (int i = 0; i < arregloR[3].length; i++) {
            series3.add(i+1, arregloT[2][i]);
        }

        dataset.addSeries(series1);
        dataset.addSeries(series2);
        dataset.addSeries(series3);

        return dataset;
    }

}
