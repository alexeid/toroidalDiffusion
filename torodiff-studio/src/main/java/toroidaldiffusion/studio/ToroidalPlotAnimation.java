package toroidaldiffusion.studio;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import toroidaldiffusion.WrappedBivariateDiffusion;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;


public class ToroidalPlotAnimation {

    final static double[] muarr = {Math.PI * 0.65, Math.PI * 0.8}; // mean of the diffusion
    final static double[] sigmaarr = {0.5, 0.75}; // variance term
    final static double[] alphaarr = {1.0, 1.0, 0.5}; // drift term

    static XYSeries createDataSet() {
        WrappedBivariateDiffusion diff = new WrappedBivariateDiffusion();
        diff.setParameters(muarr, alphaarr, sigmaarr); // set the diffusion parameters

        ToroidalPlot toroidalPlot = new ToroidalPlot(diff, 1000);

        // Dataset for path points
        Point2D[] path = toroidalPlot.path;
        // cancel auto sort
        XYSeries pathSeries = new XYSeries("Path", false);
        for (Point2D point2D : path) {
            pathSeries.add(point2D.getX(), point2D.getY());
        }

        return pathSeries;
    }

    public static void main(String[] args) {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        // Dataset for path points
        XYSeries pathSeries = createDataSet();
        // cancel auto sort
        XYSeries animatedSeries = new XYSeries("Animated Path", false);

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(pathSeries); // this shows all paths after finish
        dataset.addSeries(animatedSeries);

        // Create a scatter plot
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Wrapped Normal Diffusion (" + dataset.getItemCount(0) + " States)",
                "Phi", "Psi",
                dataset
        );
        chart.removeLegend();

        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setRangeGridlinesVisible(false);
        // transparency
        plot.setForegroundAlpha(0.5f);

        // draw mean
        BasicStroke dashLine = new BasicStroke(
                2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
                1.0f, new float[] {6.0f, 6.0f}, 0.0f
        );

        XYLineAnnotation mean1 = new XYLineAnnotation(
                muarr[0], 0, muarr[0], Math.PI * 2, dashLine, Color.lightGray);
        XYLineAnnotation mean2 = new XYLineAnnotation(
                0, muarr[1], Math.PI * 2, muarr[1], dashLine, Color.lightGray);
        plot.addAnnotation(mean1);
        plot.addAnnotation(mean2);

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, true);

        /**
         * Two Series: 0 contains all paths after animation finishes, 1 is Animated path
         */

        renderer.setSeriesPaint(0, Color.white); // Original path in light gray
        renderer.setSeriesPaint(1, Color.RED); // Animated path in red

        renderer.setSeriesShapesVisible(0, false);
        renderer.setSeriesLinesVisible(1, true);
        renderer.setSeriesShapesVisible(1, true);

        // Thicker line
//        renderer.setSeriesStroke(1, new BasicStroke(2.0f));
        // Circle shape
        renderer.setSeriesShape(1, new Ellipse2D.Double(-3, -3, 6, 6));

        renderer.setDrawSeriesLineAsPath(false);

        plot.setRenderer(renderer);

        ChartPanel chartPanel = new ChartPanel(chart);

        frame.add(chartPanel);
        frame.pack();
        frame.setVisible(true);

        // Start Animation
        Timer timer = new Timer(100, e -> {
            int currentSize = animatedSeries.getItemCount();
            if (currentSize < pathSeries.getItemCount()) {
                animatedSeries.add(
                        pathSeries.getX(currentSize).doubleValue(),
                        pathSeries.getY(currentSize).doubleValue()
                );
                chartPanel.repaint(); // Trigger the chart to repaint
            } else {
                ((Timer) e.getSource()).stop();
            }
        });
        timer.start();

    }

}
