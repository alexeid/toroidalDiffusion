package toroidaldiffusion.studio;

import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRendererState;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.util.LineUtils;
import org.jfree.data.xy.XYDataset;

import java.awt.*;
import java.awt.geom.Rectangle2D;

public class WrappedXYLineRenderer extends XYLineAndShapeRenderer {
    private final double xMax;
    private final double yMax;

    private ToroidalPlot toroidalPlot;

    public WrappedXYLineRenderer(double xMax, double yMax, ToroidalPlot toroidalPlot) {
        super(true, true);
        this.xMax = xMax;
        this.yMax = yMax;
        this.toroidalPlot = toroidalPlot;
    }

    @Override
    protected void drawPrimaryLine(XYItemRendererState state,
                                   Graphics2D g2,
                                   XYPlot plot,
                                   XYDataset dataset,
                                   int pass,
                                   int series,
                                   int item,
                                   ValueAxis domainAxis,
                                   ValueAxis rangeAxis,
                                   Rectangle2D dataArea) {
        if (item == 0) {
            return;
        }

        // get the data point...
        double x1 = dataset.getXValue(series, item);
        double y1 = dataset.getYValue(series, item);
        if (Double.isNaN(y1) || Double.isNaN(x1)) {
            return;
        }

        double x0 = dataset.getXValue(series, item - 1);
        double y0 = dataset.getYValue(series, item - 1);
        if (Double.isNaN(y0) || Double.isNaN(x0)) {
            return;
        }

        RectangleEdge xAxisLocation = plot.getDomainAxisEdge();
        RectangleEdge yAxisLocation = plot.getRangeAxisEdge();

        double transX0 = domainAxis.valueToJava2D(x0, dataArea, xAxisLocation);
        double transY0 = rangeAxis.valueToJava2D(y0, dataArea, yAxisLocation);

        double transX1 = domainAxis.valueToJava2D(x1, dataArea, xAxisLocation);
        double transY1 = rangeAxis.valueToJava2D(y1, dataArea, yAxisLocation);

        // only draw if we have good values
        if (Double.isNaN(transX0) || Double.isNaN(transY0)
                || Double.isNaN(transX1) || Double.isNaN(transY1)) {
            return;
        }

        double euclideanDist = calculateDistance(x0, y0, x1, y1);
        double toroidalDistance = toroidalPlot.toroidalDistance(x0, y0, x1, y1);

        PlotOrientation orientation = plot.getOrientation();
        boolean visible;
        if (orientation == PlotOrientation.HORIZONTAL) {
            state.workingLine.setLine(transY0, transX0, transY1, transX1);
        }
        else if (orientation == PlotOrientation.VERTICAL) {
            state.workingLine.setLine(transX0, transY0, transX1, transY1);
        }
        visible = LineUtils.clipLine(state.workingLine, dataArea);
        if (visible) {
            if (Math.abs(toroidalDistance - euclideanDist) > 1e-8) {
                // draw wrapped line
// TODO               Line2D l = wrappedLine();
//                drawFirstPassShape(g2, pass, series, item, l);
            } else // draw normally
                drawFirstPassShape(g2, pass, series, item, state.workingLine);
        }
    }

    public static double calculateDistance(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }
}
