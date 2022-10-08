package toroidalDiffusion;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Line2D;

public class DihedralAngle3DViewer extends JComponent {

    Double[][] structure;
    double bondLength = 10;

    public DihedralAngle3DViewer(Double[][] structure){
        this.structure = structure;
    }

    public void paintComponent(Graphics g){
        Graphics2D g2D = (Graphics2D)g;
        double x = 0;
        double y = 0;
        for (Double[] anglePair : structure) {
            double phi = anglePair[0];
            double psi = anglePair[1];
            double xnew = x + Math.sin(phi) * bondLength;
            double ynew = y + Math.cos(phi) * bondLength;
            Line2D line = new Line2D.Double(x, y, xnew, ynew);
            g2D.draw(line);
        }
    }
}
