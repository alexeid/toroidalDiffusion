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

        int width = getWidth();
        int height = getHeight();

        Graphics2D g2D = (Graphics2D)g;
        double x = width/2.0;
        double y = height/2.0;
        for (Double[] anglePair : structure) {
            double phi = anglePair[0];
            double psi = anglePair[1];
            double xnew = x + Math.sin(phi) * bondLength;
            double ynew = y + Math.cos(phi) * bondLength;
            Line2D line = new Line2D.Double(x, y, xnew, ynew);
            g2D.draw(line);
            x = xnew;
            y = ynew;
        }
    }

    public static void main(String[] args ) {

        Double[][] structure = new Double[20][2];

        for (int i = 0; i < structure.length; i++) {
            structure[i][0] = Math.random()*Math.PI-Math.PI/2.0;
            structure[i][1] = Math.random()*Math.PI-Math.PI/2.0;
        }

        DihedralAngle3DViewer viewer = new DihedralAngle3DViewer(structure);

        JFrame frame = new JFrame("DihedralAngle3DViewer");
        frame.setSize(1000,800);
        frame.getContentPane().add(viewer);
        frame.setVisible(true);
    }
}
