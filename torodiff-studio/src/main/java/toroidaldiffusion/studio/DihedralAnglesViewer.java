package toroidaldiffusion.studio;

import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;
import toroidaldiffusion.DihedralAngleAlignment;

import javax.swing.*;

public class DihedralAnglesViewer implements Viewer {

    /**
     * Required by ServiceLoader.
     */
    public DihedralAnglesViewer() {
    }

    @Override
    public boolean match(Object value) {
        return value instanceof DihedralAngleAlignment ||
                (value instanceof Value && ((Value) value).value() instanceof DihedralAngleAlignment);
    }

    @Override
//    public JComponent getViewer(Object value) {
//        if (value instanceof DihedralAngleAlignment alignment) {
//            return new JLabel(alignmnetToString(alignment));
//        } else if (value instanceof Value v && v.value() instanceof DihedralAngleAlignment alignment) {
//            return new JLabel(alignmnetToString(alignment));
//        }
//        return new JLabel("Value<DihedralAngleAlignment> TODO : " + ((Value) value).value().toString());
//    }

    public JComponent getViewer(Object value) {
        String html;
        if (value instanceof DihedralAngleAlignment alignment) {
            html = alignment.toHTML();
        } else if (value instanceof Value<?> v && v.value() instanceof DihedralAngleAlignment alignment) {
            html = alignment.toHTML();
        } else {
            html = "<html>Invalid value provided.</html>";
        }
        return new JLabel(html);
    }


    @Override
    public String toString() {
        return "Dihedral Angles Viewer";
    }

    private String alignmnetToString(DihedralAngleAlignment alignment) {
        return alignment.toJSON();
    }

    private String alignmentToHTML(DihedralAngleAlignment alignment) {
        return alignment.toHTML();
    }
}



