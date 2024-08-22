package toroidaldiffusion;

import java.text.DecimalFormat;

public class Pair {

    Double one;
    Double two;

    public Pair(Double one, Double two) {
        this.one = one;
        this.two = two;
    }

    @Override
    public String toString() {
        DecimalFormat df = new DecimalFormat("#.###"); // Format to three decimal places
        return "(" + df.format(one) + ", " + df.format(two) + ")";
    }

    public Double getPhi() {
        // Round the 'phi' value (one) to three decimal places and return it
        DecimalFormat df = new DecimalFormat("#.###");
        return Double.valueOf(df.format(one));
    }

    public Double getPsi() {
        // Round the 'psi' value (two) to three decimal places and return it
        DecimalFormat df = new DecimalFormat("#.###");
        return Double.valueOf(df.format(two));
    }

    public Double[] getPair() {
        return new Double[]{one, two};
    }
}
