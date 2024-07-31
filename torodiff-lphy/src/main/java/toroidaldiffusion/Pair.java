package toroidaldiffusion;

public class Pair {

    Double one;
    Double two;

    public Pair(Double one, Double two) {
        this.one = one;
        this.two = two;
    }

    public Double[] getPair() {
        return new Double[]{one, two};
    }
}
