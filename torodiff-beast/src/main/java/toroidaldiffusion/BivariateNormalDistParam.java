package toroidaldiffusion;

public final class BivariateNormalDistParam {

    public final double meanPhi;
    public final double meanPsi;
    public final double varPhi;
    public final double varPsi;
    public final double covar;

    public BivariateNormalDistParam(double meanPhi, double meanPsi, double varPhi, double varPsi, double covar) {
        this.meanPhi = meanPhi;
        this.meanPsi = meanPsi;
        this.varPhi = varPhi;
        this.varPsi = varPsi;
        this.covar = covar;
    }
}
