package toroidaldiffusion.operator;

import org.ejml.simple.SimpleMatrix;
import toroidaldiffusion.WrappedBivariateDiffusion;

public class CovarianceMatrix {

    private WrappedBivariateDiffusion diff;
//    private final double phi0;
//    private final double psi0;
    private final double t;
    private double var_phi;
    private double var_psi;
    private double covar;

    public CovarianceMatrix(WrappedBivariateDiffusion diff, double t) {
        this.diff = diff;
//        this.phi0 = phi0;
//        this.psi0 = psi0;
        this.t = t;

        computeBivariateNormal(diff, t);
    }

    private void computeBivariateNormal(WrappedBivariateDiffusion diff, double t) {
        // Update time-dependent parameters in the diffusion model
        diff.setParameters(t);

        SimpleMatrix Gammat = diff.getGammat();

        // Get covariance matrix at time t
        this.var_phi = Gammat.get(0, 0);
        this.var_psi = Gammat.get(1, 1);
        this.covar = Gammat.get(0, 1);
    }

//    public double getPhi0() {
//        return phi0;
//    }
//
//    public double getPsi0() {
//        return psi0;
//    }

    public double getVarPhi() {
        return var_phi;
    }

    public double getVarPsi() {
        return var_psi;
    }

    public double getCovar() {
        return covar;
    }
}
