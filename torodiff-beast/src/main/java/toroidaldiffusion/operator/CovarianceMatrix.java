package toroidaldiffusion.operator;

import org.ejml.simple.SimpleMatrix;
import toroidaldiffusion.WrappedBivariateDiffusion;

public class CovarianceMatrix {

    private final WrappedBivariateDiffusion diff;

    private final double t;
    private double varPhi;
    private double varPsi;
    private double covar;

    public CovarianceMatrix(WrappedBivariateDiffusion diff, double t) {
        this.diff = diff;
        this.t = t;

        getCovarianceMatrix(diff, t);
    }

    private void getCovarianceMatrix(WrappedBivariateDiffusion diff, double t) {
        // Update time-dependent parameters in the diffusion model
        diff.setParameters(t);

        //SimpleMatrix Gammat = diff.getGammat();
        SimpleMatrix Sigmamat = diff.getSigmamat();

//         Get covariance matrix at time t
        //        this.varPsi = Gammat.get(1, 1);
//        this.covar = Gammat.get(0, 1);
//        if (this.covar != Gammat.get(1, 0))
//            throw new IllegalArgumentException("Covariance in digonal must be same !");
//    }
        this.varPhi = Sigmamat.get(0, 0) * t;
        this.varPsi = Sigmamat.get(1, 1) * t;
        this.covar = 0.001;
    }

        public double getVarPhi () {
            return varPhi;
        }

        public double getVarPsi () {
            return varPsi;
        }

        public double getCovar () {
            return covar;
        }
    }
