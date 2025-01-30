/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package toroidaldiffusion;

import lphy.core.logger.LoggerUtils;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.ejml.simple.SimpleMatrix;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * @author Michael Golden
 */
public class WrappedBivariateDiffusion {

    int lk = 1;
    SimpleMatrix x = new SimpleMatrix(4, 1);
    double t = 1.0;
    public SimpleMatrix mu = new SimpleMatrix(2, 1);
    SimpleMatrix alpha = new SimpleMatrix(3, 1);
    SimpleMatrix sigma = new SimpleMatrix(2, 1);
    int maxK = 1;
    double etrunc = 100;

    SimpleMatrix vstores;
    SimpleMatrix vstoret;
    SimpleMatrix weightswindsinitial;
    SimpleMatrix logweightswindsinitial;
    SimpleMatrix A = new SimpleMatrix(2, 2);
    SimpleMatrix Sigmamat = new SimpleMatrix(2, 2);
    SimpleMatrix oneoverSigmamat = new SimpleMatrix(2, 2);
    SimpleMatrix invGammat = new SimpleMatrix(2, 2);

    SimpleMatrix twokpi;
    SimpleMatrix twokepivec = new SimpleMatrix(2, 1);
    double penalty;
    SimpleMatrix invSigmaA;
    SimpleMatrix xmuinvSigmaA;
    double lognormconstSigmaA;
    SimpleMatrix ASigma;
    SimpleMatrix ASigmaA;
    double lognormconstGammat;
    SimpleMatrix Gammat;
    SimpleMatrix ExptA;
    double s;
    double q;

    public WrappedBivariateDiffusion() {
        setParameters(1.0, 0.0, 0.0, 0.2, 0.2, 0.05, 1.0, 1.0);
    }

    public void setParameters(double t, double muphi, double mupsi, double alphaphi, double alphapsi, double alpharho, double sigmaphi, double sigmapsi) {
        mu.set(0, 0, muphi);
        mu.set(1, 0, mupsi);
        alpha.set(0, 0, alphaphi);
        alpha.set(1, 0, alphapsi);
        alpha.set(2, 0, alpharho);
        sigma.set(0, 0, sigmaphi);
        sigma.set(1, 0, sigmapsi);
        double[] muarr = {muphi, mupsi};
        double[] alphaarr = {alphaphi, alphapsi, alpharho};
        double[] sigmaarr = {sigmaphi, sigmapsi};
        setParameters(muarr, alphaarr, sigmaarr);
        setParameters(t);
    }

    public void setParameters(Double[] muarr, Double[] alphaarr, Double[] sigmaarr) {
        setParameters(new double[] {muarr[0], muarr[1]},
                new double[] {alphaarr[0], alphaarr[1], alphaarr[2]},
                new double[] {sigmaarr[0], sigmaarr[1]});
    }


    public void setParameters(double[] muarr, double[] alphaarr, double[] sigmaarr) {
        mu.set(0, 0, muarr[0]);
        mu.set(1, 0, muarr[1]);
        alpha.set(0, 0, alphaarr[0]);
        alpha.set(1, 0, alphaarr[1]);
        alpha.set(2, 0, alphaarr[2]);
        sigma.set(0, 0, sigmaarr[0]);
        sigma.set(1, 0, sigmaarr[1]);

        double quo = Math.sqrt(sigma.get(0, 0) / sigma.get(1, 0));
        A.set(0, 0, alpha.get(0, 0));
        A.set(1, 1, alpha.get(1, 0));
        A.set(0, 1, alpha.get(2, 0) * quo);
        A.set(1, 0, alpha.get(2, 0) / quo);

        double testalpha = alpha.get(0, 0) * alpha.get(1, 0) - alpha.get(2, 0) * alpha.get(2, 0);
        penalty = 0.0;
        if (testalpha <= 0.0) {
            penalty = -testalpha * 100000.0 + 100.0;
            alpha.set(2, 0, Math.signum(alpha.get(2)) * Math.sqrt(alpha.get(0) * alpha.get(1)) * 0.9999);
            A.set(0, 1, alpha.get(2, 0) * quo);
            A.set(1, 0, alpha.get(2, 0) / quo);
        }

        Sigmamat.set(0, 0, sigma.get(0, 0));
        Sigmamat.set(1, 1, sigma.get(1, 0));
        lk = 2 * maxK + 1;
        ArrayList<Double> twokpiarr = linspace(-2.0 * maxK * Math.PI, 2.0 * maxK * Math.PI, lk);
        twokpi = new SimpleMatrix(lk, 1);
        for (int i = 0; i < lk; i++) {
            twokpi.set(i, 0, twokpiarr.get(i));
        }

        LoggerUtils.log.fine("2kπ");
//        twokpi.print(); // if debug

        oneoverSigmamat.set(0, 0, 1.0 / sigma.get(0, 0));
        oneoverSigmamat.set(1, 1, 1.0 / sigma.get(1, 0));
        invSigmaA = this.oneoverSigmamat.mult(A).scale(2.0);
        lognormconstSigmaA = -Math.log(2.0 * Math.PI) + Math.log(invSigmaA.determinant()) / 2.0;

        ASigma = A.mult(Sigmamat);
        ASigmaA = ASigma.mult(A.transpose());
        ASigma = ASigma.plus(ASigma.transpose());


        s = A.trace() / 2.0;
        q = Math.sqrt(Math.abs(A.minus(SimpleMatrix.identity(2).scale(s)).determinant()));

        weightswindsinitial = new SimpleMatrix(lk * lk, 1);
        vstores = new SimpleMatrix(lk * lk, 1);
        for (int wek1 = 0; wek1 < lk; wek1++) {
            for (int wek2 = 0; wek2 < lk; wek2++) {
                double v = ((invSigmaA.get(0, 0) * twokpi.get(wek1) + invSigmaA.get(0, 1) * twokpi.get(wek2)) * twokpi.get(wek1, 0) + (invSigmaA.get(1, 0) * twokpi.get(wek1, 0) + invSigmaA.get(1, 1) * twokpi.get(wek2, 0)) * twokpi.get(wek2)) / 2.0;
                vstores.set((wek1) * lk + wek2, 0, v);
            }
        }
        LoggerUtils.log.fine("vstores");
//        vstores.print(); // if debug

        double tempt = this.t;
        this.t = -2.0; // use this to make this.t != t, and then re-compute ExptA and invGammat, ...
        setParameters(tempt);
    }

    double log2pi = Math.log(2.0 * Math.PI);

    public void setParameters(double t) {
        if (this.t != t) {
            this.t = t;

            double q2 = q * q;
            double s2 = s * s;
            double est = Math.exp(s * t);
            double e2st = est * est;
            double inve2st = 1.0 / e2st;
            double c2st = Math.exp(2.0 * q * t);
            double s2st = (c2st - 1.0 / c2st) / 2.0;
            c2st = (c2st + 1.0 / c2st) / 2.0;

            double cte = inve2st / (4.0 * q2 * s * (s2 - q2));
            double integral1 = cte * (-s2 * (3.0 * q2 + s2) * c2st - q * s * (q2 + 3.0 * s2) * s2st - q2 * (q2 - 5.0 * s2) * e2st + (q2 - s2) * (q2 - s2));
            double integral2 = cte * s * ((q2 + s2) * c2st + 2.0 * q * s * s2st - 2.0 * q2 * e2st + q2 - s2);
            double integral3 = cte * (-s * (s * c2st + q * s2st) + (e2st - 1.0) * q2 + s2);

            Gammat = Sigmamat.scale(integral1).plus(ASigma.scale(integral2)).plus(ASigmaA.scale(integral3));

            double eqt = Math.exp(q * t);
            double cqt = (eqt + 1.0 / eqt) / 2.0;
            double sqt = (eqt - 1.0 / eqt) / 2.0;
            ExptA = SimpleMatrix.identity(2).scale(cqt + s * sqt / q).minus(A.scale((sqt / q))).scale(1.0 / est);

            double z = 1.0 / (Gammat.get(0, 0) * Gammat.get(1, 1) - Gammat.get(0, 1) * Gammat.get(1, 0));
            invGammat.set(0, 0, z * Gammat.get(1, 1));
            invGammat.set(0, 1, -z * Gammat.get(0, 1));
            invGammat.set(1, 0, -z * Gammat.get(1, 0));
            invGammat.set(1, 1, z * Gammat.get(0, 0));
            lognormconstGammat = -log2pi + Math.log(invGammat.determinant()) / 2.0;

            weightswindsinitial = new SimpleMatrix(lk * lk, 1);
            vstoret = new SimpleMatrix(lk * lk, 1);
            for (int wek1 = 0; wek1 < lk; wek1++) {
                for (int wek2 = 0; wek2 < lk; wek2++) {
                    // ((node.invGammat[1,1]*node.twokpi[wek1] + node.invGammat[1,2]*node.twokpi[wek2])*node.twokpi[wek1] + (node.invGammat[2,1]*node.twokpi[wek1] + node.invGammat[2,2]*node.twokpi[wek2])*node.twokpi[wek2])/2.0
                    double v = ((invGammat.get(0, 0) * twokpi.get(wek1) + invGammat.get(0, 1) * twokpi.get(wek2)) * twokpi.get(wek1, 0) + (invGammat.get(1, 0) * twokpi.get(wek1, 0) + invGammat.get(1, 1) * twokpi.get(wek2, 0)) * twokpi.get(wek2)) / 2.0;

                    vstoret.set((wek1) * lk + wek2, 0, v);
                }
            }
        }
    }

    // rejection sampling (only work for 1 or 2 dim) from TPD
    public double[][] sampleByRejection(double phi0, double psi0, int nsamples) {

        double[][] samples = new double[nsamples][2];
        int count = 0;
        int rejection = 0;
        WrappedBivariateDiffusionTransitionProbabilityFunction func = new WrappedBivariateDiffusionTransitionProbabilityFunction(this, phi0, psi0);

        BOBYQAOptimizer optimizer = new BOBYQAOptimizer(5);

        PointValuePair pvp = optimizer.optimize(
                new MaxEval(200),
                GoalType.MAXIMIZE,
                new InitialGuess(new double[]{phi0, psi0}),
                new ObjectiveFunction(func),
                new SimpleBounds(new double[]{0, 0}, new double[]{2.0 * Math.PI, 2.0 * Math.PI}));

        //double[] p = pvp.getPoint();

        double maxP = Math.exp(pvp.getValue()) * 1.01;

        RandomGenerator random = RandomUtils.getRandom();
        while (count < samples.length) {

            double ph = random.nextDouble() * 2.0 * Math.PI;
            double ps = random.nextDouble() * 2.0 * Math.PI;
            double density = random.nextDouble() * maxP;

            if (density < Math.exp(loglikwndtpd(phi0, psi0, ph, ps))) {
                samples[count][0] = ph;
                samples[count][1] = ps;
                count += 1;
            } else {
                rejection += 1;
            }
        }
        //System.out.println(rejection + " rejections to sample " + count + " points.");
        return samples;
    }

    // WN wrapped normal, TPD transition probability density
    public double loglikwndtpd(double phi0, double psi0, double phit, double psit) {
        if (penalty > 0.0) {
            return Double.NEGATIVE_INFINITY;
        }

        x.set(0, phi0);
        x.set(1, psi0);
        x.set(2, phit);
        x.set(3, psit);

        SimpleMatrix xmu = new SimpleMatrix(2, 1);
        xmu.set(0, 0, x.get(0, 0) - mu.get(0, 0));
        xmu.set(1, 0, x.get(1, 0) - mu.get(1, 0));

        xmuinvSigmaA = invSigmaA.mult(xmu);
        double xmuinvSigmaAxmudivtwo = (xmuinvSigmaA.get(0) * xmu.get(0) + xmuinvSigmaA.get(1) * xmu.get(1)) / 2.0;

        double logtpdfinal = Double.NEGATIVE_INFINITY;
        SimpleMatrix x0 = new SimpleMatrix(2, 1);
        x0.set(0, 0, x.get(0, 0));
        x0.set(1, 0, x.get(1, 0));
        logweightswindsinitial = new SimpleMatrix(lk * lk, 1);
        for (int wek1 = 0; wek1 < lk; wek1++) {
            twokepivec.set(0, 0, twokpi.get(wek1, 0));
            for (int wek2 = 0; wek2 < lk; wek2++) {
                int index = wek1 * lk + wek2;

                double exponent = xmuinvSigmaAxmudivtwo + (xmuinvSigmaA.get(0, 0) * twokpi.get(wek1, 0) + xmuinvSigmaA.get(1, 0) * twokpi.get(wek2, 0) + vstores.get(index, 0) - lognormconstSigmaA);

                if (exponent <= etrunc) {
                    logweightswindsinitial.set(index, 0, -exponent);
                } else {
                    logweightswindsinitial.set(index, 0, Double.NEGATIVE_INFINITY);
                }

                if (logweightswindsinitial.get(index, 0) > Double.NEGATIVE_INFINITY) {
                    twokepivec.set(1, 0, twokpi.get(wek2, 0));
                    SimpleMatrix mut = mu.plus(ExptA.mult(x0.plus(twokepivec).minus(mu)));
                    SimpleMatrix xmut = new SimpleMatrix(2, 1);
                    xmut.set(0, 0, x.get(2, 0) - mut.get(0, 0));
                    xmut.set(1, 0, x.get(3, 0) - mut.get(1, 0));

                    SimpleMatrix xmutinvGammat = invGammat.mult(xmut);
                    double xmutinvGammatxmutdiv2 = (xmutinvGammat.get(0, 0) * xmut.get(0, 0) + xmutinvGammat.get(1, 0) * xmut.get(1, 0)) / 2.0;

                    double logtpdintermediate = Double.NEGATIVE_INFINITY;
                    for (int wak1 = 0; wak1 < lk; wak1++) {
                        for (int wak2 = 0; wak2 < lk; wak2++) {
                            exponent = xmutinvGammatxmutdiv2 + (xmutinvGammat.get(0, 0) * twokpi.get(wak1, 0) + xmutinvGammat.get(1, 0) * twokpi.get(wak2, 0)) + vstoret.get(wak1 * lk + wak2) - lognormconstGammat;

                            logtpdintermediate = logsumexp(logtpdintermediate, -exponent);
                        }
                    }

                    logtpdfinal = logsumexp(logtpdfinal, logtpdintermediate + logweightswindsinitial.get(index, 0));
                }

            }
        }

        double ll = logtpdfinal - penalty;
        if (Double.isNaN(ll) || ll == Double.NEGATIVE_INFINITY) {
            return -1e10;
        } else {
            return ll;
        }
    }

    // pdf of stationary dist
    public double loglikwndstat(double phi, double psi) {
        if (penalty > 0.0) {
            return -1.0e10;
        }

        x.set(0, 0, phi);
        x.set(1, 0, psi);

        SimpleMatrix xmu = new SimpleMatrix(2, 1);
        xmu.set(0, 0, x.get(0, 0) - mu.get(0, 0));
        xmu.set(1, 0, x.get(1, 0) - mu.get(1, 0));

        xmuinvSigmaA = invSigmaA.mult(xmu);
        double xmuinvSigmaAxmudivtwo = (xmuinvSigmaA.get(0) * xmu.get(0) + xmuinvSigmaA.get(1) * xmu.get(1)) / 2.0;

        double logweightswindsinitialsum = Double.NEGATIVE_INFINITY;
        for (int wek1 = 0; wek1 < lk; wek1++) {
            for (int wek2 = 0; wek2 < lk; wek2++) {
                double exponent = xmuinvSigmaAxmudivtwo + (xmuinvSigmaA.get(0, 0) * twokpi.get(wek1, 0) + xmuinvSigmaA.get(1, 0) * twokpi.get(wek2, 0)) + vstores.get(wek1 * lk + wek2, 0) - lognormconstSigmaA;
                logweightswindsinitialsum = logsumexp(logweightswindsinitialsum, -exponent);
            }
        }

        double ll = logweightswindsinitialsum - penalty;
        if (Double.isNaN(ll) || ll == Double.NEGATIVE_INFINITY) {
            return -1e10;
        } else {
            return ll;
        }
    }

    // TODO check
    // rejection sampling (only work for 1 or 2 dim) from Stationary Density
    public double[][] rejectSampleStationaryDensity(double phiMu, double psiMu, int nsamples) {

        double[][] samples = new double[nsamples][2];
        int count = 0;
        int rejection = 0;
        WrappedBivariateDiffusionStationaryDensityFunction func = new WrappedBivariateDiffusionStationaryDensityFunction(this);

        BOBYQAOptimizer optimizer = new BOBYQAOptimizer(5);

        PointValuePair pvp = optimizer.optimize(
                new MaxEval(200),
                GoalType.MAXIMIZE,
                new InitialGuess(new double[]{phiMu, psiMu}),
                new ObjectiveFunction(func),
                new SimpleBounds(new double[]{0, 0}, new double[]{2.0 * Math.PI, 2.0 * Math.PI}));

        //double[] p = pvp.getPoint();

        double maxP = Math.exp(pvp.getValue()) * 1.01;

        RandomGenerator random = RandomUtils.getRandom();
        while (count < samples.length) {

            double ph = random.nextDouble() * 2.0 * Math.PI;
            double ps = random.nextDouble() * 2.0 * Math.PI;
            double density = random.nextDouble() * maxP;

            if (density < Math.exp(loglikwndstat(ph, ps))) {
                samples[count][0] = ph;
                samples[count][1] = ps;
                count += 1;
            } else {
                rejection += 1;
            }
        }
        //System.out.println(rejection + " rejections to sample " + count + " points.");
        return samples;
    }


    public static ArrayList<Double> linspace(double start, double stop, int n) {
        ArrayList<Double> result = new ArrayList<>();

        double step = (stop - start) / (n - 1);

        for (int i = 0; i <= n - 2; i++) {
            result.add(start + (i * step));
        }
        result.add(stop);

        return result;
    }
    
    /*
    public SimpleMatrix sampstat(Random rng)
    {
        SimpleMatrix invASigma = A.invert().mult(sigma.diag()).scale(0.5);
        SimpleMatrix x = new SimpleMatrix(2,1);
        x.set(0,0, rng.nextGaussian());
        x.set(1,0, rng.nextGaussian());
        CholeskyDecomposition_F64<DMatrixRMaj> chol = DecompositionFactory_DDRM.chol(invASigma.numRows(),true);
        x = x.transpose().plus(SimpleMatrix.wrap(chol.getT(null)));
        x = x.plus(mu.transpose());
        
        SimpleMatrix temp = x.plus(Math.PI).scale(1.0/(2.0*Math.PI));
        temp.set(0,0, Math.floor(temp.get(0,0)));
        temp.set(1,0, Math.floor(temp.get(1,0)));
        x = x.minus(temp);
        return x;
    }*/


    public static double logsumexp(double a, double b) {
        if (b == Double.NEGATIVE_INFINITY) {
            return a;
        } else if (a == Double.NEGATIVE_INFINITY) {
            return b;
        } else if (a < b) {
            return b + Math.log1p(Math.exp(a - b));
        } else {
            return a + Math.log1p(Math.exp(b - a));
        }
    }

    public double[][] simulatePath(double phi0, double psi0, double dt, int steps) {
        double[][] path = new double[steps + 1][2];

        path[0][0] = phi0;
        path[0][1] = psi0;

        setParameters(dt);
        for (int i = 1; i < path.length; i++) {
            double[][] point = sampleByRejection(phi0, psi0, 1);
            phi0 = point[0][0];
            psi0 = point[0][1];

            path[i][0] = phi0;
            path[i][1] = psi0;
        }
        return path;
    }

    // return mu_t[2], the mean of WN after time t from [phi0, psi0]
    public double[] wntpdMu(double phi0, double psi0, double t) {
        setParameters(t);
        SimpleMatrix x0 = new SimpleMatrix(2, 1);
        x0.set(0, 0, phi0);
        x0.set(1, 0, psi0);
        SimpleMatrix mut = mu.plus(ExptA.mult(x0.plus(twokepivec).minus(mu)));

        // wrapped
        double wrapPhit = ToroidalUtils.wrapToMaxAngle(mut.get(0, 0), WrappedNormalConst.MAX_ANGLE_VALUE);
        double wrapPsit = ToroidalUtils.wrapToMaxAngle(mut.get(1, 0), WrappedNormalConst.MAX_ANGLE_VALUE);
        return new double[]{wrapPhit, wrapPsit};
    }


    // test how the params affect log likelihood and stationary dist
    public static void main(String[] args) throws IOException {
        WrappedBivariateDiffusion diff = new WrappedBivariateDiffusion();
        double[] muarr = {Math.PI, Math.PI * 0.5}; // mean of the diffusion
        double[] sigmaarr = {1, 1}; // variance term
        double[] alphaarr = {1.0, 1.0, 0.7}; // drift term
        diff.setParameters(muarr, alphaarr, sigmaarr); // set the diffusion parameters
        System.out.println(diff.loglikwndstat(0.0, 0.0)); // calculate the stationary density of the point (0.0, 0.0)

        final double MaxDegrees = 2 * Math.PI;

        String filename = "parameters.txt";
        PrintWriter writer = new PrintWriter(new FileWriter(filename));
        writer.println("mu1\tmu2\tsigma1\tsigma2\talpha1\talpha2\talpha3\tlogP(0,0)");
        writer.println(muarr[0] + "\t" + muarr[1] + "\t" + sigmaarr[0] + "\t" + sigmaarr[1] + "\t" +
                alphaarr[0] + "\t" + alphaarr[1] + "\t" + alphaarr[2] + "\t" +
                diff.loglikwndstat(0.0, 0.0));
        writer.flush();
        writer.close();


        filename = "wrappedNormalStationary.txt";
        int gridSize = 10;//200;
        int maxTimeInterval = 50;

        writer = new PrintWriter(new FileWriter(filename));
        writer.println("time\tphi\tpsi\tlogP\tdensity");
        for (int i = 0; i < gridSize; i++) {
            double phi = (i + 0.5) * MaxDegrees / (double) gridSize;
            for (int j = 0; j < gridSize; j++) {
                for (int t = 0; t <= 2; t++) {
                    // Stationary dist should not change for diff time
                    double time = t * 0.5;
                    diff.setParameters(time); // set the time parameter

                    double psi = (j + 0.5) * MaxDegrees / (double) gridSize;
                    double logP = diff.loglikwndstat(phi, psi);
                    writer.println(time + "\t" + phi + "\t" + psi + "\t" + logP + "\t" + Math.exp(logP)); // calculate the transition density of the point (0.0, 0.0) transitioning to (1.0, 1.0) in time t=1.0
                }
            }
        }
        writer.flush();
        writer.close();


        filename = "wrappedNormal.txt";
        double phi0 = Math.PI / 4;
        double psi0 = Math.PI;

        writer = new PrintWriter(new FileWriter(filename));
        writer.println("phit\tpsit\tlogP\tdensity");
        for (int i = 0; i < gridSize; i++) {
            double phit = (i + 0.5) * MaxDegrees / (double) gridSize;
            for (int j = 0; j < gridSize; j++) {
                double psit = (j + 0.5) * MaxDegrees / (double) gridSize;

                double logP = diff.loglikwndtpd(phi0, psi0, phit, psit);

                diff.setParameters(0.5); // set the time parameter
                writer.println(phit + "\t" + psit + "\t" + logP + "\t" + Math.exp(logP)); // calculate the transition density of the point (0.0, 0.0) transitioning to (1.0, 1.0) in time t=1.0

                //diff.setParameters(0.7); // change the time parameter
                //System.out.println(diff.loglikwndtpd(0.0, 0.0, 1.0, 1.0)); // calculate the transition density for the same points, but for a different time (t=0.7)
            }
        }
//        writer.println("time\tphi0\tpsi0\tphit\tpsit\tlogP\tdensity");
//        for (int i0 = 0; i0 < gridSize; i0++) {
//            double phi0 = (i0 + 0.5) * MaxDegrees / (double) gridSize;
//            for (int j0 = 0; j0 < gridSize; j0++) {
//                double psi0 = (j0 + 0.5) * MaxDegrees / (double) gridSize;
//                // after t
//                for (int i = 0; i < gridSize; i++) {
//                    double phit = (i + 0.5) * MaxDegrees / (double) gridSize;
//                    for (int j = 0; j < gridSize; j++) {
////                        for (int t = 0; t <= maxTimeInterval; t++) {
//                            // test time
//                            double time = 0.5;//t * 1.0 / maxTimeInterval;
//                            diff.setParameters(time); // set the time parameter
//
//                            double psit = (j + 0.5) * MaxDegrees / (double) gridSize;
//
//                            double logP = diff.loglikwndtpd(phi0, psi0, phit, psit);
////                diff.setParameters(0.5); // set the time parameter
//                            writer.println(time + "\t" + phi0 + "\t" + psi0 + "\t" + phit + "\t" + psit + "\t" + logP + "\t" + Math.exp(logP)); // calculate the transition density of the point (0.0, 0.0) transitioning to (1.0, 1.0) in time t=1.0
//
//                            //diff.setParameters(0.7); // change the time parameter
//                            //System.out.println(diff.loglikwndtpd(0.0, 0.0, 1.0, 1.0)); // calculate the transition density for the same points, but for a different time (t=0.7)
//                        }
//                    }
//                }
//            }
//        }
        writer.flush();
        writer.close();

        double[][] samples = diff.sampleByRejection(phi0, psi0, 1000);
        writer = new PrintWriter(new FileWriter("wrappedNormalSample.txt"));
        writer.println("phit\tpsit");
        for (int i = 0; i < samples.length; i++) {
            writer.println(samples[i][0] + "\t" + samples[i][1]);
        }
        writer.flush();
        writer.close();

        double[][] path = diff.simulatePath(0, Math.PI, 0.02, 1000);
        writer = new PrintWriter(new FileWriter("wrappedNormalPath.txt"));
        writer.println("phi\tpsi");
        for (int i = 0; i < path.length; i++) {
            writer.println(path[i][0] + "\t" + path[i][1]);
        }
        writer.flush();
        writer.close();

        /** after time t **/

//        filename = "muT.txt";
//        int N = 10000; // 10k
//
//        writer = new PrintWriter(new FileWriter(filename));
//        writer.println("time\tmu_t_phi\tmu_t_psi");
//        for (int t = 0; t < N; t++) {
//            // time
//            double time = t / 100.0;
//            diff.setParameters(time); // set the time parameter
//
//            double[] muT = diff.wntpdMu(0, 0, t);
//            //TODO std?
//            writer.println(time + "\t" + muT[0] + "\t" + muT[1]); // calculate the transition density of the point (0.0, 0.0) transitioning to (1.0, 1.0) in time t=1.0
//        }
//
//        writer.flush();
//        writer.close();

        /**
         * store trajectory at fixed time interval
         */

        filename = "muT.txt";
        int steps = 1000;  // Number of steps
        double[][] trajectory = new double[steps + 1][2];  // +1 to include initial point

        writer = new PrintWriter(new FileWriter(filename));
        writer.println("time\tmu_t_phi\tmu_t_psi");

// Store initial point
        final double dt = 0.5;  // Fixed time interval
        double phi = 0;
        double psi = Math.PI;
        trajectory[0][0] = phi;
        trajectory[0][1] = psi;
        writer.println(0 + "\t" + trajectory[0][0] + "\t" + trajectory[0][1]);

// Calculate trajectory at fixed time intervals
        for (int i = 1; i <= steps; i++) {
            double[] muT = diff.wntpdMu(trajectory[i-1][0], trajectory[i-1][1], dt);

            trajectory[i][0] = muT[0];
            trajectory[i][1] = muT[1];

            writer.println(i + "\t" + muT[0] + "\t" + muT[1]);
        }

//// Calculate trajectory at fixed time intervals
//        for (int i = 1; i <= steps; i++) {
//            time = i * dt;
//            diff.setParameters(time);
//            muT = diff.wntpdMu(0, Math.PI, time);
//
//            trajectory[i][0] = muT[0];
//            trajectory[i][1] = muT[1];
//
//            writer.println(time + "\t" + muT[0] + "\t" + muT[1]);
//        }

        writer.flush();
        writer.close();

    }

}
