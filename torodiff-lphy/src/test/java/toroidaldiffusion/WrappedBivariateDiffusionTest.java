package toroidaldiffusion;

import org.junit.jupiter.api.Test;

class WrappedBivariateDiffusionTest {

    final double[] mu = new double[]{4.0509262973708315, 5.1236852445169685};
    final double[] sigma = new double[]{1.0298711760259887, 0.1271923963524306};
    final double[] alpha = new double[]{0.75761827331078, 0.8251637431436857, 0.14751750525982674};


    @Test
    void testLikelihood() {
        WrappedBivariateDiffusion diffusion = new WrappedBivariateDiffusion();

        diffusion.setParameters(mu, alpha, sigma);
        double t = 0.04184533374639288;
        diffusion.setParameters(t);

        // TODO: mu too high ?
        /**
         * org.apache.commons.math3.exception.MathIllegalStateException: trust region step has failed to reduce Q
         * 	at commons.math3@3.6.1/org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer.bobyqb(BOBYQAOptimizer.java:880)
         * 	at commons.math3@3.6.1/org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer.bobyqa(BOBYQAOptimizer.java:340)
         * 	at commons.math3@3.6.1/org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer.doOptimize(BOBYQAOptimizer.java:252)
         * 	at commons.math3@3.6.1/org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer.doOptimize(BOBYQAOptimizer.java:49)
         * 	at commons.math3@3.6.1/org.apache.commons.math3.optim.BaseOptimizer.optimize(BaseOptimizer.java:153)
         * 	at commons.math3@3.6.1/org.apache.commons.math3.optim.BaseMultivariateOptimizer.optimize(BaseMultivariateOptimizer.java:65)
         * 	at commons.math3@3.6.1/org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer.optimize(MultivariateOptimizer.java:63)
         * 	at toroidaldiffusion@0.0.1-SNAPSHOT/toroidaldiffusion.WrappedBivariateDiffusion.sampleByRejection(WrappedBivariateDiffusion.java:209)
         */
        double[][] samples = diffusion.sampleByRejection(3.0875357245670547, 3.2389546454267055, 1);

        //vquad represents the quadratic model prediction of the change in the objective function due to the step.
        // diff = f − f_opt - v_quad
        //If vquad is positive or zero, it indicates that the quadratic model failed to predict improvement.
    }
}