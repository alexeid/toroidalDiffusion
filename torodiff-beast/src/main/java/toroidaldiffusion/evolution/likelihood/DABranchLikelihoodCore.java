package toroidaldiffusion.evolution.likelihood;

import toroidaldiffusion.WrappedBivariateDiffusion;

/**
 * data augmentation likelihood core based on a branch for multithreading.
 * the branch is defined above the selected child node.
 */
public class DABranchLikelihoodCore extends AbstrDALikelihoodCore {
    final private int branchNr; // the node Nr of child node below the branch

    final protected int nrOfSites;

    // to store branch likelihood calculation per site:
    // 1st dimension index is used to determine the current values or stored values,
    // 2nd dimension is to store likelihood at each site.
    protected double[][] siteLogLd;
//    protected double[] storedSiteLogLd;

    // caching mechanism : index (0 or 1) for 1st dimension
    protected int currentIndex = 0;
    protected int storedIndex = 0;

    final WrappedBivariateDiffusion diff;

    /**
     * no initialization, for calculating site likelihoods at the root.
     * @param branchNr       the node Nr of child node below the branch
     * @param nrOfSites      number of sites (codon)
     */
    public DABranchLikelihoodCore(int branchNr, int nrOfSites, WrappedBivariateDiffusion diff) {
        this.branchNr = branchNr;
        this.nrOfSites = nrOfSites; // TODO impl data range to multithreading by sites
        this.diff = diff;

        initialize();
    }


    /**
     * initializes states, likelihood arrays.
     */
    @Override
    protected void initialize() {
        // the branch above the node
        // 1st dimension is matrix index (current, stored),
        siteLogLd = new double[2][nrOfSites];
//        storedSiteLogLd = new double[nrOfSites];
        //TODO diff.setParameters(muarr, alphaarr, sigmaarr); ?
    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        storedIndex = currentIndex;
//        System.arraycopy(siteLogLd, 0, storedSiteLogLd, 0, nrOfSites);
    }

    /**
     * Store current state
     */
    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        int tmp = currentIndex;
        currentIndex = storedIndex;
        storedIndex = tmp;
//        System.arraycopy(storedSiteLogLd, 0, siteLogLd, 0, nrOfSites);
    }

    @Override
    public void unstore() {
//        currentBrLdIndex = storedBrLdIndex;
    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        siteLogLd = null;
//        storedSiteLogLd = null;
        currentIndex = -1;
        storedIndex = -1;
    }


    //============ branch likelihood ============

    //TODO need to cache per site to avoid recalculation, when only the sequence at a site is changed
//    public void calculateBranchLd(final double[] parentNodeValues, final double[] childNodeValues) {
//
//        //todo: bug: k+2, check the length and final index within the bound
//
//        assert parentNodeValues.length == nrOfSites * 2;
//        for (int k = 0; k < nrOfSites; k++) {
//            // pairs of values, dimension is 2 (angles) * N_sites
//
//            double phi0 = parentNodeValues[k];
//            double psi0 = parentNodeValues[k + 1];
//
//            double phit = childNodeValues[k];
//            double psit = childNodeValues[k + 1];
//
//            // diff.setParameters(muarr, alphaarr, sigmaarr), only once in the init method
//            branchLogLd[currentBrLdIndex][k] = diff.loglikwndtpd(phi0, psi0, phit, psit);
//
//            if (branchLogLd[currentBrLdIndex][k] == 0) {
//                throw new RuntimeException("\nBranch above node " + getBranchNr() + " likelihood = 0 !\n" +
//                        "At site " + k); //+ ", child node = " + childNode + ", parent node = " + parentNode);
//            }
//
//        } // end k  nrOfSites
//
//    }

    /**
     * set mu, sigma, alpha, before this method
     * @param parentNodeValues   pairs of angles at the parent node, length is nrOfSites * 2.
     * @param childNodeValues    pairs of angles at the child node, length is nrOfSites * 2.
     * @param branchTime         branch time (length) from parent to child
     */
    public void calculateBranchLd(final double[] parentNodeValues, final double[] childNodeValues, double branchTime) {

        assert parentNodeValues.length == nrOfSites * 2;
        assert childNodeValues.length == nrOfSites * 2;

        /**
         * switch the current index, before the likelihood calculation.
         * because store() makes storedIndex = currentIndex; when moving to the next MCMC state.
         */
        currentIndex = 1 - currentIndex; // stay 0 or 1

        //Require to set time before loglikwndtpd
        diff.setParameters(branchTime);
        // k is site index
        for (int k = 0; k < nrOfSites; k++) {
            // two angles, phi and psi, at the parent node and the kth site.
            double phi0 = parentNodeValues[k * 2];
            double psi0 = parentNodeValues[k * 2 + 1];
            // two angles, phi and psi, at the child node and the kth site.
            double phit = childNodeValues[k * 2];
            double psit = childNodeValues[k * 2 + 1];

            // diff.setParameters(muarr, alphaarr, sigmaarr), only once per tree likelihood calculation
            siteLogLd[currentIndex][k] = diff.loglikwndtpd(phi0, psi0, phit, psit);

            if (siteLogLd[currentIndex][k] == 0) {
                throw new RuntimeException("\nBranch above node " + getBranchNr() +
                        ", siteLogLd[" + k + "] = " + siteLogLd[currentIndex][k] + ", branchTime = " + branchTime +
                        "\nphi0 = " + phi0 + ", psi0 = " + psi0 + ", phit = " + phit + ", psit = " + psit);
                //+ ", child node = " + childNode + ", parent node = " + parentNode);
            }

        } // end k  nrOfSites

    }


    /**
     * Calculates log likelihood at this branch.
     * Multiple given site likelihoods, and then log the product.
     * The likelihood input has been integrated across categories.
     * It can be also used for the site likelihoods at the root.
     * @return           logged likelihood
     */
    public double calculateBranchLogLikelihood() {
        double logP = 0;
        // siteLd[].length == nrOfSites
        for (int k = 0; k < nrOfSites; k++) {
            // diff.loglikwndtpd returns log likelihood
            logP += siteLogLd[currentIndex][k];
        } // end k
        return logP;
    }

    // suppose only used by unit test
//    public void getBranchLikelihoods(double[] branchLdOut) {
//        System.arraycopy(branchLogLd[currentBrLdIndex], 0, branchLdOut, 0, branchLdOut.length);
//    }

    // ======= getters =======

//    public int getNrOfSites() {
//        return nrOfSites;
//    }

    public int getBranchNr() {
        return branchNr;
    }

    // = nrOfSites, for unit test
//    public int getBranchLdSize() {
//        // branchLd[][root] == null
//        return branchLogLd[0].length;
//    }
} // class