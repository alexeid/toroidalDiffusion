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
    // dimension is nrOfSites, Ld across categories are integrated
    protected double[] branchLogLd;
    protected double[] storedBranchLogLd;

//    protected int currentBrLdIndex = 0;
//    protected int storedBrLdIndex = 0;

    /**
     * caching likelihood per site at one node.
     * size = siteCount
     */
//    protected double[] siteLd;

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
        branchLogLd = new double[nrOfSites];
        storedBranchLogLd = new double[nrOfSites];
//        siteLd = new double[nrOfSites];
        //TODO diff.setParameters(muarr, alphaarr, sigmaarr); ?
    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
//        storedBrLdIndex = currentBrLdIndex;
        System.arraycopy(branchLogLd, 0, storedBranchLogLd, 0, nrOfSites);
    }

    /**
     * Store current state
     */
    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
//        int tmp = currentBrLdIndex;
//        currentBrLdIndex = storedBrLdIndex;
//        storedBrLdIndex = tmp;
        System.arraycopy(storedBranchLogLd, 0, branchLogLd, 0, nrOfSites);
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
        branchLogLd = null;
        storedBranchLogLd = null;
//        currentBrLdIndex = -1;
//        storedBrLdIndex = -1;
    }


    //============ branch likelihood ============

    /**
     * use before {@link #calculateBranchLd(double[], double[], double)}
     */
//    public void setBranchLdForUpdate() {
//        currentBrLdIndex = 1 - currentBrLdIndex; // 0 or 1
//    }


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

            // diff.setParameters(muarr, alphaarr, sigmaarr), only once in the init method
            branchLogLd[k] = diff.loglikwndtpd(phi0, psi0, phit, psit);

            if (branchLogLd[k] == 0) {
                throw new RuntimeException("\nBranch above node " + getBranchNr() + " likelihood = 0 !\n" +
                        "At site " + k); //+ ", child node = " + childNode + ", parent node = " + parentNode);
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
        // siteLd.length == nrOfSites
        for (int k = 0; k < branchLogLd.length; k++) {
            // diff.loglikwndtpd returns log likelihood
            logP += branchLogLd[k];
        } // end k
        return logP;
    }

    // suppose only used by unit test
//    public void getBranchLikelihoods(double[] branchLdOut) {
//        System.arraycopy(branchLogLd[currentBrLdIndex], 0, branchLdOut, 0, branchLdOut.length);
//    }

    // ======= getters =======

    public int getNrOfSites() {
        return nrOfSites;
    }

    public int getBranchNr() {
        return branchNr;
    }

    // = nrOfSites, for unit test
//    public int getBranchLdSize() {
//        // branchLd[][root] == null
//        return branchLogLd[0].length;
//    }
} // class