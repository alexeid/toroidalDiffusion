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
    // 1st dimension is matrix index (current, stored),
    // 2nd is nrOfSites, Ld across categories are integrated
    protected double[][] branchLogLd;

    protected int currentBrLdIndex = 0;
    protected int storedBrLdIndex = 0;

    /**
     * caching likelihood per site at one node.
     * size = siteCount
     */
    protected double[] siteLd;

    protected boolean useScaling = false;

    protected double[][] scalingFactors; //TODO

    static final public double SCALING_THRESHOLD = 1.0E-150; // MAX_VALUE 1.7*10^308
    double SCALE = 2;

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
        branchLogLd = new double[2][nrOfSites];
        siteLd = new double[nrOfSites];
        //TODO diff.setParameters(muarr, alphaarr, sigmaarr); ?
    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        storedBrLdIndex = currentBrLdIndex;
    }

    /**
     * Store current state
     */
    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        int tmp = currentBrLdIndex;
        currentBrLdIndex = storedBrLdIndex;
        storedBrLdIndex = tmp;
    }

    @Override
    public void unstore() {
        currentBrLdIndex = storedBrLdIndex;
    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        branchLogLd = null;
        currentBrLdIndex = -1;
        storedBrLdIndex = -1;

        scalingFactors = null;
    }

    @Override
    public void setUseScaling(double scale) {
        useScaling = (scale != 1.0);

        if (useScaling) {
            scalingFactors = new double[2][nrOfSites];
        }
    }

    //============ branch likelihood ============

    /**
     * use before {@link #calculateBranchLd(double[], double[])}
     */
    public void setBranchLdForUpdate() {
        currentBrLdIndex = 1 - currentBrLdIndex; // 0 or 1
    }


    //TODO need to cache per site to avoid recalculation, when only the sequence at a site is changed
    public void calculateBranchLd(final double[] parentNodeValues, final double[] childNodeValues) {

        for (int k = 0; k < nrOfSites; k++) {
            // pairs of values, dimension is 2 (angles) * N_sites
            double phi0 = parentNodeValues[k];
            double psi0 = parentNodeValues[k + 1];

            double phit = childNodeValues[k];
            double psit = childNodeValues[k + 1];

            // diff.setParameters(muarr, alphaarr, sigmaarr), only once in the init method
            branchLogLd[currentBrLdIndex][k] = diff.loglikwndtpd(phi0, psi0, phit, psit);

            if (branchLogLd[currentBrLdIndex][k] == 0) {
                throw new RuntimeException("\nBranch above node " + getBranchNr() + " likelihood = 0 !\n" +
                        "At site " + k); //+ ", child node = " + childNode + ", parent node = " + parentNode);
            }

        } // end k  nrOfSites


        if (useScaling) {
            throw new UnsupportedOperationException("in dev");
//            scaleBranchLds(nodeParent);
        }

    }

    /**
     * Calculates log likelihood at this branch.
     * Multiple given site likelihoods, and then log the product.
     * The likelihood input has been integrated across categories.
     * It can be also used for the site likelihoods at the root.
     * @return           logged likelihood
     */
    public double calculateBranchLogLikelihood() {
        final double[] siteLd = branchLogLd[currentBrLdIndex];

        double logP = 0;

        // siteLd.length == nrOfSites
        for (int k = 0; k < siteLd.length; k++) {
            // diff.loglikwndtpd returns log likelihood
            logP += siteLd[k];
        } // end k

        return logP;
    }

    // log likelihood at root given codon frequencies
    public double calculateRootLogLikelihood(double[] rootValues, WrappedBivariateDiffusion diff) {
        // TODO
        throw new UnsupportedOperationException("in dev");
    }

    /**
     * This function returns the scaling factor for that pattern by summing over
     * the log scalings used at each node. If scaling is off then this just returns
     * a 0.
     *
     * @return the log scaling factor
     */
    @Override
    public double getLogScalingFactor(int patternIndex_) {
//    	if (m_bUseScaling) {
//    		return -(m_nNodeCount/2) * Math.log(SCALE);
//    	} else {
//    		return 0;
//    	}
        double logScalingFactor = 0.0;
        if (useScaling) {
            throw new UnsupportedOperationException("in development");
//            for (int i = 0; i < nrOfNodes; i++) {
//                logScalingFactor += scalingFactors[currentBrLdIndex[i]][i][patternIndex_];
//            }
        }
        return logScalingFactor;
    }


    // suppose only used by unit test
    public void getBranchLikelihoods(double[] branchLdOut) {
        System.arraycopy(branchLogLd[currentBrLdIndex], 0, branchLdOut, 0, branchLdOut.length);
    }

    // ======= getters =======

    public int getNrOfSites() {
        return nrOfSites;
    }

    public int getBranchNr() {
        return branchNr;
    }

    // = nrOfSites, for unit test
    public int getBranchLdSize() {
        // branchLd[][root] == null
        return branchLogLd[0].length;
    }
} // class