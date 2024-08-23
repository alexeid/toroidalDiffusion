package toroidaldiffusion.evolution.likelihood;

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
    protected double[][] branchLd;

    protected int currentBrLdIndex = 0;
    protected int storedBrLdIndex = 0;

    /**
     * caching probability tables obtained from substitutionModel,
     * size = stateCount * stateCount
     */
    protected double[] probabilities;
    /**
     * caching likelihood per site at one node.
     * size = siteCount
     */
    protected double[] siteLd;

    protected boolean useScaling = false;

    protected double[][] scalingFactors; //TODO

    static final public double SCALING_THRESHOLD = 1.0E-150; // MAX_VALUE 1.7*10^308
    double SCALE = 2;

    /**
     * no initialization, for calculating site likelihoods at the root.
     * @param branchNr       the node Nr of child node below the branch
     * @param nrOfSites      number of sites (codon)
     */
    public DABranchLikelihoodCore(int branchNr, int nrOfSites) {
        this.branchNr = branchNr;
        this.nrOfSites = nrOfSites; // TODO impl data range to multithreading by sites

        siteLd = new double[nrOfSites]; // need init siteLd[] here

        initialize();
    }


    /**
     * initializes states, likelihood arrays.
     */
    @Override
    protected void initialize() {
        // the branch above the node
        // merged likelihood for all categories
        branchLd = new double[2][nrOfSites];


        // used to cache transition probability matrix P
        probabilities = new double[matrixSize];
        Arrays.fill(probabilities, 1.0);

        siteLd = new double[nrOfSites];
    }

    /**
     * @return the reference of array caching probability,
     *         which has to be here to make thread safe.
     */
    public double[] getProbRef() {
        return this.probabilities;
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
        branchLd = null;
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
     * use before {@link #calculateBranchLd(int[], int[], double[])}
     */
    public void setBranchLdForUpdate() {
        currentBrLdIndex = 1 - currentBrLdIndex; // 0 or 1
    }


    /**
     * Branch likelihood calculation with site model categories from state i to j.
     * The flattened transition probability matrix P(t) index : n = w + i * state + j.
     * @param i   state i in parent node.
     * @param j    state j in child node.
     * @param proportions       the proportions of sites in each category. length = nrOfCategories.
     * @return   the branch likelihood at one site
     */
    public double calculateBranchLdAtSite(final int i, final int j, final double[] proportions) {
        final double[] matrices = getCurrentMatrix();

        //branchLd is defined as the branch above the child node 1
        int w = i * nrOfStates + j;
        double branchLdSite = matrices[w] * proportions[0]; // i * nrOfStates + j

        for (int l = 1; l < nrOfCategories; l++) {
            w += matrixSize; // l * matrixSize + i * nrOfStates + j;
            branchLdSite += matrices[w] * proportions[l];
        } // end l nrOfCategories

        return branchLdSite;
    }

    /**
     * It calculates branch likelihoods at a branch,
     * and integrates branch likelihood across categories.
     *  @param statesParentNode  array of states at 'parent' node
     * @param statesChildNode    array of states at one 'child' node
     * @param proportions the proportions of sites in each category. length = nrOfCategories.
     */
    //TODO need to cache per site to avoid recalculation, when only the sequence at a site is changed
    public void calculateBranchLd(final int[] statesParentNode, final int[] statesChildNode,
                                  double[] proportions) {

        /**
         * Calculate DA branch likelihood per site per node.
         * matrix P(t) is flattened to n = w + i * state + j,
         * where n is index of flattened transition probability matrix (double[] matrices?),
         * i is parent state, j is child state, w is the category index.
         */

//        assert statesParentNode.length == nrOfSites && statesChildNode.length == nrOfSites;

        int i; // parent state
        int j; // child state
        // integrate over categories
        for (int k = 0; k < nrOfSites; k++) {

            i = statesParentNode[k]; // parent
            j = statesChildNode[k];

            if (j < nrOfStates) { // && state3 < nrOfStates && state2 < nrOfStates

                // branch likelihoods branchLd[currentBrLdIndex][], len is nrOfSites
                // this has been integrated over categories
                branchLd[currentBrLdIndex][k] = calculateBranchLdAtSite(i, j, proportions);

            } else {
                throw new UnsupportedOperationException("State out of range at site " + k +
                        ", child node state = " + j + " parent node state = " + i);

                //TODO child 1 has a gap or unknown state so treat it as unknown

//                    for (int j = 0; j < nrOfStates; j++) {
//                        branchLd3[v] = 1.0;
//                        v++;
//                    }
            }

            if (branchLd[currentBrLdIndex][k] == 0) {
                // transition probability matrix
                double[] matrices = getCurrentMatrix();
                for (int l = 0; l < nrOfCategories; l++)
                    System.err.println("Category " + l +  ": transition probability = " +
                            matrices[l * matrixSize + j * nrOfStates + i]);
//                System.err.println
                throw new RuntimeException("\nBranch above node " + getBranchNr() + " likelihood = 0 !\n" +
                        "At site " + k + ", child node state = " + j + ", parent node state = " + i);
            }

        } // end k  nrOfSites


//        } else {
//            throw new IllegalArgumentException("Every node must has states and length of " + nrOfSites + " ! \n" +
//                    "child 1 " + childNodeStates.length + ", parent " + node3States.length);
//        }


//        final double[] tmp = branchLd[currentBrLdIndex[nodeChild1]][nodeChild1];
//        for (int i=0; i < tmp.length; i++) {
//            if (tmp[i] == 0) {
//                // branchLd.length = nrOfSites * nrOfCategories, states.length = nrOfSites
//                int site = i % nrOfSites;
//                System.err.println("i = " + i + " site = " + site + " : " +
//                        "child 1 id = " + nodeChild1 + " state = " + stateIndex1[site] +
//                        //", child 2 id = " + nodeIndex2 + " state = " + stateIndex2[site] +
//                        ", parent id = " + nodeParent + " state = " + stateIndex3[site] +
//                        ", branchLd[" + i + "] = " + tmp[i]);
//            }
//        }

        if (useScaling) {
            throw new UnsupportedOperationException("in dev");
//            scaleBranchLds(nodeParent);
        }

//
//        int k =0;
//        for (int i = 0; i < patternCount; i++) {
//            double f = 0.0;
//
//            for (int j = 0; j < stateCount; j++) {
//                f += branchLd[currentPartialsIndices[nodeParent]][nodeParent][k];
//                k++;
//            }
//            if (f == 0.0) {
//                Logger.getLogger("error").severe("A branch likelihood (node index = " + nodeParent + ", pattern = "+ i +") is zero for all states.");
//            }
//        }
    }

    /**
     * Log likelihood calculation engine.
     * Multiple given site likelihoods, and then log the product.
     * The likelihood input has been integrated across categories.
     * It can be also used for the site likelihoods at the root.
     * @param siteLd     likelihoods (not logged) by sites (codons).
    //     * @param scalingThreshold   threshold for the product to call {@link Math#log(double)}.
     * @return           logged likelihood
     */
    public double logIntegratedLikelihood(double[] siteLd) {

        double product = 1.0;
        double logP = 0;

        // siteLd.length == nrOfSites
        for (int k = 0; k < siteLd.length; k++) {
            // multiple (not logged) site likelihoods
            product *= siteLd[k];

            // hard code to log when product is too small, Double.MAX_VALUE 1.79...e+308
            if (product < SCALING_THRESHOLD ) { // || siteLd[k] < SCALING_THRESHOLD
                // important check before implement log scaling
                if (product == 0)
                    throw new RuntimeException("Likelihood product -Inf ! " +
                            "\nlikelihood = " + siteLd[k] + " at site " + k);

                logP += Math.log(product); //+ getLogScalingFactor(k); TODO
                product = 1.0;
            } // if

        } // end k

        // log the rest
        if (product < 1.0)
            logP += Math.log(product); //+ getLogScalingFactor(k); TODO

        return logP;
    }


    /**
     * Calculates log likelihood at this branch, excluding root frequency prior.
     * The input branch likelihoods here have been integrated across categories.
     * frequency[] need to be added later in DATreeLikelihood
     */
    public double calculateBranchLogLikelihood() {
        return logIntegratedLikelihood(branchLd[currentBrLdIndex]);
    }

    // log likelihood at root given codon frequencies
    public double calculateRootLogLikelihood(NodeStates rootStates, double[] frequencies) {
//        int siteCount = rootStates.getSiteCount();
//        double[] siteLdAtRoot = new double[siteCount];

        int state;
        for (int k = 0; k < nrOfSites; k++) {
            // hard code for root node
            state = rootStates.getState(k); // 0-63
            siteLd[k] = frequencies[state];
        }

        return logIntegratedLikelihood(siteLd); //+ getLogScalingFactor(k); TODO
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
        System.arraycopy(branchLd[currentBrLdIndex], 0, branchLdOut, 0, branchLdOut.length);
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
        return branchLd[0].length;
    }
} // class