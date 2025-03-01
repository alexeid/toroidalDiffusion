package toroidaldiffusion.evolution.likelihood;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.State;
import toroidaldiffusion.ToroidalUtils;
import toroidaldiffusion.WrappedBivariateDiffusion;
import toroidaldiffusion.WrappedBivariateNormal;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import static java.lang.Math.sqrt;

public class PhyloWrappedBivariateDiffusion extends GenericDATreeLikelihood {

    // 2 values
    final public Input<Function> muInput = new Input<>("mu", "the mean of the stationary distribution.");
    // 2 values
    final public Input<Function> sigmaInput = new Input<>("sigma", "the two variance terms.");
    // 3 values
//    final public Input<Function> alphaInput = new Input<>("alpha", "the three drift terms.");
    final public Input<Function> driftInput = new Input<>("drift", "the two drift terms.");
    final public Input<Function> driftCorrInput = new Input<>("driftCorr", "the correlation of two drift terms, " +
            "ranged within (-1, 1), so that it always satisfies alpha1*alpha2 > alpha3^2.");

    //boolean used to switch between recalculate_branchLogP & recalculate_total_LogP
    final public Input<Boolean> sitewiseUpdateInput = new Input<>("useSitewiseUpdate",
            "If true, update likelihood site by site instead of full branch recalculation", true);

    protected boolean sitewiseUpdate;


    /****** calculation engine ******/
    protected WrappedBivariateDiffusion diff = new WrappedBivariateDiffusion();
//    protected double[] muArr; // stationary mean of the diffusion
//    protected double[] sigmaArr; // diffusion coefficient
//    protected double[] driftArr; // drift
//    protected double driftCorr; // ranged within (-1, 1), so that it always satisfies alpha1*alpha2 > alpha3^2.

//    protected BeagleTreeLikelihood beagle;
    /**
     * calculation engine for each branch, excl. root index, nrOfNodes-1
     */
    protected DABranchLikelihoodCore[] daBranchLdCores;
//    protected DABranchLikelihoodCore daRootLdCores;

    /**
     * multi-threading {@link DABranchLikelihoodCallable}
     */
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<>();


    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
//    protected int hasDirt; // TODO not need?

    /****** caching ******/

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] branchLengths; // length = nodeCount-1
    protected double[] storedBranchLengths;


    /**
     * caching log likelihoods for each of branch
     * root index is used to store frequencies prior at root
     */
    protected double[] branchLogLikelihoods; // length = nodeCount
    protected double[] storedBranchLogLikelihoods;

    /**
     * Total likelihoods
     */
    protected double totalLK;
    protected double storedTotalLK;


    public PhyloWrappedBivariateDiffusion(){}

    @Override
    public void initAndValidate() {
        // tree, and all nodes values
        super.initAndValidate();

        if (muInput.get().getDimension() != 2)
            throw new IllegalArgumentException("Expected two mu in 'muInput'. Found: " + muInput.get().getDimension());
        if (sigmaInput.get().getDimension() != 2)
            throw new IllegalArgumentException("Expected two sigma in 'sigmaInput'. Found: " + sigmaInput.get().getDimension());
        if (driftInput.get().getDimension() != 2)
            throw new IllegalArgumentException("Expected two drift terms in 'driftInput'. Found: " + driftInput.get().getDimension());

        sitewiseUpdate = sitewiseUpdateInput.get();
        // set mu, sigma, alpha
        setDiffusionParams();

        // no pattern, use getSiteCount()
        final int siteCount = daTreeModel.getSiteCount();

        // nodeCount = 2 * ntips - 1
        final int nodeCount = tree.getNodeCount();
        // caching branch lengths and log likelihoods for each of branch,
        branchLengths = new double[nodeCount-1];
        storedBranchLengths = new double[nodeCount-1];
        // root index is used to store frequencies prior at root
        branchLogLikelihoods = new double[nodeCount];
        storedBranchLogLikelihoods = new double[nodeCount];

        // branch Likelihood excl. root index
        daBranchLdCores = new DABranchLikelihoodCore[nodeCount-1];
        // init likelihood core using branch index (child node index below the branch)
        for (int n = 0; n < getRootIndex(); n++) {
            final Node node = tree.getNode(n);
            assert node.getNr() == n;
            // make every node dirty
            node.makeDirty(Tree.IS_FILTHY);
            // init by the node below the branch
            daBranchLdCores[n] = new DABranchLikelihoodCore(n, siteCount, diff);
        }
        tree.getRoot().makeDirty(Tree.IS_FILTHY); //init: recalculate the whole tree

        totalLK = noCachLogP();

        //bug:
//        storedTotalLK = totalLK;

        // root special
//        daRootLdCores = new DABranchLikelihoodCore(getRootIndex(), siteCount, diff);

        //TODO  multi-threading uses likelihoodCallers

//        for (int n = 0; n < getRootIndex(); n++) { // n is branchNr
//            List<DABranchLikelihoodCore> cores = new ArrayList<>();
//            cores.add(daBranchLdCores[n]);
//            likelihoodCallers.add(new DABranchLikelihoodCallable(cores));
//        }
//        // add root special calculation method
//        List<DABranchLikelihoodCore> cores = Collections.singletonList(daRootLdCores);
//        likelihoodCallers.add(new DABranchLikelihoodCallable(cores));
    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */

    @Override
    public double calculateLogP() {
        if (sitewiseUpdate) {
            logP = recalculateTotalLogP();
        } else {
            logP = noCachLogP();
        }
        return logP;
    }


//code style: changes to capital
    public double noCachLogP() {
        //TODO beagle
//        if (beagle != null) {
//            logP = beagle.calculateLogP();
//            return logP;
//        }

        // set diffusion params before computing likelihood
        setDiffusionParams();

        // rootIndex = getNodeCount - 1
        final int rootIndex = getRootIndex();
        // branch likelihoods indexes excludes root index
        for (int n = 0; n < rootIndex-1; n++) {
            final Node node = tree.getNode(n);
            DABranchLikelihoodCore daBranchLdCore = daBranchLdCores[n];
                try {
                        final Node parent = node.getParent();
    //                    double branchTime = this.branchLengths[nodeNr];
                        double branchTime = (parent.getHeight() - node.getHeight());

                    // Check for invalid branch times
                    if (branchTime == 0)
                        throw new UnsupportedOperationException("0 branch length, such as SA, not supported!");
                    if (branchTime < 1e-10)
                        throw new ArithmeticException("Reject proposal: time is 0 at the branch!");

                    try {
                        // must call updateBranch(node) before use computeBranchLK  //get branchTime first
                        this.branchLogLikelihoods[n] = daBranchLdCore.computeBranchLK(daTreeModel, parent, node, branchTime);
                    } catch (Exception e) {
//                        final Node parent = node.getParent();
                        System.err.println("parent (" + parent.getNr() +") = " + parent.getID() +
                                ", node (" + node.getNr() + ") = " + node.getID() + ", branchTime = " + branchTime);
                        System.err.println("mu = " + Arrays.toString(getMuArr()) +
                                ", sigma = " + Arrays.toString(sigmaInput.get().getDoubleValues()) +
                                ", alpha = " + Arrays.toString(getAlphaArr(driftInput.get().getDoubleValues(), getDriftCorr())));
                        throw new RuntimeException(e);
                    }

//            System.out.println("logP = " + logP);
            } catch (ArithmeticException e) {
                Log.err.println(e.getMessage());
                return Double.NEGATIVE_INFINITY;
            }
        } // end n loop

        // caching root likelihood branchLogLikelihoods[rootIndex]
            Node root = tree.getRoot();

            //remove IF tree.isDirty
                // the pairs of angles flatten to 1d, length = nsite * 2;
            double[] rootValues = daTreeModel.getNodeValue(root);
                // rootIndex = nNode - 1
            this.branchLogLikelihoods[rootIndex] = calculateRootLogLikelihood(rootValues, diff);


        // sum logP
        logP =0;
        // sum up all branch ld, plus frequencies at root
        for (int i = 0; i < branchLogLikelihoods.length; i++) {
            logP += branchLogLikelihoods[i];
        }

        //TODO Scaling

//        System.out.println("tree logP = " + logP);
        return logP;
    }

    public double recalculateTotalLogP() {
        setDiffusionParams();

        double logP = this.totalLK;
        final int rootIndex = getRootIndex();

        //three branches are affected by each single changed site
        int affectedCount = 3;
        double[] oldTerms = new double[affectedCount];
        double[] newTerms = new double[affectedCount];

        // Process all branches except root
        for (int n = 0; n < rootIndex-1; n++) {
            final Node node = tree.getNode(n);
            // Get changed sites
            int[] changedSites = determineChangedSites(node);

            if (updateBranch(node) != Tree.IS_CLEAN) {
                final int nodeNr = node.getNr();
                final Node parent = node.getParent();
                double branchTime = this.branchLengths[nodeNr];

                // Get the core for this branch
                DABranchLikelihoodCore daBranchLdCore = daBranchLdCores[n];

                // Calculate old and new terms for changed sites
                for (int i = 0; i < changedSites.length; i++) {
                    int siteIndex = changedSites[i];
                    oldTerms[i] = daBranchLdCore.computeSingleSiteLK(
                            daTreeModel.getNodeValue(parent),
                            daTreeModel.getNodeValue(node),
                            storedBranchLengths[nodeNr],
                            siteIndex
                    );
                    newTerms[i] = daBranchLdCore.computeSingleSiteLK(
                            daTreeModel.getNodeValue(parent),
                            daTreeModel.getNodeValue(node),
                            branchTime,
                            siteIndex
                    );
                }
            }

        }

        // Update total likelihood by site
        logP = update_log_LK(logP, oldTerms, newTerms);

//        // Handle root
//        Node root = tree.getRoot();
//        if (updateRoot(root) != Tree.IS_CLEAN) {
//            double[] rootValues = daTreeModel.getNodeValue(root);
//            this.branchLogLikelihoods[rootIndex] = calculateRootLogLikelihood(rootValues, diff);
//            logP += this.branchLogLikelihoods[rootIndex];
//        }

        this.totalLK = logP;
        return logP;
    }

    protected class UpdateStats {
        protected int updateCount = 0;
        protected double maxDelta = 0;
        protected int RECOMPUTE_THRESHOLD = 1000;
        protected static final int MIN_THRESHOLD = 100;
        protected static final int MAX_THRESHOLD = 10000;
        protected static final double ERROR_TOLERANCE = 1e-6;
    }

    protected UpdateStats updateStats = new UpdateStats();


    protected double update_log_LK(double current_total, double[] oldTerms, double[] newTerms){
        // Calculate sums
        double oldSum = 0;
        for (double term : oldTerms) oldSum += term;

        double newSum = 0;
        for (double term : newTerms) newSum += term;

        // Update total
        double updated = current_total - oldSum + newSum;

        // Track updates
        updateStats.updateCount++;
        double delta = Math.abs(newSum - oldSum);
        updateStats.maxDelta = Math.max(updateStats.maxDelta, delta);

        // Check if recomputation needed
        if (updateStats.updateCount >= updateStats.RECOMPUTE_THRESHOLD) {
            // Recompute from scratch
            double exact = noCachLogP();
            double error = Math.abs(updated - exact);

            // Adjust threshold based on error
            if (error > UpdateStats.ERROR_TOLERANCE) {
                // If error too large, reduce threshold to check more often
                updateStats.RECOMPUTE_THRESHOLD = Math.max(
                        updateStats.RECOMPUTE_THRESHOLD / 2,
                        UpdateStats.MIN_THRESHOLD
                );
                updated = exact;  // Use exact value if error too large
            } else {
                // If error acceptable, allow more updates before next check
                updateStats.RECOMPUTE_THRESHOLD = Math.min(
                        updateStats.RECOMPUTE_THRESHOLD * 2,
                        UpdateStats.MAX_THRESHOLD
                );
            }

            // Reset counters
            updateStats.updateCount = 0;
            updateStats.maxDelta = 0;
        }

        return updated;

    }

    /**
     * Returns an array of site indices that have changed for the given node.
     * If no sites changed, returns an empty array.
     */
    public int[] determineChangedSites(final Node node) {
        // Tips never change
        if (node.isLeaf()) return new int[0];

        // Calculate array index for this node
        int nTips = tree.getLeafNodeCount();
        final int nodeNr = node.getNr();
        int arrI = nodeNr - nTips;

        // Validate index
        if (arrI >= internalNodesParam.getMinorDimension2())
            throw new IllegalArgumentException("Internal node index (" + nodeNr + ") - nTips (" +
                    tree.getLeafNodeCount() + ") cannot > the internal node sequences parameter 2nd minor dim (" +
                    internalNodesParam.getMinorDimension2() + ") ! ");

        // Get number of sites
        int nSites = daTreeModel.getSiteCount();

        // Use a list to collect changed site indices
        List<Integer> changedSites = new ArrayList<>();

        // For each site, check if either phi or psi is dirty
        int start = arrI * internalNodesParam.getMinorDimension1();
        for (int k = 0; k < nSites; k++) {
            int phiIndex = start + (k * 2);
            int psiIndex = phiIndex + 1;

            if (internalNodesParam.isDirty(phiIndex) || internalNodesParam.isDirty(psiIndex)) {
                changedSites.add(k);
            }
        }


        // Convert to primitive int array
        int[] result = new int[changedSites.size()];
        for (int i = 0; i < changedSites.size(); i++) {
            result[i] = changedSites.get(i);
        }

        return result;
    }


//    protected int[] determineChangedSites(Node parent, Node node) { //find out changed site index
//        // For storing changed site indices
//        boolean[] changedFlags = new boolean[daTreeModel.getSiteCount()];
//        int changedCount = 0;
//
//        // Check both node and parent if they're internal nodes
//        Node[] nodesToCheck = {parent, node};
//        for (Node currentNode : nodesToCheck) {
//            if (!currentNode.isLeaf()) {
//                int arrI = currentNode.getNr() - tree.getLeafNodeCount();
//                int start = arrI * internalNodesParam.getMinorDimension1();
//
//                // Check each site
//                for (int k = 0; k < daTreeModel.getSiteCount(); k++) {
//                    if (!changedFlags[k]) {  // Only check if not already marked changed
//                        int phiIndex = start + (k * 2);
//                        int psiIndex = phiIndex + 1;
//
//                        if (internalNodesParam.isDirty(phiIndex) ||
//                                internalNodesParam.isDirty(psiIndex)) {
//                            changedFlags[k] = true;
//                            changedCount++;
//                        }
//                    }
//                }
//            }
//        }
//
//        // Create result array with exact size needed
//        int[] changeSites = new int[changedCount];
//        int pos = 0;
//        for (int k = 0; k < changedFlags.length; k++) {
//            if (changedFlags[k]) {
//                changeSites[pos++] = k;
//            }
//        }
//
//        return changeSites;
//    }




        /**
         * Main method to check if this branch LK requires update or not
         * @param node  the child node used to identify the branch
         * @return  the integer to represent if it is dirty using {@link Tree#IS_CLEAN} and so on.
         */
    protected int updateBranch(final Node node) {
        // the branch between node and parent
        // root is excluded from node when creating DABranchLikelihoodCore
        final Node parent = node.getParent();
        final int nodeNr = node.getNr();
        final int parentNum = parent.getNr();

        // if tips, always false
        //TODO the method has assumptions, check the detail before use
        boolean seqUpdate = isInternalNodeSeqDirty(node) || isInternalNodeSeqDirty(parent); //return i
//        boolean seqUpdate = false;

        int nodeUpdate = node.isDirty() | parent.isDirty();
//        int nodeUpdate = Tree.IS_DIRTY;

//        final double branchRate = 1.0; //TODO branchRateModel.getRateForBranch(node);
        // do not use getLength, code below to save time
        final double branchTime = (parent.getHeight() - node.getHeight()); //* branchRate;

        if (branchTime == 0)
            throw new UnsupportedOperationException("0 branch length, such as SA, not supported !");
        if (branchTime < 1e-10)
            throw new ArithmeticException("Reject proposal : " +
                    "time is 0 at the branch between parent node " + parentNum + " and node " + nodeNr +
                    " !\n" + "branch length = " + node.getLength() );//+ ", branchRate = " + branchRate);

        //TODO how to distinguish branch len change and internal node seq change, when topology is same

        // ====== 1. update the transition probability matrix(ices) if the branch len changes ======
        if (seqUpdate || nodeUpdate != Tree.IS_CLEAN || branchTime != branchLengths[nodeNr]) {
            // signal to recalculate
            nodeUpdate |= Tree.IS_DIRTY;

            // update branch length here
//            this.branchLengths[nodeNr] = branchTime;
        }

        return nodeUpdate;  //return tree_dirty: no caching likelihood
    }

    //update_total logLK of changed sites
//    protected int update_total_lk (double){
//
//    }





    // implement caching
    protected int updateRoot(final Node root) {
        //TODO the method has assumptions, check the detail before use
        boolean seqUpdate = isInternalNodeSeqDirty(root);
//        boolean seqUpdate = false;
        int nodeUpdate = root.isDirty();
//        int nodeUpdate = Tree.IS_DIRTY;

        if (seqUpdate || nodeUpdate != Tree.IS_CLEAN) {
            // signal to recalculate
            nodeUpdate |= Tree.IS_DIRTY;
        }
        return nodeUpdate;
    }

    // the stationary density on the root
    protected double calculateRootLogLikelihood(double[] rootValues, WrappedBivariateDiffusion diff) {
        if (rootValues.length != daTreeModel.getSiteCount() * 2)
            throw new IllegalArgumentException("Incorrect length on root sequences ! " + rootValues.length);
        double logP = 0;
        // s is site index
        for (int s = 0; s < daTreeModel.getSiteCount(); s++) {
            // two angles, phi and psi, at the parent node and the sth site.
            double phi0 = rootValues[s * 2];
            double psi0 = rootValues[s * 2 + 1];
            // TODO seem not working ?
            // stationary density
            logP += diff.loglikwndstat(phi0, psi0);
        }
        return logP;
    }

    @Override
    public void init(PrintStream out) {
        super.init(out);
        out.print("alpha3\t");
        out.print("WNsd1\tWNsd2_3\tWNsd4\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        super.log(sample, out);
        double[] driftArr = driftInput.get().getDoubleValues(); // drift
        double driftCorr = getDriftCorr();
        double[] alpha = getAlphaArr(driftArr, driftCorr);
        // alpha3
        out.print(alpha[2] + "\t");
        double[] sigma = sigmaInput.get().getDoubleValues();
        // WN sd, [sd1, sd2_3; sd2_3, sd4], so arr = [sd1, sd2_3, sd4]
        double[] wnsd = WrappedBivariateNormal.getWNSd(alpha, sigma);
        out.print(wnsd[0] + "\t" + wnsd[1] + "\t" + wnsd[2] + "\t");
    }

    // refresh muarr, alphaarr, sigmaarr for computing likelihood
    protected void setDiffusionParams() {
        double[] muArr = getMuArr(); // stationary mean of the diffusion
        double[] sigmaArr = sigmaInput.get().getDoubleValues(); // diffusion coefficient
        double[] driftArr = driftInput.get().getDoubleValues(); // drift
        double driftCorr = getDriftCorr();

        double[] alphaarr = getAlphaArr(driftArr, driftCorr);
        // init WrappedBivariateDiffusion here, setParameters(muarr, alphaarr, sigmaarr) once.
        // use diff.loglikwndtpd(phi0, psi0, phit, psit) later when compute likelihood
        diff.setParameters(muArr, alphaarr, sigmaArr);
    }

    protected double[] getMuArr() {
        double[] muArr = muInput.get().getDoubleValues();
        if (muArr[0] < 0 || muArr[0] > ToroidalUtils.MAX_ANGLE_VALUE ||
                muArr[1] < 0 || muArr[1] > ToroidalUtils.MAX_ANGLE_VALUE)
            throw new IllegalArgumentException("Mu should be [0, 2PI], but muArr[0] = " + muArr[0] +
                    ", muArr[1] = " + muArr[1] + " !");
        return muArr;
    }
    protected double getDriftCorr() {
        double driftCorr = driftCorrInput.get().getArrayValue();
        if (driftCorr <= -1.0 || driftCorr >= 1.0)
            throw new IllegalArgumentException("Drifts correlation must be within (-1, 1). Found: " + driftCorr);
        return driftCorr;
    }

    // compute A given two drifts and their correlation
    protected double[] getAlphaArr(double[] driftArr, double driftCorr) {
        // corr within (-1, 1), so that always alpha1*alpha2 > alpha3^2
        double alpha3 = sqrt(driftArr[0]*driftArr[1]) * driftCorr;
        double[] alphaarr = new double[]{driftArr[0], driftArr[1], alpha3};

        if (alphaarr[0] * alphaarr[1] <= alphaarr[2] * alphaarr[2])
            throw new IllegalArgumentException("alpha1 * alpha2 must > alpha3 * alpha3 ! But alpha = {" +
                    alphaarr[0] + ", " + alphaarr[1] + ", " + alphaarr[2] +  "} is invalid.");
        return alphaarr;
    }

    //****** multi-threading ******//


    // ArithmeticException if branch time < 1e-10
    class DABranchLikelihoodCallable implements Callable<Double> {

        private final List<DABranchLikelihoodCore> brLDCores;

        // by group
        public DABranchLikelihoodCallable(List<DABranchLikelihoodCore> brLDCores) {
            throw new UnsupportedOperationException("in dev !");
//            this.brLDCores = brLDCores;
        }

        @Override
        public Double call() throws Exception {
//            try {
            double logP = 0;
            for (DABranchLikelihoodCore core : brLDCores) {
                final int branchNr = core.getBranchNr();
                final Node node = tree.getNode(branchNr);

                // caching branchLogLikelihoods[nodeNr]
                if (updateBranch(node) != Tree.IS_CLEAN)
                    branchLogLikelihoods[branchNr] = core.sumSiteLogLikelihood();

                logP += branchLogLikelihoods[branchNr];
            }
//            } catch (Exception e) {
//                System.err.println("Something wrong to calculate branch likelihood above node " +
//                        branchNr + " during multithreading !");
//                e.printStackTrace();
//                System.exit(0);
//            }
//            System.out.println("Branch likelihood logP = " + branchLogLikelihoods[branchNr] + " above node " + branchNr);
            return logP;
        }

    }


    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
//        if (beagle != null) {
//            return beagle.requiresRecalculation();
//        }

        // must be called before super
        setDiffusionParams();
        return super.requiresRecalculation();

// TODO check isDirtyCalculation()?
//        if (daTreeModel.getTipsValuesParam().somethingIsDirty()) {
//            return true;
//        }
//        if (daTreeModel.getInternalNodesValuesParam().somethingIsDirty()) {
//            return true;
//        }
//        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
//            //m_nHasDirt = Tree.IS_DIRTY;
//            return true;
//        }
        // TODO check
//        for (Node node : tree.getNodesAsArray())
//            if (node.isDirty() != Tree.IS_CLEAN)
//                return true;
//        return tree.somethingIsDirty();
    }

    @Override
    public void store() {
//        if (beagle != null) {
//            beagle.store();
//            super.store();
//            return;
//        }

        for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
            daBrLdCore.store();
        //TODO estimating internal node sequences needs store/restore?
        // param store is protected?

        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, getNrOfBranches());
        System.arraycopy(branchLogLikelihoods, 0, storedBranchLogLikelihoods, 0, getNrOfBranches()+1);

        storedTotalLK = totalLK;
        // store logP after all stored
        super.store(); // storedLogP = logP; isDirty = false
    }

    @Override
    public void restore() {
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }

        // must refresh diffusion param
        setDiffusionParams();

        for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
            daBrLdCore.restore();
        //TODO estimating internal node sequences needs store/restore?
//        daTreeModelInput.get().getInternalNodesValuesParam().restore();

        // pass reference in restore, but have to copy array in store.
//        double[] tmp1 = branchLengths;
//        branchLengths = storedBranchLengths;
//        storedBranchLengths = tmp1;
//
//        double[] tmp2 = branchLogLikelihoods;
//        branchLogLikelihoods = storedBranchLogLikelihoods;
//        storedBranchLogLikelihoods = tmp2;

        System.arraycopy(storedBranchLengths, 0, branchLengths, 0, getNrOfBranches());
        System.arraycopy(storedBranchLogLikelihoods, 0, branchLogLikelihoods, 0, getNrOfBranches()+1);

        totalLK = storedTotalLK;

        // store logP after all stored
        super.restore(); // logP = storedLogP; isDirty = false
    }

    // for testing, nodeIndex is the child node below this branch
//    public void getBranchLdFromCore(int nodeIndex, double[] branchLdOut) {
//        daBranchLdCores[nodeIndex].getBranchLikelihoods(branchLdOut);
//    }

    public int getRootIndex() {
        int rootIndex = tree.getNodeCount() - 1;
        assert rootIndex == tree.getRoot().getNr();
        return rootIndex;
    }

    public int getNrOfBranches() {
        int branchCount = tree.getNodeCount() - 1;
        assert branchCount == branchLengths.length;
        return branchCount;
    }

    /**
     * @param branchNr the node Nr of child node below the branch
     * @return {@link DABranchLikelihoodCore} containing cached values.
     */
    public DABranchLikelihoodCore getDaBranchLdCores(int branchNr) {
        return daBranchLdCores[branchNr];
    }


} // class DATreeLikelihood

