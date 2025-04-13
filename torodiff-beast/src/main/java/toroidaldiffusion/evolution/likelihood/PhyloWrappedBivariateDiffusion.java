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
        tree.getRoot().makeDirty(Tree.IS_FILTHY);
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
                // caching branchLogLikelihoods[nodeNr]
                if (updateBranch(node) != Tree.IS_CLEAN) {
//                    this.branchLogLikelihoods[n] = daBranchLdCore.calculateBranchLogLikelihood();
                    final int nodeNr = node.getNr();
                    final Node parent = node.getParent();
                    double branchTime = this.branchLengths[nodeNr];

                    try {
                        // Check what type of update is needed
                        int nodeUpdate = node.isDirty() | parent.isDirty();
                        boolean branchLengthUpdate = branchTime != storedBranchLengths[nodeNr];
                        boolean seqUpdate = isInternalNodeSeqDirty(node) || isInternalNodeSeqDirty(parent);

                        if (nodeUpdate != Tree.IS_CLEAN || branchLengthUpdate) {
                            // must call updateBranch(node) before use computeBranchLK
                            this.branchLogLikelihoods[n] = daBranchLdCore.computeBranchLK(daTreeModel, parent, node, branchTime);
                        } else if (seqUpdate) {
                            //determine changed site in the sequence: for both parent and children
//                            int changeSiteindexParent = determineChangedsiteIndex(parent);
//                            int changeSiteindexChildren = determineChangedsiteIndex(node);
//
                            //If changeSites List = empty, boolean NoChangedSites = True
                            List changeSiteindexParent = determineChangedsiteIndex(parent);
                            List changeSiteindexChildren = determineChangedsiteIndex(node);

                            //The whole branch is changed (Gibb's sampler) - recalculate whole branch likelihood
                            //todo: check conditions
                            if (changeSiteindexParent.size() == daTreeModel.getSiteCount() || changeSiteindexChildren.size() == daTreeModel.getSiteCount()) {

                                // Full recalculation for all sites
                                this.branchLogLikelihoods[n] = daBranchLdCore.computeBranchLK(daTreeModel, parent, node, branchTime);
                            }

                            //Some sites have been changed e.g. Wrapped RandomWalk Operator
                            else {
                                boolean noChangedSitesInChildren = detectChangedSites(changeSiteindexChildren);

                                if (!noChangedSitesInChildren) { //changed sites found in Children, use children index
                                    this.branchLogLikelihoods[n] = daBranchLdCore.computeBranchLKbySite(daTreeModel, parent, node, branchTime, changeSiteindexChildren);
                                }
                                else{ //changed sites found in Parent, use parent index
                                    if (changeSiteindexParent.isEmpty()) {
                                        throw new IllegalArgumentException("No changed sites found in both parents and children");
                                    }
                                    this.branchLogLikelihoods[n] = daBranchLdCore.computeBranchLKbySite(daTreeModel, parent, node, branchTime, changeSiteindexParent);
                                }
                            }

//                            if (changeSiteindexChildren != -1) {
//                                this.branchLogLikelihoods[n] = daBranchLdCore.computeBranchLKbySite(daTreeModel, parent, node, branchTime, changeSiteindexChildren);
//                            } else {
//                                this.branchLogLikelihoods[n] = daBranchLdCore.computeBranchLKbySite(daTreeModel, parent, node, branchTime, changeSiteindexParent);
//                            }
                        }

                    } catch (Exception e) {
//                        final Node parent = node.getParent();
                        System.err.println("parent (" + parent.getNr() +") = " + parent.getID() +
                                ", node (" + node.getNr() + ") = " + node.getID() + ", branchTime = " + branchTime);
                        System.err.println("mu = " + Arrays.toString(getMuArr()) +
                                ", sigma = " + Arrays.toString(sigmaInput.get().getDoubleValues()) +
                                ", alpha = " + Arrays.toString(getAlphaArr(driftInput.get().getDoubleValues(), getDriftCorr())));
                        throw new RuntimeException(e);
                    }
                }
//            System.out.println("logP = " + logP);
            } catch (ArithmeticException e) {
                Log.err.println(e.getMessage());
                return Double.NEGATIVE_INFINITY;
            }
        } // end n loop

        // caching root likelihood branchLogLikelihoods[rootIndex]
        Node root = tree.getRoot();
        if (updateRoot(root) != Tree.IS_CLEAN) {
            // the pairs of angles flatten to 1d, length = nsite * 2;
            double[] rootValues = daTreeModel.getNodeValue(root);
            // rootIndex = nNode - 1
            this.branchLogLikelihoods[rootIndex] = calculateRootLogLikelihood(rootValues, diff);
        }

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

    /**
     *
     * @param node
     * @return int k
     */

    public List<Integer> determineChangedsiteIndex(final Node node) {

        // Create a list to store all changed site indices
        List<Integer> changedSites = new ArrayList<>();

        //tips never change
        if(node.isLeaf()) return changedSites;

        //get the numb of tips
        int nTips = tree.getLeafNodeCount();
        //get the node number of the current node
        final int nodeNr = node.getNr();
        //the index of the current node in internalNodesParam
        int arrI = nodeNr - nTips;

        // minor dim 2 = number of internal nodes
        if (arrI >= internalNodesParam.getMinorDimension2())
            throw new IllegalArgumentException("Internal node index (" + nodeNr + ") - nTips (" +
                    tree.getLeafNodeCount() + ") cannot > the internal node sequences parameter 2nd minor dim (" +
                    internalNodesParam.getMinorDimension2() + ") ! ");
        //get number of sites
        int nSites = daTreeModel.getSiteCount();

        // minor dim 1 = nsite * 2
        int start = arrI * internalNodesParam.getMinorDimension1();

        for (int k = 0; k < nSites; k++) {
            int phiIndex = start + (k * 2);
            int psiIndex = phiIndex + 1;

            if (internalNodesParam.isDirty(phiIndex) || internalNodesParam.isDirty(psiIndex)) {
                changedSites.add(k);
            }
        }

        return changedSites;
    }


    /**
     * Detect if any sites have been changed
     */
    public boolean detectChangedSites(List<Integer> changedSites) {
        boolean noChangeSitesdFound = changedSites.isEmpty();
     return noChangeSitesdFound;
    }

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
        boolean seqUpdate = isInternalNodeSeqDirty(node) || isInternalNodeSeqDirty(parent);
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
            this.branchLengths[nodeNr] = branchTime;
        }

        return nodeUpdate;
    }

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

    /**
     * Get diff
     * @return
     */
    public WrappedBivariateDiffusion getDiff() {
        return diff;
    }


} // class DATreeLikelihood

