package toroidaldiffusion.evolution.likelihood;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.State;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

public class PhyloWrappedBivariateDiffusion extends GenericDATreeLikelihood {

//    final public Input<TreeLikelihood.Scaling> scaling = new Input<>("scaling",
//            "type of scaling to use, one of " + Arrays.toString(TreeLikelihood.Scaling.values()) +
//                    ". If not specified, the -beagle_scaling flag is used.",
//            TreeLikelihood.Scaling._default, TreeLikelihood.Scaling.values());

    // 2 values
    final public Input<Function> muInput = new Input<>("mu", "the mean of the stationary distribution.");
    // 2 values
    final public Input<Function> sigmaInput = new Input<>("sigma", "the two variance terms.");
    // 3 values
    final public Input<Function> alphaInput = new Input<>("sigma", "the three drift terms.");


    /****** calculation engine ******/
//    protected BeagleTreeLikelihood beagle;
    /**
     * calculation engine for each branch, excl. root index, nrOfNodes-1
     */
    protected DABranchLikelihoodCore[] daBranchLdCores;
    protected DABranchLikelihoodCore daRootLdCores;

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

    /****** TODO rest ******/
    // dealing with proportion of site being invariant
    double proportionInvariant = 0;


    public PhyloWrappedBivariateDiffusion(){}

    @Override
    public void initAndValidate() {
        // tree, and all nodes values
        super.initAndValidate();

        // no pattern, use getSiteCount()
        final int siteCount = daTreeModel.getSiteCount();

        // init DALikelihoodCore
        final int leafNodeCount = tree.getLeafNodeCount();
        // caching branch lengths and log likelihoods for each of branch,
        branchLengths = new double[leafNodeCount-1];
        storedBranchLengths = new double[leafNodeCount-1];
        // root index is used to store frequencies prior at root
        branchLogLikelihoods = new double[leafNodeCount];
        storedBranchLogLikelihoods = new double[leafNodeCount];

        // branch Likelihood excl. root index
        daBranchLdCores = new DABranchLikelihoodCore[leafNodeCount-1];
        // init likelihood core using branch index (child node index below the branch)
        for (int n = 0; n < getRootIndex(); n++) {
            final Node node = tree.getNode(n);
            assert node.getNr() == n;
            // make every node dirty
            node.makeDirty(Tree.IS_FILTHY);
            // init by the node below the branch
            daBranchLdCores[n] = new DABranchLikelihoodCore(n, siteCount);
        }
        tree.getRoot().makeDirty(Tree.IS_FILTHY);
        // root special
        daRootLdCores = new DABranchLikelihoodCore(getRootIndex(), siteCount);

        for (int n = 0; n < getRootIndex(); n++) { // n is branchNr
            List<DABranchLikelihoodCore> cores = new ArrayList<>();
            cores.add(daBranchLdCores[n]);
            likelihoodCallers.add(new DABranchLikelihoodCallable(cores));
        }
        // add root special calculation method
        List<DABranchLikelihoodCore> cores = Collections.singletonList(daRootLdCores);
        likelihoodCallers.add(new DABranchLikelihoodCallable(cores));
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
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {
        //TODO beagle
//        if (beagle != null) {
//            logP = beagle.calculateLogP();
//            return logP;
//        }

        // exclude root node, branches = nodes - 1
        final int rootIndex = getRootIndex();

        final double[] frequencies = substitutionModel.getFrequencies();
        NodeStates rootStates = nodesStates.getNodeStates(rootIndex);

        // branch likelihoods indexes excludes root index
        for (int n = 0; n < rootIndex; n++) {
            final Node node = tree.getNode(n);
            DABranchLikelihoodCore daBranchLdCore = daBranchLdCores[n];
            try {
                // caching branchLogLikelihoods[nodeNr]
                if (updateBranch(daBranchLdCore, node) != Tree.IS_CLEAN) {
                    this.branchLogLikelihoods[n] = daBranchLdCore.calculateBranchLogLikelihood();
                }
//            System.out.println("logP = " + logP);
            } catch (ArithmeticException e) {
                Log.err.println(e.getMessage());
                return Double.NEGATIVE_INFINITY;
            }
        } // end n loop
        // root special
        this.branchLogLikelihoods[rootIndex] =
                daRootLdCores.calculateRootLogLikelihood(rootStates, frequencies);

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


    protected int updateBranch(final DABranchLikelihoodCore daBranchLdCore, final Node node) {
        // the branch between node and parent
        // root is excluded from node when creating DABranchLikelihoodCore
        final Node parent = node.getParent();
        final int nodeNr = node.getNr();
        final int parentNum = parent.getNr();

        // if tips, always false
        //TODO cache per site to avoid recalculation, when only sequence at a site is changed
        boolean seqUpdate = nodesStates.isNodeStatesDirty(nodeNr) || nodesStates.isNodeStatesDirty(parentNum);

        int nodeUpdate = node.isDirty() | parent.isDirty();

        final double branchRate = branchRateModel.getRateForBranch(node);
        // do not use getLength, code below to save time
        final double branchTime = (parent.getHeight() - node.getHeight()) * branchRate;

//TODO deal with 0 branch length, such as SA
        if (branchTime == 0)
            throw new UnsupportedOperationException("0 branch length, such as SA, not supported !");
        if (branchTime < 1e-10)
            throw new ArithmeticException("Reject proposal : " +
                    "time is 0 at the branch between parent node " + parentNum + " and node " + nodeNr +
                    " !\n" + "branch length = " + node.getLength() + ", branchRate = " + branchRate);

        //TODO how to distinguish branch len change and internal node seq change, when topology is same

        // ====== 1. update the transition probability matrix(ices) if the branch len changes ======
        if (seqUpdate || nodeUpdate != Tree.IS_CLEAN || branchTime != branchLengths[nodeNr]) {
            this.branchLengths[nodeNr] = branchTime;
            daBranchLdCore.setNodeMatrixForUpdate(); // TODO review the index
            // rate category
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                // pass the reference of array for caching probability,
                // Note: cannot move it in this class because of multithreading by nodes.
                double[] probabilities = daBranchLdCore.getProbRef();
                double[] iexp = daBranchLdCore.getIexpRef();
                // this new code is faster
                substitutionModel.getTransiProbs(parent.getHeight(), node.getHeight(),
                        jointBranchRate, iexp, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(probabilities));

//                for (int j=0; j < probabilities.length; j++)
//                    //TODO P(t) cannot be 0, but short branch causes numeric precision error.
//                    if (probabilities[j] <= 0) {
////                        System.err.println(Arrays.toString(probabilities));
//                        throw new ArithmeticException("Select " + probabilities[j] + " probability in P(t) matrix " +
//                                ", possibly caused by a short branch ! " + "matrix index = " + j +
//                                ", branch Nr = " + nodeNr + ", branchTime = " + branchTime);
//                    }

                daBranchLdCore.setNodeMatrix(i, probabilities); //cannot rm arraycopy
            }
            nodeUpdate |= Tree.IS_DIRTY;
        }
// TODO only some sites are changed
//       else if (seqUpdate) { }

        // ====== 2. recalculate likelihood if either child node wasn't clean ======
        if (nodeUpdate != Tree.IS_CLEAN) {
            // code in SiteModel, node is not used
            final double[] proportions = siteModel.getCategoryProportions(node);
            final int[] nodeStates = nodesStates.getStates(nodeNr);
            final int[] parentNodeStates = nodesStates.getStates(parentNum);

            // brLD is linked to the child node index down
            daBranchLdCore.setBranchLdForUpdate();
            // populate branchLd[][excl. root], nodeIndex is child
            daBranchLdCore.calculateBranchLd(parentNodeStates, nodeStates, proportions);
        }

        return nodeUpdate;
    }


    //****** multi-threading ******//


    // ArithmeticException if branch time < 1e-10
    class DABranchLikelihoodCallable implements Callable<Double> {
//        private final DABranchLikelihoodCore brLDCore;
//        private final int branchNr; // used to make thread safe
//
//        // per branch
//        public DABranchLikelihoodCallable(DABranchLikelihoodCore brLDCore, int branchNr) {
//            this.brLDCore = brLDCore;
//            this.branchNr = branchNr;
//        }

        private final List<DABranchLikelihoodCore> brLDCores;

        // by group
        public DABranchLikelihoodCallable(List<DABranchLikelihoodCore> brLDCores) {
            this.brLDCores = brLDCores;
        }

        @Override
        public Double call() throws Exception {
//            try {
            double logP = 0;
            for (DABranchLikelihoodCore core : brLDCores) {
                final int branchNr = core.getBranchNr();
                final Node node = tree.getNode(branchNr);

                final int rootIndex = getRootIndex();
                if (branchNr == rootIndex) {
                    final double[] freqs = substitutionModel.getFrequencies();
                    NodeStates rootStates = nodesStates.getNodeStates(rootIndex);

                    branchLogLikelihoods[rootIndex] = core.calculateRootLogLikelihood(rootStates, freqs);
                    logP += branchLogLikelihoods[rootIndex];
                } else {
                    // caching branchLogLikelihoods[nodeNr]
                    if (updateBranch(core, node) != Tree.IS_CLEAN)
                        branchLogLikelihoods[branchNr] = core.calculateBranchLogLikelihood();

                    logP += branchLogLikelihoods[branchNr];
                }
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
//        hasDirt = Tree.IS_CLEAN;
// TODO check isDirtyCalculation()?
        if (nodesStates.somethingIsDirty()) {
//            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (siteModel.isDirtyCalculation()) {
//            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return tree.somethingIsDirty();
    }

    @Override
    public void store() {
//        if (beagle != null) {
//            beagle.store();
//            super.store();
//            return;
//        }
        super.store(); // storedLogP = logP; isDirty = false

        for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
            daBrLdCore.store();

        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, getNrOfBranches());
        System.arraycopy(branchLogLikelihoods, 0, storedBranchLogLikelihoods, 0, getNrOfBranches()+1);
    }

    @Override
    public void restore() {
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        super.restore(); // logP = storedLogP; isDirty = false

        for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
            daBrLdCore.restore();

        // pass reference in restore, but have to copy array in store.
        double[] tmp1 = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp1;

        double[] tmp2 = branchLogLikelihoods;
        branchLogLikelihoods = storedBranchLogLikelihoods;
        storedBranchLogLikelihoods = tmp2;
    }

    // for testing, nodeIndex is the child node below this branch
    public void getBranchLdFromCore(int nodeIndex, double[] branchLdOut) {
        daBranchLdCores[nodeIndex].getBranchLikelihoods(branchLdOut);
    }

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

