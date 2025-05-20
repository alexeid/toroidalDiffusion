package toroidaldiffusion.operator;

import beast.base.evolution.tree.Node;
import toroidaldiffusion.BivariateNormalDistParam;
import toroidaldiffusion.Pair;
import toroidaldiffusion.ToroidalUtils;
import toroidaldiffusion.WrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;

import java.util.Random;

public class DihedralAngleGibbsSampler2 {

    private WrappedBivariateDiffusion diff; //calculation engine
//    private double logForwardDensity;
//    private double logBackwardDensity;
    public double logHastingsratio;

    private final Random random = new Random();

    /**
     * Set the diff for calculation
     *
     * @param diff The wrapped bivariate diffusion model with parameters
     */
    public void setDiff(WrappedBivariateDiffusion diff) {
        this.diff = diff;
    }


    public double[] gibbsSampling(Node node, DihedralAngleTreeModel daTreeModel) {
        int siteCount = daTreeModel.getSiteCount();
        double[] newAngles = new double[siteCount * 2];
        double logBackwardDensity = 0;
        double logForwardDensity = 0;
        this.logHastingsratio = 0;

        // 1. Verify fields are initialized
        if (node == null || diff == null) {
            throw new IllegalStateException("Sampler not properly initialized. Call update() first.");
        }

        double[] nodeValues = daTreeModel.getNodeValue(node);

        // Get 2 child nodes for both root and internalnodes
        final Node parent = node.getParent(); //if node is root return null
        final Node child1 = node.getChild(0);
        final Node child2 = node.getChild(1);

        //if node is root, node.getParent() = null
        boolean isRoot = parent == null;

        //branch time
        //final double parentBranchTime = parent.getHeight() - node.getHeight();
        final double parentBranchTime = isRoot ? 0.0 : parent.getHeight() - node.getHeight();
        final double child1BranchTime = node.getHeight() - child1.getHeight();
        final double child2BranchTime = node.getHeight() - child2.getHeight();

        // nodes values at all sites
        double[] parentNodeValues = isRoot ? null : daTreeModel.getNodeValue(parent);
        double[] child1NodeValues = daTreeModel.getNodeValue(child1);
        double[] child2NodeValues = daTreeModel.getNodeValue(child2);

        //check each site, construct distribution
        for (int i = 0; i < daTreeModel.getSiteCount(); i++) {

            int phiIndex = i * 2;
            int psiIndex = i * 2 + 1;

            Pair parentPair = isRoot ? null : new Pair(parentNodeValues[phiIndex], parentNodeValues[psiIndex]);
            Pair c1Pair = new Pair(child1NodeValues[phiIndex], child1NodeValues[psiIndex]);
            Pair c2Pair = new Pair(child2NodeValues[phiIndex], child2NodeValues[psiIndex]);

            // get offset
            Pair offset = getOffsetMean(parentPair, c1Pair, c2Pair);
            // Apply offset to mean
            Pair offsetMeanParent = isRoot ? null : offsetAnglePairAndWrap(parentPair, offset);
            Pair offsetMeanChild1 = offsetAnglePairAndWrap(c1Pair, offset);
            Pair offsetMeanChild2 = offsetAnglePairAndWrap(c2Pair, offset);

            //For root: parentDist = null
            CovarianceMatrix parentCM = new CovarianceMatrix(diff, parentBranchTime);
            CovarianceMatrix child1CM = new CovarianceMatrix(diff, child1BranchTime);
            CovarianceMatrix child2CM = new CovarianceMatrix(diff, child2BranchTime);

            BivariateNormalDistParam combinedDistParam = combineBivariateNormals(offsetMeanParent, parentCM, offsetMeanChild1, child1CM,
                    offsetMeanChild2, child2CM);

            Pair proposedAnglePair = sample(combinedDistParam, offset);

            newAngles[phiIndex] = proposedAnglePair.getPhi();
            newAngles[psiIndex] = proposedAnglePair.getPsi();

            // Calculate density of the proposed angles under the combined distribution
            logForwardDensity += logBivariateNormalWrappedPDF(proposedAnglePair.getPhi(), proposedAnglePair.getPsi(), combinedDistParam);
            // Calculate density of the current angles under the combined distribution
            logBackwardDensity += logBivariateNormalWrappedPDF(nodeValues[phiIndex], nodeValues[psiIndex], combinedDistParam);

        }
        //propose - current
        logHastingsratio = logForwardDensity - logBackwardDensity;
        // all sites
        return newAngles;
    }

    public Pair sample(BivariateNormalDistParam combinedDist, Pair offset) {

        // Calculate Cholesky decomposition
        double a11 = Math.sqrt(combinedDist.varPhi);
        double a21 = combinedDist.covar / a11;
        double a22 = Math.sqrt(combinedDist.varPsi - a21 * a21);

        // Generate two independent standard normal samples: z1, z2 ∼N(0,1)
        double z1 = random.nextGaussian();
        double z2 = random.nextGaussian();

        // Transform to correlated bivariate normal
        double sampledPhi = combinedDist.meanPhi + a11 * z1;
        double sampledPsi = combinedDist.meanPsi + a21 * z1 + a22 * z2;

        // Extract the offset for sampling
        double offsetPhi = offset.getPhi();
        double offsetPsi = offset.getPsi();

        // Rotate back
        double proposedPhi = ToroidalUtils.wrapToMaxAngle((sampledPhi - offsetPhi)); //+ or - offset?
        double proposedPsi = ToroidalUtils.wrapToMaxAngle((sampledPsi - offsetPsi));

        return new Pair(proposedPhi, proposedPsi);
    }

    /**
     * Combines three bivariate normal distributions into a single bivariate normal.
     * This is used for Gibbs-like sampling at an internal node of a phylogenetic tree,
     * combining information from a parent node and two child nodes.
     *
     * @param parent An array containing [mean1, mean2, variance1, variance2, covariance] for the parent node
     * @param child1 An array containing [mean1, mean2, variance1, variance2, covariance] for the first child node
     * @param child2 An array containing [mean1, mean2, variance1, variance2, covariance] for the second child node
     * @return double[] An array containing [mean1, mean2, variance1, variance2, covariance] for the combined distribution
     */
    public BivariateNormalDistParam combineBivariateNormals(Pair parent, CovarianceMatrix parentMat, Pair child1, CovarianceMatrix child1Mat, Pair child2, CovarianceMatrix child2Mat) {

        boolean isRoot = parent == null || parentMat == null;

        double[] prec_p;
        double[] pm_p;
        if (isRoot){
            prec_p = new double[]{0, 0, 0}; //if isRoot, covariance matrices = 0
            pm_p = new double[]{0, 0, 0};
        } else {
            // For parent
            prec_p = calculatePrecisionMatrix(parentMat.getVarPhi(), parentMat.getVarPsi(), parentMat.getCovar());
            pm_p = calculatePrecisionTimesMean(prec_p, parent);

        }

        // Calculate precision matrices (inverse of covariance matrices)
        // For parent
        //double[] prec_p = calculatePrecisionMatrix(parentMat.getVarPhi(), parentMat.getVarPsi(), parentMat.getCovar());
        double prec11_p = prec_p[0];
        double prec22_p = prec_p[1];
        double prec12_p = prec_p[2];

        //For child1
        double[] prec_c1 = calculatePrecisionMatrix(child1Mat.getVarPhi(), child1Mat.getVarPsi(), child1Mat.getCovar());
        double prec11_c1 = prec_c1[0];
        double prec22_c1 = prec_c1[1];
        double prec12_c1 = prec_c1[2];

        //For child2
        double[] prec_c2 = calculatePrecisionMatrix(child2Mat.getVarPhi(), child2Mat.getVarPsi(), child2Mat.getCovar());
        double prec11_c2 = prec_c2[0];
        double prec22_c2 = prec_c2[1];
        double prec12_c2 = prec_c2[2];

        // Calculate precision * mean for each distribution
        //double[] pm_p = calculatePrecisionTimesMean(prec_p, parent);
        double[] pm_c1 = calculatePrecisionTimesMean(prec_c1, child1);
        double[] pm_c2 = calculatePrecisionTimesMean(prec_c2, child2);

        // Sum these products
        double sum_pm1 = pm_p[0] + pm_c1[0] + pm_c2[0];
        double sum_pm2 = pm_p[1] + pm_c1[1] +  pm_c2[1];

        // Combine precision matrices (add them)
        double comb_prec11 = prec11_p + prec11_c1 + prec11_c2;
        double comb_prec22 = prec22_p + prec22_c1 + prec22_c2;
        double comb_prec12 = prec12_p + prec12_c1 + prec12_c2;

        // Calculate combined covariance matrix (inverse of combined precision matrix)
        double comb_det = comb_prec11 * comb_prec22 - comb_prec12 * comb_prec12;
        double comb_v1 = comb_prec22 / comb_det;
        double comb_v2 = comb_prec11 / comb_det;
        double comb_cv = -comb_prec12 / comb_det;

        // Calculate combined means
        double comb_m1 = (comb_prec22 * sum_pm1 - comb_prec12 * sum_pm2) / comb_det;
        double comb_m2 = (comb_prec11 * sum_pm2 - comb_prec12 * sum_pm1) / comb_det;

        // Return the parameters of the combined distribution
        return new BivariateNormalDistParam(comb_m1, comb_m2, comb_v1, comb_v2, comb_cv);
    }


    /**
     * Calculates the precision matrix (inverse of covariance matrix) for a bivariate normal distribution.
     *
     * @param v1 Variance of the first variable
     * @param v2 Variance of the second variable
     * @param cv Covariance between the two variables
     * @return double[] Array containing [precision11, precision22, precision12]
     */
    private double[] calculatePrecisionMatrix(double v1, double v2, double cv) {
        double det = v1 * v2 - cv * cv;
        double prec11 = v2 / det;
        double prec22 = v1 / det;
        double prec12 = -cv / det;

        return new double[]{prec11, prec22, prec12};
    }

    /**
     * Calculates the product of precision matrix and mean vector.
     *
     * @param prec Array containing precision matrix elements [prec11, prec22, prec12]
     * @param mean Array containing means [mean1, mean2]
     * @return double[] Array containing [prec*mean]_1, [prec*mean]_2
     */
    private double[] calculatePrecisionTimesMean(double[] prec, Pair mean) {
        double prec11 = prec[0], prec22 = prec[1], prec12 = prec[2];

        double pm1 = prec11 * mean.getPhi() + prec12 * mean.getPsi();
        double pm2 = prec12 * mean.getPhi() + prec22 * mean.getPsi();

        return new double[]{pm1, pm2};
    }

    /**
     * For a set of angle pairs {(φ₁,ψ₁), (φ₂,ψ₂), (φ₃,ψ₃)}, find offsets (δφ,δψ) such that
     * ∑ₖ,ⱼ ||(φₖ+δφ,ψₖ+δψ) - (φⱼ+δφ,ψⱼ+δψ)||² is minimised
     *
     * @param angles
     * @return bestOffsets for phi and psi
     */
    public static Pair findOptimalRotation(double[][] angles) {
        double bestDistance = Double.MAX_VALUE;
        double bestOffset1 = 0;
        double bestOffset2 = 0;
        double PI2 = 2 * Math.PI;

        // Try these offsets (in radians)
        double[] possibleOffsets = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4, Math.PI,
                5 * Math.PI / 4, 3 * Math.PI / 2, 7 * Math.PI / 4};

        // Try all combinations of offsets
        for (double offset1 : possibleOffsets) {
            for (double offset2 : possibleOffsets) {
                // Apply offset and measure total distance
                double[][] rotatedAngles = new double[angles.length][2];
                for (int i = 0; i < angles.length; i++) {
                    rotatedAngles[i][0] = (angles[i][0] + offset1) % PI2;
                    rotatedAngles[i][1] = (angles[i][1] + offset2) % PI2;
                }

                // Calculate total pairwise distances
                double totalDistance = 0;
                for (int i = 0; i < rotatedAngles.length; i++) {
                    for (int j = i + 1; j < rotatedAngles.length; j++) {
                        double dx = rotatedAngles[i][0] - rotatedAngles[j][0];
                        double dy = rotatedAngles[i][1] - rotatedAngles[j][1];
                        totalDistance += dx * dx + dy * dy;
                    }
                }

                // Keep track of best offset
                if (totalDistance < bestDistance) {
                    bestDistance = totalDistance;
                    bestOffset1 = offset1;
                    bestOffset2 = offset2;
                }
            }
        }

        return new Pair(bestOffset1, bestOffset2); // a pair of offsets for phi and psi
    }

    /**
     * Rotates a bivariate normal distribution by applying offsets to the anglePair components.
     * This handles the circular nature of angular data.
     *
     * @param offset            An array containing [offset1, offset2] to be added to the means
     * @return A new array containing the rotated distribution parameters
     */
    private Pair offsetAnglePairAndWrap(Pair anglePair, Pair offset) {
        // offset means only
        double newPhi = ToroidalUtils.wrapToMaxAngle(anglePair.getPhi() + offset.getPhi()); // offset mean1 (phi)
        double newPsi = ToroidalUtils.wrapToMaxAngle(anglePair.getPsi() + offset.getPsi()); // offset mean2 (psi)
        return new Pair(newPhi, newPsi);
    }

    public Pair sampleFromBivariateNormal(Pair mean, CovarianceMatrix covMat) {
        // Calculate Cholesky decomposition
        double a11 = Math.sqrt(covMat.getVarPhi());
        double a21 = covMat.getCovar() / a11;
        double a22 = Math.sqrt(covMat.getVarPsi() - a21 * a21);

        // Generate two independent standard normal samples
        double z1 = random.nextGaussian();
        double z2 = random.nextGaussian();

        // Transform to correlated bivariate normal
        double x1 = mean.getPhi() + a11 * z1;
        double x2 = mean.getPsi() + a21 * z1 + a22 * z2;

        return new Pair(x1, x2);
    }


    private Pair getOffsetMean(Pair parent, Pair child1, Pair child2) {
        boolean isRoot = parent == null;

        double[][] meanAngles;
        Pair offset;

        if (isRoot) {
            meanAngles = new double[][]{
                    {child1.getPhi(), child1.getPsi()},
                    {child2.getPhi(), child2.getPsi()}
            };
            offset = findOptimalRotation(meanAngles);
        } else{
            meanAngles = new double[][]{
                    {parent.getPhi(), parent.getPsi()},
                    {child1.getPhi(), child1.getPsi()},
                    {child2.getPhi(), child2.getPsi()}
            };
            offset = findOptimalRotation(meanAngles);
        }
//        // Extract mean angles for optimal rotation
//        double[][] meanAngles = {
//                {parent.getPhi(), parent.getPsi()},
//                {child1.getPhi(), child1.getPsi()},
//                {child2.getPhi(), child2.getPsi()}
//        };

        // Find optimal offset for phi and psi
        return offset;
    }

    public double logBivariateNormalWrappedPDF(double phi, double psi, BivariateNormalDistParam combinedDistParam) {
        double meanPhi = combinedDistParam.meanPhi;
        double meanPsi = combinedDistParam.meanPsi;
        double varPhi = combinedDistParam.varPhi;
        double varPsi = combinedDistParam.varPsi;
        double covar = combinedDistParam.covar;

        Pair angle = new Pair(phi, psi);
        Pair mean = new Pair(meanPhi, meanPsi);

        // Create angle arrays for finding optimal rotation
        double[][] angles = {
                {phi, psi},
                {meanPhi, meanPsi}
        };

        // TODO bad coding to use double[][] angles
        //         todo: do we need to rotate sampled angles and current angles closet to the combined distribution mean?
        Pair offset = findOptimalRotation(angles);

        // Apply rotation to make angles closest to the mean
        Pair offsetAngle = offsetAnglePairAndWrap(angle, offset);

        // Apply same rotation to the mean
        Pair offsetMean = offsetAnglePairAndWrap(mean, offset);

        // Calculate precision matrix elements
        double[] prec = calculatePrecisionMatrix(varPhi, varPsi, covar);
        double prec11 = prec[0];
        double prec22 = prec[1];
        double prec12 = prec[2];

        // Calculate determinant of covariance matrix
        double det = varPhi * varPsi - covar * covar;

        // Calculate deviation from mean
        double dx = offsetAngle.getPhi() - offsetMean.getPhi();
        double dy = offsetAngle.getPsi() - offsetMean.getPsi();

        // Calculate the exponent part (Mahalanobis distance)
        double exponent = -0.5 * (prec11*dx*dx + 2*prec12*dx*dy + prec22*dy*dy);

        // Calculate log normalization constant
        double logNormConst = -Math.log(2 * Math.PI) - 0.5 * Math.log(det);

        return logNormConst + exponent;
    }



    public double getLogHastingsratio() {
        return logHastingsratio;
    }

}