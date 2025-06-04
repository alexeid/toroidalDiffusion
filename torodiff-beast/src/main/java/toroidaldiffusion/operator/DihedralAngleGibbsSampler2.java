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
    public double logHastingsratio;
    public double varianceInflation = 100;
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

            // in original frame of reference
            Pair proposedAnglePair = sample(combinedDistParam, offset);

            newAngles[phiIndex] = proposedAnglePair.getPhi();
            newAngles[psiIndex] = proposedAnglePair.getPsi();

            // Calculate density of the proposed angles under the combined distribution
            logForwardDensity += logBivariateNormalWrappedPDF(proposedAnglePair.getPhi(), proposedAnglePair.getPsi(), offset, combinedDistParam);
            // Calculate density of the current angles under the combined distribution
            logBackwardDensity += logBivariateNormalWrappedPDF(nodeValues[phiIndex], nodeValues[psiIndex], offset, combinedDistParam);

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
     * This is used for Gibbs-like sampling at an internal node of a tree,
     * combining bivariate distribution from a parent node and two child nodes.
     * e.g.
     * @param parent = a mean vector(μ) of biviriate distribution: [meanPhi, meanPsi]
     * @param parentcovMatrix =  a covariance matrix (Σ) containing [invVarPsi, invVarPhi, invCovar]
     * Combining Precision Matrices Σ-1 (Ω): Ω_Combined = Ω_Parent + Ω_Ch1 + Ω_Ch2
     * Combining Precision * Mean (Ω * μ): Ω_Combined * μ_Combined = Ω_Parent * μ_Parent + Ω_Ch1 * μ_Ch1 + Ω_Ch2 * μ_Ch2
     * Recover combined covariance matrix: Σ_Combined = Ω_Combined^-1
     * Recover combined mean: μ_Combined = Ω_Combined * μ_Combined * Ω_Combined^-1
     * @return
     */
    public BivariateNormalDistParam combineBivariateNormals(Pair parent, CovarianceMatrix parentcovMatrix, Pair child1, CovarianceMatrix ch1covMatrix, Pair child2, CovarianceMatrix ch2covMatrix) {

        boolean isRoot = parent == null || parentcovMatrix == null;

        double[] precMatrixParent;
        double[] prectimesMeanParent;

        //For parent
        //IF isRoot = TRUE, covariance matrix = 0, mean vector = 0
        if (isRoot){
            precMatrixParent = new double[]{0, 0, 0};
            prectimesMeanParent = new double[]{0, 0, 0};
        } else {
            // Calculate precision matrices (Σ-1)
            precMatrixParent = calculatePrecisionMatrix(parentcovMatrix.getVarPhi(), parentcovMatrix.getVarPsi(), parentcovMatrix.getCovar());
            //Calculate Precision * Mean
            prectimesMeanParent = calculatePrecisionTimesMean(precMatrixParent, parent);
        }

        // Calculate precision matrices for child1 + child2
        double[] precMatrixCh1 = calculatePrecisionMatrix(ch1covMatrix.getVarPhi(), ch1covMatrix.getVarPsi(), ch1covMatrix.getCovar());
        double[] precMatrixCh2 = calculatePrecisionMatrix(ch2covMatrix.getVarPhi(), ch2covMatrix.getVarPsi(), ch2covMatrix.getCovar());

        // Calculate precision * mean for child1 + child2
        double[] prectimesMeanCh1 = calculatePrecisionTimesMean(precMatrixCh1, child1);
        double[] prectimesMeanCh2 = calculatePrecisionTimesMean(precMatrixCh2, child2);

//        //invVarPsi, invVarPhi, invCovar
//        double invVarPsiParent = precMatrixParent[0];
//        double invVarPhiParent = precMatrixParent[1];
//        double invCovarParent = precMatrixParent[2];

        //For child1
//        double invVarPsiCh1 = precMatrixCh1[0];
//        double invVarPhiCh1 = precMatrixCh1[1];
//        double invCovarCh1 = precMatrixCh1[2];

        //For child2
//        double invVarPsiCh2 = precMatrixCh2[0];
//        double invVarPhiCh2 = precMatrixCh2[1];
//        double invCovarCh2 = precMatrixCh2[2];

        // Sum these products
        double sumPrectimesMean1 = prectimesMeanParent[0] + prectimesMeanCh1[0] + prectimesMeanCh2[0];
        double sumPrectimesMean2 = prectimesMeanParent[1] + prectimesMeanCh1[1] +  prectimesMeanCh2[1];

        // Combine precision matrices (add them)
        double combInvVarPsi = precMatrixParent[0] + precMatrixCh1[0] + precMatrixCh2[0];
        double combInvVarPhi = precMatrixParent[1] + precMatrixCh1[1] + precMatrixCh2[1];
        double combInvCovar = precMatrixParent[2] + precMatrixCh1[2] + precMatrixCh2[2];

        // Recover combined covariance matrix (inverse of combined precision matrix)
        double comb_det = combInvVarPsi * combInvVarPhi - combInvCovar * combInvCovar;  //det(Ω)=v1*v2 − cv2
//        double combVarPhi = combInvVarPhi / comb_det;
//        double combVarPsi = combInvVarPsi / comb_det;
        double combVarPhi = (combInvVarPhi / comb_det) * varianceInflation;
        double combVarPsi = (combInvVarPsi / comb_det) * varianceInflation;

        double combCovar = -combInvCovar / comb_det;

        // Calculate combined means
        double combMeanPhi = (combInvVarPhi * sumPrectimesMean1 - combInvCovar * sumPrectimesMean2) / comb_det;
        double combMeanPsi = (combInvVarPsi * sumPrectimesMean2 - combInvCovar * sumPrectimesMean1) / comb_det;

        // Return the parameters of the combined distribution
        return new BivariateNormalDistParam(combMeanPhi, combMeanPsi, combVarPhi, combVarPsi, combCovar);
    }


    /**
     * Calculates the precision matrix (inverse of covariance matrix) for a bivariate normal distribution.
     *
     * @param varPhi Variance of phi
     * @param varPsi Variance of psi
     * @param covar Covariance between phi and psi
     * Σ = [[varPhi covar]],
            [covar varPsi]]
     * Σ-1 = 1/det * [[varPsi -covar]], = [[invVarPsi invCovar]],
     *                [-covar  varPhi]]    [invCovar  invVarPhi]]
     * @return double[] Array containing [invVarPsi, invVarPhi, invCovar]
     */
    public double[] calculatePrecisionMatrix(double varPhi, double varPsi, double covar) {
        double det = varPhi * varPsi - covar * covar; // det = 1 / (varPsi*varPsi - covar^2)
        double invVarPsi = varPsi / det;
        double invVarPhi = varPhi / det;
        double invCovar = -covar / det;

        return new double[]{invVarPsi, invVarPhi, invCovar};
    }

    /**
     * Calculates the product of precision matrix and mean vector.
     *
     * @param precMatrix Σ-1(Ω) = [[invVarPsi invCovar]],
     *                             [invCovar  invVarPhi]]
     * @param mean mean vector(μ) of biviriate distribution [meanPhi, meanPsi]
     *
     * precMatrix * mean = [[invVarPsi invCovar]], * [meanPhi,
     *                      [invCovar  invVarPhi]]    meanPsi]
     * @return a vector = [prectimesMean1, prectimesMean2]
     */
    public double[] calculatePrecisionTimesMean(double[] precMatrix, Pair mean) {
        double invVarPsi = precMatrix[0], invVarPhi = precMatrix[1], invCovar = precMatrix[2];

        double prectimesMean1 = invVarPsi * mean.getPhi() + invCovar * mean.getPsi(); //(Ω * μ)1 = invVarPsi * meanPhi + invCovar * meanPsi
        double prectimesMean2 = invCovar * mean.getPhi() + invVarPhi * mean.getPsi(); //(Ω * μ)2 = invCovar * meanPhi + invVarPhi * meanPsi

        return new double[]{prectimesMean1, prectimesMean2};
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
    public Pair offsetAnglePairAndWrap(Pair anglePair, Pair offset) {
        // offset means only
        double newPhi = ToroidalUtils.wrapToMaxAngle(anglePair.getPhi() + offset.getPhi()); // offset mean1 (phi)
        double newPsi = ToroidalUtils.wrapToMaxAngle(anglePair.getPsi() + offset.getPsi()); // offset mean2 (psi)
        return new Pair(newPhi, newPsi);
    }


    public Pair getOffsetMean(Pair parent, Pair child1, Pair child2) {
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

        return offset;
    }


    /**
     *
     * @param phi original frame of reference
     * @param psi original frame of reference
     * @param combinedDistParam offset frame of reference
     * Multivariate Normal Distribution PDF: p(x)= 1/((2π)^d/2 * ∣Σ∣^0.5) * exp(−0.5(x−μ)T * Σ−1(x−μ))
     * @return the density of the given angle pair after offsetting to match distribution frame of reference
     */
    public double logBivariateNormalWrappedPDF(double phi, double psi, Pair offset, BivariateNormalDistParam combinedDistParam) {
        double meanPhi = combinedDistParam.meanPhi;
        double meanPsi = combinedDistParam.meanPsi;
        double varPhi = combinedDistParam.varPhi;
        double varPsi = combinedDistParam.varPsi;
        double covar = combinedDistParam.covar;

        Pair angle = new Pair(phi, psi);
        Pair mean = new Pair(meanPhi, meanPsi);

        // Apply rotation to make angles closest to the mean
        Pair offsetAngle = offsetAnglePairAndWrap(angle, offset);

        // Apply the same rotation to the mean
        Pair offsetMean = offsetAnglePairAndWrap(mean, offset);

        // Calculate precision matrix elements: Σ−1
        double[] precMatrix = calculatePrecisionMatrix(varPhi, varPsi, covar);
        double invVarPsi = precMatrix[0];
        double invVarPhi = precMatrix[1];
        double invCovar = precMatrix[2];
        // Calculate determinant of covariance matrix
        double det = varPhi * varPsi - covar * covar;

        // Calculate deviation from mean: (x − μ)
        double dx = offsetAngle.getPhi() - offsetMean.getPhi();
        double dy = offsetAngle.getPsi() - offsetMean.getPsi();

        // Calculate the exponent part (Mahalanobis distance): (x−μ)T * Σ−1 * (x−μ)
        // multivariate normal scalling factor = -0.5
        double exponent = -0.5 * (invVarPsi*dx*dx + 2*invCovar*dx*dy + invVarPhi*dy*dy);

        // Calculate log normalization constant: 1/((2π)^d/2 * ∣Σ∣^0.5)
        double logNormConst = -Math.log(2 * Math.PI) - 0.5 * Math.log(det);

        return logNormConst + exponent;
    }



    public double getLogHastingsratio() {
        return logHastingsratio;
    }

}