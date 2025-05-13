package toroidaldiffusion.operator;

import beast.base.evolution.tree.Node;
import org.ejml.simple.SimpleMatrix;
import toroidaldiffusion.ToroidalUtils;
import toroidaldiffusion.WrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;

import java.util.Random;

public class DihedralAngleGibbsSampler2 {

    private WrappedBivariateDiffusion diff; //calculation engine
    private double logForwardDensity;
    private double logBackwardDensity;
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
        this.logBackwardDensity = 0;
        this.logForwardDensity = 0;
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

        //branch time
        final double parentBranchTime = parent.getHeight() - node.getHeight();
        final double child1BranchTime = node.getHeight() - child1.getHeight();
        final double child2BranchTime = node.getHeight() - child2.getHeight();

        //nodes values
        double[] parentNodeValues = daTreeModel.getNodeValue(parent);
        double[] child1NodeValues = daTreeModel.getNodeValue(child1);
        double[] child2NodeValues = daTreeModel.getNodeValue(child2);


        //check each site, construct distribution
        for (int i = 0; i < daTreeModel.getSiteCount(); i++) {
            int phiIndex = i * 2;
            int psiIndex = i * 2 + 1;

            double parentPhi = parentNodeValues[phiIndex];
            double parentPsi = parentNodeValues[psiIndex];

            double c1Phi = child1NodeValues[phiIndex];
            double c1Psi = child1NodeValues[psiIndex];

            double c2Phi = child2NodeValues[phiIndex];
            double c2Psi = child2NodeValues[psiIndex];

            //For root: parentDist = null
            double[] parentDist = getBivariateNormalParameters(diff,
                    parentPhi, parentPsi, parentBranchTime);
            double[] child1Dist = getBivariateNormalParameters(diff,
                    c1Phi, c1Psi, child1BranchTime);
            double[] child2Dist = getBivariateNormalParameters(diff,
                    c2Phi, c2Psi, child2BranchTime);

            double[] combinedDistWithOffset = getRotatedCombinedDistribution(parentDist, child1Dist, child2Dist);

            double[] proposedAnglePair = sample(combinedDistWithOffset);

            newAngles[phiIndex] = proposedAnglePair[0];
            newAngles[psiIndex] = proposedAnglePair[1];

            // Calculate density of the proposed angles under the combined distribution
            logForwardDensity += logBivariateNormalWrappedPDF(proposedAnglePair[0], proposedAnglePair[1], combinedDistWithOffset);
            // Calculate density of the current angles under the combined distribution
            logBackwardDensity += logBivariateNormalWrappedPDF(nodeValues[phiIndex], nodeValues[psiIndex], combinedDistWithOffset);
        }

        return newAngles;
    }

    private double[] sample(double[] combinedDistWithOffset) {

//        double[] parentDist = getBivariateNormalParameters(diff,
//                parentPhi, parentPsi, parentBranchTime);
//        double[] child1Dist = getBivariateNormalParameters(diff,
//                c1Phi, c1Psi, child1BranchTime);
//        double[] child2Dist = getBivariateNormalParameters(diff,
//                c2Phi, c2Psi, child2BranchTime);
//
//        double[] combinedDistWithOffset = getRotatedCombinedDistribution(parentDist, child1Dist, child2Dist);

        // Save the combined distribution (without phiOffset and psiOffset)

        // Extract the offset for sampling
        double offsetPhi = combinedDistWithOffset[5];
        double offsetPsi = combinedDistWithOffset[6];

        //sample one pair of angles
        double[] sample = sampleFromBivariateNormal(combinedDistWithOffset);

        // Rotate back
        double proposedPhi = ToroidalUtils.wrapToMaxAngle((sample[0] - offsetPhi)); //is it necessary - offset?
        double proposedPsi = ToroidalUtils.wrapToMaxAngle((sample[1] - offsetPsi));
        double[] proposedPair = new double[2];
        proposedPair[0] = proposedPhi;
        proposedPair[1] = proposedPsi;

        return proposedPair;
    }


    /**
     * Gets the parameters of the stationary wrapped bivariate normal distribution
     * @return Array containing [mean_phi, mean_psi, variance_phi, variance_psi, covariance]
    //     */
//    public double[] getStationaryDistributionParameters() {
//        // The mean of the stationary distribution is simply mu
//        double mean_phi = diff.mu.get(0, 0);
//        double mean_psi = diff.mu.get(1, 0);
//
//        // invSigmaA which already contains 2 * Σ⁻¹ * A
//        // invert to get A⁻¹ * Σ / 2
//        SimpleMatrix invSigmaA = diff.getInvSigmaA();
//        SimpleMatrix statCovar = invSigmaA.invert();
//
//        // Extract components
//        double var_phi = statCovar.get(0, 0);
//        double var_psi = statCovar.get(1, 1);
//        double covar = statCovar.get(0, 1);
//
//        return new double[] {
//                mean_phi,
//                mean_psi,
//                var_phi,
//                var_psi,
//                covar
//        };
//    }


    /**
     * Extracts parameters of the bivariate normal distribution
     *
     * @param diff The WrappedBivariateDiffusion instance with the model parameters
     * @param phi0 Initial phi angle
     * @param psi0 Initial psi angle
     * @param t    Time
     * @return Array containing [mean_phi, mean_psi, variance_phi, variance_psi, covariance]
     */
    public double[] getBivariateNormalParameters(WrappedBivariateDiffusion diff,
                                                 double phi0, double psi0,
                                                 double t) {
        // Update time-dependent parameters in the diffusion model
        diff.setParameters(t);

        SimpleMatrix Gammat = diff.getGammat();

        // Get covariance matrix at time t
        double var_phi = Gammat.get(0, 0);
        double var_psi = Gammat.get(1, 1);
        double covar = Gammat.get(0, 1);

        // Return all parameters as an array
        return new double[]{
                phi0,
                psi0,
                var_phi,        // variance_phi
                var_psi,        // variance_psi
                covar           // covariance
        };
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
    public double[] combineBivariateNormals(double[] parent, double[] child1, double[] child2) {
        // Extract parameters for each distribution
        // For parent
        double m1_p = parent[0];
        double m2_p = parent[1];
        double v1_p = parent[2];
        double v2_p = parent[3];
        double cv_p = parent[4];

        // For child1
        double m1_c1 = child1[0];
        double m2_c1 = child1[1];
        double v1_c1 = child1[2];
        double v2_c1 = child1[3];
        double cv_c1 = child1[4];

        // For child2
        double m1_c2 = child2[0];
        double m2_c2 = child2[1];
        double v1_c2 = child2[2];
        double v2_c2 = child2[3];
        double cv_c2 = child2[4];

        // Calculate precision matrices (inverse of covariance matrices)
        // For parent
        double[] prec_p = calculatePrecisionMatrix(v1_p, v2_p, cv_p);
        double prec11_p = prec_p[0];
        double prec22_p = prec_p[1];
        double prec12_p = prec_p[2];

        double[] prec_c1 = calculatePrecisionMatrix(v1_c1, v2_c1, cv_c1);
        double prec11_c1 = prec_c1[0];
        double prec22_c1 = prec_c1[1];
        double prec12_c1 = prec_c1[2];

        double[] prec_c2 = calculatePrecisionMatrix(v1_c2, v2_c2, cv_c2);
        double prec11_c2 = prec_c2[0];
        double prec22_c2 = prec_c2[1];
        double prec12_c2 = prec_c2[2];

        // Combine precision matrices (add them)
        double comb_prec11 = prec11_p + prec11_c1 + prec11_c2;
        double comb_prec22 = prec22_p + prec22_c1 + prec22_c2;
        double comb_prec12 = prec12_p + prec12_c1 + prec12_c2;

        // Calculate precision * mean for each distribution
        double[] mean_p = {m1_p, m2_p};
        double[] mean_c1 = {m1_c1, m2_c1};
        double[] mean_c2 = {m1_c2, m2_c2};

        double[] pm_p = calculatePrecisionTimesMean(prec_p, mean_p);
        double[] pm_c1 = calculatePrecisionTimesMean(prec_c1, mean_c1);
        double[] pm_c2 = calculatePrecisionTimesMean(prec_c2, mean_c2);

        //for parent
        double pm1_p = pm_p[0];
        double pm2_p = pm_p[1];

        //for child1
        double pm1_c1 = pm_c1[0];
        double pm2_c1 = pm_c1[1];

        //for child2
        double pm1_c2 = pm_c2[0];
        double pm2_c2 = pm_c2[1];

        // Sum these products
        double comb_pm1 = pm1_p + pm1_c1 + pm1_c2;
        double comb_pm2 = pm2_p + pm2_c1 + pm2_c2;

        // Calculate combined covariance matrix (inverse of combined precision matrix)
        double comb_det = comb_prec11 * comb_prec22 - comb_prec12 * comb_prec12;
        double comb_v1 = comb_prec22 / comb_det;
        double comb_v2 = comb_prec11 / comb_det;
        double comb_cv = -comb_prec12 / comb_det;

        // Calculate combined means
        double comb_m1 = (comb_prec22 * comb_pm1 - comb_prec12 * comb_pm2) / comb_det;
        double comb_m2 = (comb_prec11 * comb_pm2 - comb_prec12 * comb_pm1) / comb_det;

        // Return the parameters of the combined distribution
        return new double[]{comb_m1, comb_m2, comb_v1, comb_v2, comb_cv};
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
    private double[] calculatePrecisionTimesMean(double[] prec, double[] mean) {
        double prec11 = prec[0], prec22 = prec[1], prec12 = prec[2];
        double m1 = mean[0], m2 = mean[1];

        double pm1 = prec11 * m1 + prec12 * m2;
        double pm2 = prec12 * m1 + prec22 * m2;

        return new double[]{pm1, pm2};
    }

    /**
     * For a set of angle pairs {(φ₁,ψ₁), (φ₂,ψ₂), (φ₃,ψ₃)}, find offsets (δφ,δψ) such that
     * ∑ₖ,ⱼ ||(φₖ+δφ,ψₖ+δψ) - (φⱼ+δφ,ψⱼ+δψ)||² is minimised
     *
     * @param angles
     * @return bestOffsets for phi and psi
     */
    public static double[] findOptimalRotation(double[][] angles) {
        double bestDistance = Double.MAX_VALUE;
        double[] bestOffset = {0.0, 0.0};
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
                    bestOffset[0] = offset1;
                    bestOffset[1] = offset2;
                }
            }
        }

        return bestOffset;
    }

    /**
     * Rotates a bivariate normal distribution by applying offsets to the mean components.
     * This handles the circular nature of angular data.
     *
     * @param distributionParam An array containing [mean1, mean2, variance1, variance2, covariance]
     * @param offset            An array containing [offset1, offset2] to be added to the means
     * @return A new array containing the rotated distribution parameters
     */
    private double[] rotateDistribution(double[] distributionParam, double[] offset) {
        // Create a copy of the distribution array
        double[] rotatedDist = new double[distributionParam.length];

        // Copy all parameters
        System.arraycopy(distributionParam, 0, rotatedDist, 0, distributionParam.length);

        // Apply rotation to mean components only
        rotatedDist[0] = ToroidalUtils.wrapToMaxAngle(distributionParam[0] + offset[0]); // Rotate mean1 (phi)
        rotatedDist[1] = ToroidalUtils.wrapToMaxAngle(distributionParam[1] + offset[1]); // Rotate mean2 (psi)

        // Variances and covariance remain unchanged
        // rotatedDist[2] = distribution[2]; // variance1
        // rotatedDist[3] = distribution[3]; // variance2
        // rotatedDist[4] = distribution[4]; // covariance

        return rotatedDist;
    }

    public double[] sampleFromBivariateNormal(double[] params) {
        double mean1 = params[0];
        double mean2 = params[1];
        double var1 = params[2];
        double var2 = params[3];
        double covar = params[4];

        // Calculate Cholesky decomposition
        double a11 = Math.sqrt(var1);
        double a21 = covar / a11;
        double a22 = Math.sqrt(var2 - a21 * a21);

        // Generate two independent standard normal samples
        double z1 = random.nextGaussian();
        double z2 = random.nextGaussian();

        // Transform to correlated bivariate normal
        double x1 = mean1 + a11 * z1;
        double x2 = mean2 + a21 * z1 + a22 * z2;

        return new double[]{x1, x2};
    }

//    // getBivariateNormalParameters: 1. Extract all distributions as [mean1, mean2, var1, var2, covar]
//    public double[] sampleNodesAngles(double[] parentDist, double[] child1Dist, double[] child2Dist) {
//
//        // Get combined distribution with rotation information
//        double[] combinedDistWithOffset = getRotatedCombinedDistribution(parentDist, child1Dist, child2Dist);
//
//        // Extract the distribution parameters and offset
//        double[] combinedDist = new double[5]; // First 5 elements are the distribution parameters
//        System.arraycopy(combinedDistWithOffset, 0, combinedDist, 0, 5);
//
//        double offsetPhi = combinedDistWithOffset[5];
//        double offsetPsi = combinedDistWithOffset[6];
//
//        // 6. Sample from combined distribution
//        double[] sample = sampleFromBivariateNormal(combinedDist);
//
//        // 7. Rotate back
//        double phi = (sample[0] - offsetPhi) % (2*Math.PI);
//        double psi = (sample[1] - offsetPsi) % (2*Math.PI);
//
//        // Ensure positive angles
//        if (phi < 0) phi += 2*Math.PI;
//        if (psi < 0) psi += 2*Math.PI;
//
//        return new double[] {phi, psi};
//    }

    /**
     * Combines three bivariate normal distributions after finding and applying optimal rotation
     *
     * @param parentDist Parent distribution parameters
     * @param child1Dist Child 1 distribution parameters
     * @param child2Dist Child 2 distribution parameters
     * @return Combined distribution parameters after optimal rotation
     */
    public double[] getRotatedCombinedDistribution(double[] parentDist, double[] child1Dist, double[] child2Dist) {
        // Extract mean angles for optimal rotation
        double[][] meanAngles = {
                {parentDist[0], parentDist[1]},
                {child1Dist[0], child1Dist[1]},
                {child2Dist[0], child2Dist[1]}
        };

        // Find optimal rotation
        double[] offset = findOptimalRotation(meanAngles);

        // Apply rotation to distributions
        double[] rotatedParent = rotateDistribution(parentDist, offset);
        double[] rotatedChild1 = rotateDistribution(child1Dist, offset);
        double[] rotatedChild2 = rotateDistribution(child2Dist, offset);

        // Combine distributions
        double[] combinedDist = combineBivariateNormals(rotatedParent, rotatedChild1, rotatedChild2);

        // Store the rotation offset with the combined distribution
        // Append it to the end of the array for use when calculating densities
        double[] combinedDistWithOffset = new double[combinedDist.length + 2];
        System.arraycopy(combinedDist, 0, combinedDistWithOffset, 0, combinedDist.length);
        combinedDistWithOffset[combinedDist.length] = offset[0];
        combinedDistWithOffset[combinedDist.length + 1] = offset[1];

        return combinedDistWithOffset;
    }

    /**
     * Find pdf for constructed bivariate normal distribution
     * PDF of a bivariate normal distribution for a point (x₁, x₂) with means (μ₁, μ₂), variances (σ₁², σ₂²), and covariance σ₁₂:
     * f(x₁, x₂) = (1/(2π|Σ|^(1/2))) * exp(-0.5 * (x-μ)ᵀ Σ⁻¹ (x-μ))
     * Find the best offset to minimise distance between  [phi, psi] and distribution mean
     */
    public double logBivariateNormalWrappedPDF(double phi, double psi, double[] params) {
        double mean1 = params[0];
        double mean2 = params[1];
        double var1 = params[2];
        double var2 = params[3];
        double covar = params[4];

        // Create angle arrays for finding optimal rotation
        double[][] angles = {
                {phi, psi},
                {mean1, mean2}
        };

//         Find the optimal rotation
//         todo: do we need to rotate sampled angles and current angles closet to the combined distribution mean?
        double[] offset = findOptimalRotation(angles);

        // Apply rotation to make angles closest to the mean
        double rotatedPhi = ToroidalUtils.wrapToMaxAngle(phi + offset[0]);
        double rotatedPsi = ToroidalUtils.wrapToMaxAngle(psi + offset[1]);

        // Apply same rotation to the mean
        double rotatedMean1 = ToroidalUtils.wrapToMaxAngle(mean1 + offset[0]);
        double rotatedMean2 = ToroidalUtils.wrapToMaxAngle(mean2 + offset[1]);

        // Calculate precision matrix elements
        double[] prec = calculatePrecisionMatrix(var1, var2, covar);
        double prec11 = prec[0];
        double prec22 = prec[1];
        double prec12 = prec[2];

        // Calculate determinant of covariance matrix
        double det = var1 * var2 - covar * covar;

        // Calculate deviation from mean
        double dx = rotatedPhi - rotatedMean1;
        double dy = rotatedPsi - rotatedMean2;

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