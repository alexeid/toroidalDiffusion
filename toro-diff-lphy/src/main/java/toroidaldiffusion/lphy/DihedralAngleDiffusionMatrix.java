package toroidaldiffusion.lphy;

import lphy.graphicalModel.DeterministicFunction;
import lphy.graphicalModel.GeneratorInfo;
import lphy.graphicalModel.ParameterInfo;
import lphy.graphicalModel.Value;
import lphy.graphicalModel.types.DoubleArray2DValue;

import java.util.Map;

/**
 * Created by adru001 on 2/02/20.
 */
public class DihedralAngleDiffusionMatrix extends DeterministicFunction<Double[][]> {

    public static final String lengthParamName = "length";
    public static final String phiVarParamName = "phiVariance";
    public static final String psiVarParamName = "psiVariance";
    public static final String covarParamName = "covariance";

    public DihedralAngleDiffusionMatrix(@ParameterInfo(name = lengthParamName, description = "the length of the peptide backbone to model the angular diffusion of.") Value<Integer> length,
                                        @ParameterInfo(name = phiVarParamName, description = "the variance of the phi angles.") Value<Double> phiVariance,
                                        @ParameterInfo(name = psiVarParamName, description = "the variance of the psi angles.") Value<Double> psiVariance,
                                        @ParameterInfo(name = covarParamName, description = "the covariance between phi and psi angles.") Value<Double> covariance) {

        setParam(lengthParamName, length);
        setParam(phiVarParamName, phiVariance);
        setParam(psiVarParamName, psiVariance);
        setParam(covarParamName, covariance);
    }


    @GeneratorInfo(name = "dihedralAngleDiffusionMatrix", description = "This function constructs a variance covariance matrix for the neutral angular diffusion model.")
    public Value<Double[][]> apply() {

        Map<String, Value> params = getParams();
        double phiVar = (Double)params.get(phiVarParamName).value();
        double psiVar = (Double)params.get(psiVarParamName).value();
        double covar = (Double)params.get(covarParamName).value();
        int length = (Integer)params.get(lengthParamName).value();

        return new DoubleArray2DValue( constructMatrix(phiVar, psiVar, covar, length), this);
    }

    private Double[][] constructMatrix(double phiVar, double psiVar, double covar, int length) {

        int matSize = length*2;

        Double[][] matrix = new Double[matSize][matSize];

        // construct matrix
        for (int i = 0; i < matSize; i++) {
            for (int j = 0; j < matSize; j++) {
                if (i == j) {
                    if (i % 2 == 0) {
                        matrix[i][j] = phiVar;
                    } else {
                        matrix[i][j] = psiVar;
                    }
                } else if (Math.abs(i-j) == 1) {
                    matrix[i][j] = covar;
                } else matrix[i][j] = 0.0;
            }
        }
        return matrix;
    }
}
