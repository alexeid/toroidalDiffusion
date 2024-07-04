package toroidaldiffusion;

import lphy.base.evolution.alignment.ContinuousCharacterData;
import lphy.base.evolution.continuous.PhyloMultivariateBrownian;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

/**
 * Created by alexpopinga on 2/02/20.
 * TODO Alex P needs to check
 */
public class PhyloToroidalBrownian extends PhyloMultivariateBrownian {

    boolean anglesInRadians = true;

    // ANGLES IN RADIANS FOR THIS IMPLEMENTATIONS
    double MAX_ANGLE_VALUE = Math.PI*2.0;

    public PhyloToroidalBrownian(@ParameterInfo(name = "tree", description = "the time tree.") Value<TimeTree> tree,
                                 @ParameterInfo(name = "diffusionMatrix", description = "the multivariate diffusion rates.") Value<Double[][]> diffusionRate,
                                 @ParameterInfo(name = "y0", description = "the value of multivariate traits at the root.") Value<Double[]> y) {
        super (tree, diffusionRate, y);
    }

    protected Double[] handleBoundaries(double[] rawValues) {

        Double[] newValues = new Double[rawValues.length];

        for (int i = 0; i < rawValues.length; i++) {
            newValues[i] = ToroidalUtils.wrapToMaxAngle(rawValues[i], MAX_ANGLE_VALUE);
        }
        return newValues;
    }

    @GeneratorInfo(name = "PhyloToroidalBrownian", verbClause = "is assumed to have evolved under",
            narrativeName = "phylogenetic toroidal Brownian motion process",
            category = GeneratorCategory.PHYLO_LIKELIHOOD, examples = {"simplePhyloToroidalBrownian.lphy"},
            description = "The phylogenetic toroidal Brownian motion distribution.")
    public RandomVariable<ContinuousCharacterData> sample() {
        return super.sample();
    }
}
