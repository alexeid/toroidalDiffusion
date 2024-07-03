package toroidalDiffusion.spi;

import jebl.evolution.sequences.SequenceType;
import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;
import lphy.core.spi.LPhyCoreImpl;
import toroidalDiffusion.DihedralAngleDiffusionMatrix;
import toroidalDiffusion.PhyloToroidalBrownian;
import toroidalDiffusion.PhyloWrappedBivariateDiffusion;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * The provider of SPI which is an implementation of a service.
 * It requires a public no-args constructor.
 * @author Walter Xie
 */
public class ToroidalDiffusion extends LPhyCoreImpl {

    /**
     * Required by ServiceLoader.
     */
    public ToroidalDiffusion() {
        //TODO print package or classes info here?
    }

    @Override
    public List<Class<? extends GenerativeDistribution>> declareDistributions() {
        return Arrays.asList(PhyloToroidalBrownian.class,
                PhyloWrappedBivariateDiffusion.class);
    }

    @Override
    public List<Class<? extends BasicFunction>> declareFunctions() {
        return Collections.singletonList(DihedralAngleDiffusionMatrix.class);
    }

    public String getExtensionName() {
        return "Toroidal diffusion";
    }
}
