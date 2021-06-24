package toroidalDiffusion.spi;

import jebl.evolution.sequences.SequenceType;
import lphy.graphicalModel.Func;
import lphy.graphicalModel.GenerativeDistribution;
import lphy.spi.LPhyExtension;
import toroidalDiffusion.DihedralAngleDiffusionMatrix;
import toroidalDiffusion.PhyloToroidalBrownian;
import toroidalDiffusion.PhyloWrappedBivariateDiffusion;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * The provider of SPI which is an implementation of a service.
 * It requires a public no-args constructor.
 * @author Walter Xie
 */
public class ToroidalDiffusion implements LPhyExtension {

    /**
     * Required by ServiceLoader.
     */
    public ToroidalDiffusion() {
        //TODO print package or classes info here?
    }

    @Override
    public List<Class<? extends GenerativeDistribution>> getDistributions() {
        return Arrays.asList(PhyloToroidalBrownian.class,
                PhyloWrappedBivariateDiffusion.class);
    }

    @Override
    public List<Class<? extends Func>> getFunctions() {
        return Collections.singletonList(DihedralAngleDiffusionMatrix.class);
    }

    @Override
    public Map<String, ? extends SequenceType> getSequenceTypes() {
        return new ConcurrentHashMap<>();
    }
}
