package toroidaldiffusion.spi;

import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;
import lphy.core.spi.LPhyCoreImpl;
import toroidaldiffusion.DihedralAngleDiffusionMatrix;
import toroidaldiffusion.PhyloCircularBrownian;
import toroidaldiffusion.PhyloToroidalBrownian;
import toroidaldiffusion.PhyloWrappedBivariateDiffusion;

import java.util.Arrays;
import java.util.List;

/**
 *
 * The provider of SPI which is an implementation of a service.
 * It requires a public no-args constructor.
 * @author Walter Xie
 */
public class ToroidalDiffusionImpl extends LPhyCoreImpl {

    /**
     * Required by ServiceLoader.
     */
    public ToroidalDiffusionImpl() {
        //TODO print package or classes info here?
    }

    @Override
    public List<Class<? extends GenerativeDistribution>> declareDistributions() {
        return Arrays.asList(PhyloCircularBrownian.class, PhyloToroidalBrownian.class,
                PhyloWrappedBivariateDiffusion.class);
    }

    @Override
    public List<Class<? extends BasicFunction>> declareFunctions() {
        return Arrays.asList(DihedralAngleDiffusionMatrix.class);
    }

    public String getExtensionName() {
        return "Toroidal diffusion";
    }
}
