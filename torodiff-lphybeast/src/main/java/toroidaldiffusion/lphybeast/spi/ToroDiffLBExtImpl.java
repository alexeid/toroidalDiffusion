package toroidaldiffusion.lphybeast.spi;

import beast.base.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.base.function.tree.InternalNodesID;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;
import toroidaldiffusion.lphybeast.tobeast.generator.PhyloWrappedBivariateDiffusionToBeast;
import toroidaldiffusion.lphybeast.tobeast.values.DihedralAnglesToBeast;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The "Container" provider class of SPI
 * which include a list of {@link ValueToBEAST},
 * {@link GeneratorToBEAST}, and {@link DataType}
 * to extend.
 * @author Walter Xie
 */
public class ToroDiffLBExtImpl implements LPhyBEASTExt {

    // the first matching converter is used.
    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return Arrays.asList( DihedralAnglesToBeast.class );
    }

    // the first matching converter is used.
    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList( PhyloWrappedBivariateDiffusionToBeast.class );
    }

    // LPhy SequenceType => BEAST DataType
    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        return new ConcurrentHashMap<>();
    }

    //*** these below are extra from Exclusion, only implemented in extensions ***//

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return List.of(InternalNodesID.class);
    }

    @Override
    public List<Class> getExcludedValueType() {
        // For a complex logic, or arrays, use isExcludedValue
        return List.of( );
    }


}
