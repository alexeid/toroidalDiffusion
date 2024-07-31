
/**
 * @author Walter Xie
 */
module toroidaldiffusion {
    requires lphy.base;

    // TODO review, some libs seem not used
    requires ejml.core;
    requires ejml.simple;
    requires biojava.core;
    requires biojava.alignment;
    requires biojava.structure;
    requires org.glassfish.jaxb.core;

    exports toroidaldiffusion;
    exports toroidaldiffusion.spi;

    // LPhy extensions
    uses lphy.core.spi.Extension;
    // declare what service interface the provider intends to use
    provides lphy.core.spi.Extension with toroidaldiffusion.spi.ToroidalDiffusionImpl;
}