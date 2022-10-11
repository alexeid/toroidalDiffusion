/**
 * @author Walter Xie
 */
module toroidaldiffusion.lphy {
    requires transitive lphy;

    exports toroidaldiffusion.lphy;

    // declare what service interface the provider intends to use
    provides lphy.spi.LPhyExtension with toroidaldiffusion.lphy.spi.ToroidalDiffusion;
}