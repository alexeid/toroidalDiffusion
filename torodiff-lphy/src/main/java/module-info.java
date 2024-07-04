/**
 * @author Walter Xie
 */
module toroidaldiffusion {
//    requires lphy.base;
    requires lphystudio;

    // TODO review, some libs seem not used
    requires ejml.core;
    requires ejml.simple;
    requires biojava.core;
    requires biojava.alignment;
    requires biojava.structure;
    requires jaxb.core;

    exports toroidaldiffusion;
    exports toroidaldiffusion.spi;

    // ViewerRegister loads all Viewers
//    uses lphystudio.app.graphicalmodelpanel.viewer.Viewer;
    // declare what service interface the provider intends to use

    // the core uses hard core to register all internal Viewers,
    // but extensions need to declare what service interface the provider intends to use for new Viewers.
//    provides lphystudio.app.graphicalmodelpanel.viewer.Viewer with lphystudio.viewer.PopSizeFuncViewer;

    // Note: to adapt with the system not using Java module but using class path,
    // they need to be declared inside META-INF/services/lphystudio.app.graphicalmodelpanel.viewer.Viewer as well.
}