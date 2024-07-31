
module toroidaldiffusion.studio {

    requires transitive lphystudio;
    requires toroidaldiffusion;

    exports toroidaldiffusion.studio;

    uses toroidaldiffusion.studio.DihedralAngle3DViewer;
    uses toroidaldiffusion.studio.ToroidalPlot;

    // Viewer SPI
    uses lphystudio.app.graphicalmodelpanel.viewer.Viewer;
    // declare what service interface the provider intends to use
    provides lphystudio.app.graphicalmodelpanel.viewer.Viewer with toroidaldiffusion.studio.DihedralAnglesViewer;

    // Note: to adapt with the system not using Java module but using class path,
    // they need to be declared inside META-INF/services/lphystudio.app.graphicalmodelpanel.viewer.Viewer as well.

}