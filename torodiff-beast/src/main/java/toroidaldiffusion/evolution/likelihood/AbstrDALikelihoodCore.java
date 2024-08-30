package toroidaldiffusion.evolution.likelihood;

/**
 * framework of general data augmentation likelihood core
 */
abstract public class AbstrDALikelihoodCore {

    /**
     * reserve memory for branchLd, indices and other
     * data structures required by the core *
     * call it inside constructor
     */
    protected abstract void initialize();

    /**
     * clean up after last likelihood calculation, if at all required *
     */
    @Override
    public abstract void finalize() throws Throwable;


    /**
     * store current state *
     */
    abstract public void store();

    /**
     * restore state *
     */
    abstract public void restore();

    /**
     * reset current state to stored state, only used when switching from non-scaled to scaled or vice versa *
     */
    abstract public void unstore();

}
