package lychi;

import java.util.Enumeration;
import chemaxon.struc.Molecule;

/**
 * Generic interface for generating tautomers
 */
public interface TautomerGenerator {
    // Return the number of tautomer generated
    public int getTautomerCount ();
    // Return the maximum number of tautomers to generate
    public int getMaxTautomers (); 
    // Set the maximum number of tautomers to generate
    public void setMaxTautomers (int maxcount);

    /**
     * Generate tautomers and return the actual number of
     * tautomers generated. If the returned value is < 0, then
     * the input molecule contains features that aren't 
     * compatible (e.g., query atoms/bonds, ambiguous aromatic
     * bonds, etc.)
     */
    public int generate (Molecule mol);

    /**
     * Return generated tautomers
     */
    public Enumeration<Molecule> tautomers ();

    /**
     * Return the canonical tautomer if the implementation supports 
     * such feature
     */
    public Molecule getCanonicalTautomer ();
}
