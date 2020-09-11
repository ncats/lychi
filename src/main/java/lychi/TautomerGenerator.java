package lychi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import chemaxon.struc.Molecule;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.GraphInvariant;
import lychi.util.Base32;

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
    public int generate (Chemical mol);
    public int generate (Molecule mol);

    /**
     * Return generated tautomers
     */
    public Enumeration<Chemical> tautomers ();

    /**
     * Return the canonical tautomer if the implementation supports 
     * such feature
     */
    public Chemical getCanonicalTautomer ();
    public Molecule getCanonicalTautomerRefactor ();

    static int[] getGrinv(GraphInvariant gi) {
        ArrayList<Long> giv = new ArrayList<Long>();
        for (int i=0; i<gi.getNumberOfAtoms(); i++)
            giv.add(gi.getAtomInvariantValue(i));
        Collections.sort(giv);

        int[] rank = new int[gi.getNumberOfAtoms()];
        for (int i=0; i<rank.length; i++)
            for (int j=0; j<giv.size(); j++) {
                if (gi.getAtomInvariantValue(i) == giv.get(j)) {
                    rank[i] = j;
                    break;
                }
            }
        return rank;
    }

    /*
    private static ByteBuffer buffer = ByteBuffer.allocate(Long.BYTES);
    private static MessageDigest md;
    static {
        try {
            md = MessageDigest.getInstance("SHA1");
        } catch (NoSuchAlgorithmException e) {
            e.printStackTrace();
        }
    }

    static String getGrinvHash(GraphInvariant gi) {
        md.reset();
        ArrayList<Long> giv = new ArrayList<Long>();
        for (int i=0; i<gi.getNumberOfAtoms(); i++)
            giv.add(gi.getAtomInvariantValue(i));
        Collections.sort(giv);
        for (int i=0; i<giv.size(); i++)
            md.update(giv.get(i).byteValue());
        String h = Base32.encode45(md.digest());
        return h;
    }*/
}
