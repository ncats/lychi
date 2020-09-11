package lychi.tautomers;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.struc.Molecule;
//import chemaxon.struc.StereoConstants;
//import chemaxon.struc.Molecule;
import gov.nih.ncats.molwitch.*;
//import chemaxon.struc.MolAtom;
//import chemaxon.struc.MolBond;
import chemaxon.formats.MolImporter;
//import chemaxon.util.MolHandler;

import gov.nih.ncats.molwitch.io.ChemicalReader;
import gov.nih.ncats.molwitch.io.ChemicalReaderFactory;
import lychi.TautomerGenerator;
import lychi.util.Base32;
import lychi.util.GrayCode;

/**
 * The general description of the strategy used is described here:
 *
 *   http://tripod.nih.gov/?p=118
 * 
 * TODO: this code doesn't generate tautomers in the same order, i.e., it's
 * input dependent.
 */
public class NCGCTautomerGenerator implements TautomerGenerator, Observer {
    private static final Logger logger = 
        Logger.getLogger(NCGCTautomerGenerator.class.getName());

    static final int ACC = 1;
    static final int DON = 2;
    static final int DEFAULT_MAXCOUNT = 1024;
    static final int DEFAULT_MAXDEPTH = 10; // 1,11 is the default max

    static boolean debug = false;
    static {
        try {
            debug = Boolean.getBoolean("tautomer.debug");
        }
        catch (Exception ex) {
        }
    }

    static class Pair implements Comparable<Pair> {
        int index;
        int rank;

        Pair (int index, int rank) {
            this.index = index;
            this.rank = rank;
        }

        public int compareTo (Pair p) {
            return rank - p.rank;
        }
    }

    static class MolPath implements Comparable<MolPath> {
        int[] path;
        Set<Atom> atoms = new HashSet<Atom>();
        BitSet urank = new BitSet ();
        Chemical mol;
        Atom start, end;
        
        MolPath (Chemical mol, Bond[] path, int[] rank) {
            this.mol = mol;
            this.path = new int[path.length];

            for (int i = 0; i < path.length; ++i) {
                Bond b = path[i];

                if (start == null) {
                    start = b.getAtom1();
                    end = b.getAtom2();
                }
                else if (start == b.getAtom1()) {
                    start = end;
                    end = b.getAtom2();
                }
                else if (start == b.getAtom2()) {
                    start = end;
                    end = b.getAtom1();
                }
                else if (end == b.getAtom1()) {
                    end = b.getAtom2();
                }
                else if (end == b.getAtom2()) {
                    end = b.getAtom1();
                }
                else { // 
                    logger.warning
                        ("something's rotten in the state of maryland");
                }
                logger.info("start="+(mol.indexOf(start)+1)
                            +" end="+(mol.indexOf(end)+1)
                            +" ["+(mol.indexOf(b.getAtom1())+1)+","
                            +(mol.indexOf(b.getAtom2())+1)+"]");

                atoms.add(b.getAtom1());
                atoms.add(b.getAtom2());
                this.path[i] = mol.indexOf(b);
            }

            for (Atom a : atoms) {
                int i = mol.indexOf(a);
                if (i < 0) {
                    throw new IllegalArgumentException
                        ("Path contains atoms from different molecules!");
                }
                else {
                    urank.set(rank[i]);
                }
            }

            // order path based on the rank
            int head = Math.min(rank[mol.indexOf(path[0].getAtom1())],
                                rank[mol.indexOf(path[0].getAtom2())]);
            int tail = Math.min(rank[mol.indexOf
                                     (path[path.length-1].getAtom1())],
                                rank[mol.indexOf
                                     (path[path.length-1].getAtom2())]);
            if (tail < head) { // flip
                for (int i = 0, j = path.length-1; i < j; ++i, --j) {
                    int b = this.path[i];
                    this.path[i] = this.path[j];
                    this.path[j] = b;
                }

                // swap the two ends
                Atom a = start;
                start = end;
                end = a;
            }
        }

        public void flip (Chemical m) {
            // flip the parity for the bond designated by this path against
            // the given molecule. this implies that the input molecule must
            // be a copy 
            for (int i = 0; i < path.length; ++i) {
                Bond ref = mol.getBond(path[i]);
                m.getBond(path[i]).setBondType(ref.getBondType().switchParity());
            }
        }

        public int length () { return path.length; }
        public int[] getPath () { return path; }
        public boolean intersects (MolPath p) {
            for (Atom a : p.atoms)
                if (atoms.contains(a))
                    return true;
            return false;
        }

        // check to see if the given path has the same endpoints (starting
        // & ending) as the current path
        public boolean sameEndPoints (MolPath p) {
            return (start == p.start && end == p.end)
                || (start == p.end && end == p.start);
        }

        public int compareTo (MolPath p) {
            int i = urank.nextSetBit(0);
            int j = p.urank.nextSetBit(0);
            for (; i >= 0 && j >= 0 && i == j;) {
                i = urank.nextSetBit(i+1);
                j = p.urank.nextSetBit(j+1);
            }
            return i - j;
        }

        public String toString () {
            StringBuilder sb = new StringBuilder
                ("end: ["+(mol.indexOf(start)+1)+","
                 +(mol.indexOf(end)+1)+"] parity:");
            for (int i = 0; i < path.length; ++i) {
                Bond b = mol.getBond(path[i]);
                sb.append(" ["+(mol.indexOf(b.getAtom1())+1)
                          +","+b.getBondType()
                          +","+(mol.indexOf(b.getAtom2())+1)
                          +"]");
            }
            sb.append(" rank: "+urank);

            return sb.toString();
        }
    }

    private int maxcount;
    private int maxdepth = DEFAULT_MAXDEPTH;
    private List<Chemical> tautomers = new ArrayList<Chemical>();
    private Chemical canonicalTautomer;

    private List<MolPath> paths = new ArrayList<MolPath>();
    private int[] rank;
    private Chemical mol;

    public NCGCTautomerGenerator () {
        this (DEFAULT_MAXCOUNT);
    }

    public NCGCTautomerGenerator (int maxcount) {
        this.maxcount = maxcount;
    }

    public int getMaxDepth () { return maxdepth; }
    public void setMaxDepth (int maxdepth) { this.maxdepth = maxdepth; }

    public int getMaxTautomers () { return maxcount; }
    public void setMaxTautomers (int maxcount) { this.maxcount = maxcount; }

    public int generate (Chemical m) {
        if (!instrument (m)) {
            return -1;
        }
        tautomers.clear();

        // now generate
        GrayCode gc = GrayCode.createBinaryGrayCode(paths.size());
        gc.addObserver(this);
        gc.generate();

        return tautomers.size();
    }

    @Override
    public int generate(Molecule mol) {
        Chemical m = null;
        try {
            m = Chemical.parse(mol.toFormat("smiles"));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return generate(m);
    }

    public Enumeration<Chemical> tautomers () {
        return Collections.enumeration(tautomers);
    }

    public int getTautomerCount () { return tautomers.size(); }

    public Chemical getCanonicalTautomer () {
        throw new UnsupportedOperationException 
            ("Canonical tautomer is currently not support by "
             +"this implementation!");
    }

    @Override
    public Molecule getCanonicalTautomerRefactor() {
        throw new UnsupportedOperationException
                ("Canonical tautomer refactor is currently not support by "
                        +"this implementation!");
    }

    public void update (Observable o, Object arg) {
        int[] c = (int[])arg;
        List<MolPath> select = new ArrayList<MolPath>();
        for (int i = 0; i < c.length; ++i) {
            if (c[i] != 0) {
                MolPath pi = paths.get(i);
                for (int j = i+1; j < c.length; ++j) {
                    if (c[j] != 0) {
                        MolPath pj = paths.get(j);
                        if (pi.intersects(pj))
                            // overlapping paths, do nothing
                            return;
                    }
                }
                select.add(pi);
            }
        }

        if (tautomers.size() < maxcount) {
            // now for each path we toggle its parity
            Chemical m = mol.copy();

            if (debug) {
                StringBuilder sb = new StringBuilder ();
                for (MolPath p : select) {
                    p.flip(m);
                    sb.append(p.toString()+"\n");
                }
                m.setProperty("PARITY", sb.toString());
                try {
                    System.out.print(m.toSd());
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            else {
                for (MolPath p : select) {
                    p.flip(m);
                }
            }
            tautomers.add(m);
        }
        else { // we're done
            o.deleteObserver(this);
        }
    }

    protected boolean instrument (Chemical m) {
        mol = m.copy();
	mol.expandSGroups();
	//mol.hydrogenize(false);
	mol.makeHydrogensImplicit();//mol.implicitizeHydrogens(MolAtom.ALL_H);
        mol.kekulize();//mol.aromatize(false);

        rank = getTopologyRank (mol);
        Atom[] atoms = new Atom[mol.getAtomCount()];
        for (int i=0; i<atoms.length; i++)
            atoms[i] = mol.getAtom(i);

        BitSet don = new BitSet (atoms.length);
        BitSet acc = new BitSet (atoms.length);
        for (int i = 0; i < atoms.length; ++i) {
            int q = atoms[i].getCharge();
            int h = atoms[i].getImplicitHCount();

            int sb = h, db = 0, tb = 0; // single-, double-, triple-bond count
            for (Bond b: atoms[i].getBonds()) {
            //for (int j = 0; j < atoms[i].getBondCount(); ++j) {
                switch (b.getBondType().getOrder()) {
                case 1: ++sb; break;
                case 2: ++db; break;
                case 3: ++tb; break;
                default: // can't handle query/aromatic bond
                    return false;
                }
            }

            switch (atoms[i].getAtomicNumber()) {
            case 6: // C
                // don't consider keto-enol at the moment
                break;

            case 7: // N
                if (q == 0) {
                    if (sb == 1 && db == 1 && tb == 0) {
                        acc.set(i);
                    }
                    else if (sb == 3 && db == 0 && tb == 0) {
                        if (h > 0) {
                            don.set(i);
                        }
                    }
                    else if (sb == 1 && db == 2 && tb == 0) {
                        // XN(=O)(=O)
                    }
                }
                else if (q == 1) {
                }
                else if (q == -1) {
                    if (sb == 2 && db == 0 && tb == 0) {
                        atoms[i].setImplicitHCount(1);
                        atoms[i].setCharge(0);
                        don.set(i);
                    }
                }
                break;
                
            case  8: // O
            case 16: // S
                if (q == 0) {
                    if (sb == 2 && db == 0 && tb == 0) {
                        if (h > 0) {
                            don.set(i);
                        }
                    }
                    else if (sb == 0 && db == 1 && tb == 0) {
                        acc.set(i);
                    }
                }
                break;
            }
        }

        // generate all parity paths
        for (int i = acc.nextSetBit(0); i >= 0; i = acc.nextSetBit(i+1))
            for (int j = don.nextSetBit(0); j >= 0; j = don.nextSetBit(j+1))
                genParityPath (i, j);

        //Collections.sort(paths);

        /*
        if (debug) {
            mol.setProperty("ACC", acc.toString());
            mol.setProperty("DON", don.toString());
            System.out.print(mol.toFormat("sdf"));

            for (MolPath p : paths) {
                System.out.println(p);
            }
        }
        */

        return !paths.isEmpty();
    }

    protected void genParityPath (int start, int end) {
        LinkedList<Bond> path = new LinkedList<Bond>();
        Set<Atom> visited = new HashSet<Atom>();
        genParityPath (path, visited, 1, 
                       mol.getAtom(start), mol.getAtom(end));
    }

    protected void genParityPath 
        (LinkedList<Bond> path, Set<Atom> visited,
         int parity, Atom a, Atom end) {

        if (a == end) {
            MolPath newpath = new MolPath (mol, path.toArray(new Bond[0]), rank);
            paths.add(newpath);
        }
        else if (path.size() <= maxdepth) {
            visited.add(a);
            for (Bond b: a.getBonds()) {
                Atom xa = b.getOtherAtom(a);
                if (visited.contains(xa) || xa.getCharge() != 0) {
                }
                else if (parity + 1 == b.getBondType().getOrder()) {
                    path.push(b);
                    genParityPath (path, visited, 1-parity, xa, end);
                    path.pop();
                }
            }
            visited.remove(a);
        }
    }

    protected static int[] getTopologyRank (Chemical mol) {
        Chemical m = mol.copy();

        for (Iterator<Atom> ai = m.getAtoms().iterator(); ai.hasNext();) {
            Atom a = ai.next();
            a.setAtomicNumber(6);
            a.setRadical(0);
            a.setCharge(0);
            a.setChirality(Chirality.Non_Chiral);
	    //a.setAtno(6);
	    //a.setRadical(0);
	    //a.setCharge(0);
	    //a.setFlags(0);
	}

	for (Iterator<Bond> bi = m.getBonds().iterator(); bi.hasNext();) {
	    Bond b = bi.next();
	    b.setBondType(Bond.BondType.SINGLE);
	    b.setStereo(Bond.Stereo.NONE);
	    //b.setFlags(0);
	    //b.setType(1);
	}

	    GraphInvariant gi = m.getGraphInvariant();
	    int[] rank = TautomerGenerator.getGrinv(gi);
        return rank;
    }

    static void test () throws Exception {
        // guanine
        Chemical m = Chemical.createFromSmiles("Nc1nc2[nH]cnc2c(=O)[nH]1");
        TautomerGenerator tau = new NCGCTautomerGenerator ();
        int n = tau.generate(m);
        logger.info("guanine has "+n+" tautomers...");
        for (Enumeration<Chemical> en = tau.tautomers();
             en.hasMoreElements(); ) {
            Chemical mol = en.nextElement();
            System.out.println(mol.toSmiles());
        }
    }

    public static void main (String[] argv) throws Exception {
        TautomerGenerator tau = new NCGCTautomerGenerator ();

        if (argv.length == 0) {
            logger.info("## Reading from STDIN");
            ChemicalReader cr = ChemicalReaderFactory.newReader (System.in);
            while (cr.canRead()) {
                tau.generate(cr.read());
            }
        }
        else {
            test();
/*
            for (String a : argv) {
                Chemical mol = Chemical.parse(a);
                tau.generate(mol);
            }
*/
        }
    }
}
