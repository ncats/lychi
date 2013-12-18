package lychi.tautomers;

import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Enumeration;
import java.util.Collections;
import java.util.Arrays;
import java.util.Comparator;
import java.util.BitSet;
import java.util.Observer;
import java.util.Observable;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.struc.StereoConstants;
import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.formats.MolImporter;
import chemaxon.util.MolHandler;

import lychi.TautomerGenerator;
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
        Set<MolAtom> atoms = new HashSet<MolAtom>();
        BitSet urank = new BitSet ();
        Molecule mol;
        MolAtom start, end;
        
        MolPath (MolBond[] path, int[] rank) {
            this.path = new int[path.length];

            for (int i = 0; i < path.length; ++i) {
                MolBond b = path[i];
                if (mol == null) {
                    mol = (Molecule)b.getParent();
                }

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

            for (MolAtom a : atoms) {
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
                MolAtom a = start;
                start = end;
                end = a;
            }
        }

        public void flip (Molecule m) {
            // flip the parity for the bond designated by this path against
            // the given molecule. this implies that the input molecule must
            // be a copy 
            for (int i = 0; i < path.length; ++i) {
                MolBond ref = mol.getBond(path[i]);
                m.getBond(path[i]).setType(3 - ref.getType());
            }
        }

        public int length () { return path.length; }
        public int[] getPath () { return path; }
        public boolean intersects (MolPath p) {
            for (MolAtom a : p.atoms)
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
                MolBond b = mol.getBond(path[i]);
                sb.append(" ["+(mol.indexOf(b.getAtom1())+1)
                          +","+b.getType()
                          +","+(mol.indexOf(b.getAtom2())+1)
                          +"]");
            }
            sb.append(" rank: "+urank);

            return sb.toString();
        }
    }

    private int maxcount;
    private int maxdepth = DEFAULT_MAXDEPTH;
    private List<Molecule> tautomers = new ArrayList<Molecule>();
    private Molecule canonicalTautomer;

    private List<MolPath> paths = new ArrayList<MolPath>();
    private int[] rank;
    private Molecule mol;

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

    public int generate (Molecule m) {
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

    public Enumeration<Molecule> tautomers () {
        return Collections.enumeration(tautomers);
    }
    public int getTautomerCount () { return tautomers.size(); }

    public Molecule getCanonicalTautomer () {
        throw new UnsupportedOperationException 
            ("Canonical tautomer is currently not support by "
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
            Molecule m = mol.cloneMolecule();

            if (debug) {
                StringBuilder sb = new StringBuilder ();
                for (MolPath p : select) {
                    p.flip(m);
                    sb.append(p.toString()+"\n");
                }
                m.setProperty("PARITY", sb.toString());
                System.out.print(m.toFormat("sdf"));
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

    protected boolean instrument (Molecule m) {
        mol = m.cloneMolecule();
	mol.expandSgroups();
	mol.hydrogenize(false);
	mol.implicitizeHydrogens(MolAtom.ALL_H);
        mol.aromatize(false);

        rank = getTopologyRank (mol);
        MolAtom[] atoms = mol.getAtomArray();

        BitSet don = new BitSet (atoms.length);
        BitSet acc = new BitSet (atoms.length);
        for (int i = 0; i < atoms.length; ++i) {
            int q = atoms[i].getCharge();
            int h = atoms[i].getImplicitHcount();

            int sb = h, db = 0, tb = 0; // single-, double-, triple-bond count
            for (int j = 0; j < atoms[i].getBondCount(); ++j) {
                MolBond b = atoms[i].getBond(j);
                switch (b.getType()) {
                case 1: ++sb; break;
                case 2: ++db; break;
                case 3: ++tb; break;
                default: // can't handle query/aromatic bond
                    return false;
                }
            }

            switch (atoms[i].getAtno()) {
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
                        atoms[i].setImplicitHcount(1);
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
        LinkedList<MolBond> path = new LinkedList<MolBond>();
        Set<MolAtom> visited = new HashSet<MolAtom>();
        genParityPath (path, visited, 1, 
                       mol.getAtom(start), mol.getAtom(end));
    }

    protected void genParityPath 
        (LinkedList<MolBond> path, Set<MolAtom> visited, 
         int parity, MolAtom a, MolAtom end) {

        if (a == end) {
            MolPath newpath = new MolPath (path.toArray(new MolBond[0]), rank);
            paths.add(newpath);
        }
        else if (path.size() <= maxdepth) {
            visited.add(a);
            for (int i = 0; i < a.getBondCount(); ++i) {
                MolBond b = a.getBond(i);
                MolAtom xa = b.getOtherAtom(a);
                if (visited.contains(xa) || xa.getCharge() != 0) {
                }
                else if (parity + 1 == b.getType()) {
                    path.push(b);
                    genParityPath (path, visited, 1-parity, xa, end);
                    path.pop();
                }
            }
            visited.remove(a);
        }
    }

    protected static int[] getTopologyRank (Molecule mol) {
        Molecule m = mol.cloneMolecule();

        for (MolAtom a : m.getAtomArray()) {
	    a.setAtno(6);
	    a.setRadical(0);
	    a.setCharge(0);
	    a.setFlags(0);
	}

	for (MolBond b : m.getBondArray()) {
	    b.setFlags(0);
	    b.setType(1);
	}

        int[] rank = new int[m.getAtomCount()];
        m.getGrinv(rank);

        return rank;
    }

    static void test () throws Exception {
        // guanine
        MolHandler mh = new MolHandler ("Nc1nc2[nH]cnc2c(=O)[nH]1");
        TautomerGenerator tau = new NCGCTautomerGenerator ();
        int n = tau.generate(mh.getMolecule());
        logger.info("guanine has "+n+" tautomers...");
        for (Enumeration<Molecule> en = tau.tautomers();
             en.hasMoreElements(); ) {
            Molecule m = en.nextElement();
            System.out.println(m.toFormat("smiles:q"));
        }
    }

    public static void main (String[] argv) throws Exception {
        TautomerGenerator tau = new NCGCTautomerGenerator ();

        if (argv.length == 0) {
            logger.info("## Reading from STDIN");
            MolImporter mi = new MolImporter (System.in);
            for (Molecule mol = new Molecule (); mi.read(mol); ) {
                tau.generate(mol);
            }
        }
        else {
            for (String a : argv) {
                MolImporter mi = new MolImporter (a);
                for (Molecule mol = new Molecule (); mi.read(mol); ) {
                    tau.generate(mol);
                }
            }
        }
    }
}
