package lychi;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import chemaxon.formats.MolImporter;
import chemaxon.struc.DPoint3;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.struc.RgMolecule;
import chemaxon.util.MolHandler;
import lychi.tautomers.SayleDelanyTautomerGenerator;
import lychi.util.Base32;
import lychi.util.ChemUtil;

/**
 * A new version of CanonicalForm but without using Reactor class;
 * Please keep a detailed log of the changes in VERSION here:
 *
 * 07.22.2010: 0x05: Didn't previously handle hydrogen isotopes.
 *     Make changes to makeHydrogensImplicit from MolAtom.ALL_H to
 *     MolAtom.ALL_H & ~MolAtom.ISOTOPE_H.  Also make changes to 
 *     hashKeyArray to preserve isotopes in the last layer.
 *
 * 09.01.2010: 0x06: Fix all standardization errors based on PubChem
 *     run.  Changed NeutralizedRules, Deprotonator, 
 *     fixed NormalizedRules 0 & 35.  
 *
 * 10.28.2010: 0x07: Fix hashKeyArray to use the new AtomIterator
 *     so that a canonical ordering is guaranteed whereas the previous 
 *     version that used ChemUtil.graphInvariantOrder on the skeleton 
 *     and atom labels to break ties doesn't work too well.
 */
public class LyChIStandardizer {
    private static Logger logger = Logger.getLogger
        (LyChIStandardizer.class.getName());

    /**
     * This static version value must be updated if any changes is made
     * to this class that would be incompatible with earlier results!!!
     */
    public static final int VERSION = 0x10;

    

    /**
     * This flag, when true, checks for "deeper" symmetry by enumerating
     * unspecified stereo forms and confirming that they 
     */
	private static final boolean DEEP_SYMMETRY = true;
	
	
    static final private boolean DEBUG;
    static final private boolean UNMEX; // apply UNM extra rules
    static {
        boolean debug = false;
        boolean extra = false;
        try {
            debug = Boolean.getBoolean("lychi-debug");
            extra = Boolean.getBoolean("lychi-unm-extra");
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        DEBUG = debug;
        UNMEX = extra;
    }


    public static String[] NeutralizedRules = new String[] {
        "[*-;$([!#6;!$([-]~[+]);!$([v4;B,Fe,Al,Au,Ga,In]);!$([v5;S,Se,Te,Sn,Si,Zr]);!$([v6;P,As,Sb,Ta]);!$([v2;Ag]);!$([v3;Hg,Zn]);!$([v2;Cl,Br]);!$([v4;Cl,Br]);!$([v6;Cl,Br]);!$([v7;Te]);!$([v8;Cl,Br])]):1]>>[*+0:1][H]"
    };

    /**
     * Please be VERY CAREFUL about changing the order of this array!
     * Any rearrangement of this array will certainly produce different
     * hash keys!  Please also update the VERSION value above if any
     * changes is made here.
     */
    public static String[] NormalizedRules = new String[] {
        "[*+;$([!#8;!#16]);H0:1]([O-:2])>>[*:1](=[O:2])",
        "[N;v4:1][O-:2]>>[N+:1][O-:2]",
        //"[N;v5:1][N;v3:2]>>[N+:1][N:2]",
        "[O+:1][N-:2]>>[O:1][N:2]",
        "[O;X1:1][#6:2]=[#7;X2;v4:3]=[*:4]>>[O:1]=[#6:2][#7:3]=[*:4]",
        //"[O:1]=[#6:2][#6:3]=[#7;X3;v5:4]>>[O:1][#6:2]=[#6:3][#7+:4]",
        "[#6][NH2+:1][O;X1-:2]>>[#6][N:1]=[O:2]",
        "[#6;X3+:1]-[#7;X3:2]>>[C:1]=[N+:2]",
        "[C;X2+:1]=[N;X2:2]>>[C:1]#[N+:2]",
        "[#6;X3-:1][N;X2+:2]#[N;X1:3]>>[C:1]=[N:2]#[#7:3]",
        "[N;X2-:1][N;X2+:2]#[N;X1:3]>>[N:1]=[N:2]#[#7:3]",
        "[N+:1][C-:2]=[O:3]>>[N:1]=[C:2]=[O:3]",
        "[#6][S!@;X3+:1]([#6])[#8-:2]>>[#6][S:1]([#6])=[O:2]",
        "[#6]C([#6])=[S+:1][#8-:2]>>[#6]C([#6])=[S:1]=[O:2]",
        "[#6][S-:1]([#6])[C+:2]>>[#6][S:1]([#6])=[C:2]",
        "[*][S+:1]([*])([#8-:2])=O>>[*][S:1]([*])(=[O:2])=O",
        "[#6][S+:1]([#6])([#8-:2])=C>>[#6][S:1]([#6])(=C)=[O:2]",
        "[S!@+:1][O-;X1:2]>>[S:1]=[O:2]",
        "[#6][P+:1]([O;X2])([O;X2])[#8-:2]>>[#6][P:1](O)(O)=[O:2]",
        "[#6][P-:1]([#6])([#6])[C+:2]>>[#6][P:1]([#6])([#6])=[C:2]",
        "[O;X2][Se+:1]([O;X2])[#8-:2]>>O[Se:1](O)=[O:2]",
        "[O;X2][Si+:1]([O;X2])[#8-:2]>>O[Si:1](O)=[O:2]",
        "[N+;X1:1]=[N-:2]>>[N:1]#[N:2]",
        "[N;v5:1](=[O:2])>>[N+:1]([O-:2])",
        //"[N;v5:1](=[O:2])(=[#6:3][O:4])([*:5])>>[N;v3:1]([O:2])([#6:3]=[O:4])([*:5])",
        "[N;v5;R0:1](=[O:3])(=[O:2])>>[N+:1](=[O:3])([O-:2])",
        "[O-:1][P+;v4:2]>>[O+0:1]=[P+0;v5:2]",
        "[N-;X1:1]=[N+;X3;v5:2]>>[N+0:1]#[N+0;v5:2]",
        "[N:1]#[N:2]=[*:3]>>[N-:1]=[N+:2]=[*:3]",
        "[O-:1][*+;X3;$([#16]([#6])[#6]):2]>>[O+0:1]=[*+0;X3;$([#16]([#6])[#6]):2]",
        "[*H4+]>>*",
        "[*H3+]>>*",
        "[*H2+]>>*",
        "[*H+]>>*",

        // fixing special bogus valence
        "[N;v4;H0]>>[N;H0+]",
        "[N++;v4:1][O:2]>>[N+:1][O-:2]",
        "[S;v5:1](=[O:2])(=[O:3])([*:4])>>[S:1]([O-:2])(=[O:3])([*:4])",
        "[#6:1](=[S;v3:2])([NH:3])>>[#6:1]([S:2])(=[N;v3:3])",
        "[O-;v2:1]>>[O;+0:1]"
        // remove water
        //"[OX2;h2]>>"
    };

    public static String[] DehydratedRules = new String[] {
        "[OX2;h2]>>"
    };
    
    private static final ThreadLocal<SMIRKS[]> NEUTRALIZER = 
        new ThreadLocal<SMIRKS[]> () {
        @Override
        protected SMIRKS[] initialValue () {
            try {
                SMIRKS[] neutralizer = new SMIRKS[NeutralizedRules.length];
                for (int i = 0; i < NeutralizedRules.length; ++i) {
                    neutralizer[i] = new SMIRKS (NeutralizedRules[i]);
                }
                return neutralizer;
            }
            catch (Exception ex) {
                logger.log(Level.SEVERE, 
                           "Can't create neutralized transforms", ex);
            }
            return null;
        }
    };

    private static final ThreadLocal<SMIRKS[]> NORMALIZER = 
        new ThreadLocal<SMIRKS[]> () {
        @Override
        protected SMIRKS[] initialValue () {
            try {
                SMIRKS[] normalizer = new SMIRKS[NormalizedRules.length];
                for (int i = 0; i < NormalizedRules.length; ++i) {
                    normalizer[i] = new SMIRKS (NormalizedRules[i]);
                }
                return normalizer;
            }
            catch (Exception ex) {
                logger.log(Level.SEVERE, 
                           "Can't create normalized transforms", ex);
            }
            return null;
        }
    };

    private static final ThreadLocal<SMIRKS[]> DEHYDRATER =
        new ThreadLocal<SMIRKS[]> () {
            @Override
            protected SMIRKS[] initialValue () {
                try {
                    SMIRKS[] dehydrater = new SMIRKS[DehydratedRules.length];
                    for (int i = 0; i < DehydratedRules.length; ++i) {
                        dehydrater[i] = new SMIRKS (DehydratedRules[i]);
                    }
                    return dehydrater;
                }
                catch (Exception ex) {
                    logger.log(Level.SEVERE, 
                               "Can't create dehydrated transforms", ex);
                }
                return null;
            }
        };
    
    // normalizer for external transform rules
    private static final String EXTERNAL_TRANSFORM_RULES =
        "resources/transform_rules.smirks";
    
    private static final ThreadLocal<SMIRKS[]> EXTERNAL_NORMALIZER = 
        new ThreadLocal<SMIRKS[]> () {
        @Override
        protected SMIRKS[] initialValue () {
            BufferedReader reader = null;
            try {
                reader = new BufferedReader(new InputStreamReader(LyChIStandardizer.class.getResourceAsStream(EXTERNAL_TRANSFORM_RULES)));
                String line;
                List<String> transforms = new ArrayList<String>();
                while((line = reader.readLine()) != null) {
                    if(line.startsWith("#")) {
                        continue;
                    }
                    transforms.add(line.split("\\s+")[0]);
                }
                SMIRKS[] normalizer = new SMIRKS[transforms.size()];
                for (int i = 0; i < transforms.size(); ++i) {
                    normalizer[i] = new SMIRKS (transforms.get(i));
                }
                return normalizer;
            }
            catch (Exception ex) {
                logger.log(Level.SEVERE, 
                           "Can't create external normalized transforms", ex);
            } finally {
                if(reader != null) {
                    try {
                        reader.close();
                    } catch (IOException e) {
                        // TODO Auto-generated catch block
                        logger.log(Level.WARNING, e.getMessage());
                    }
                }
            }
            return null;
        }
    };


    static class  MolComparator implements Comparator<Molecule> {
        public int compare (Molecule m1, Molecule m2) {
            if (m1 == null && m2 == null) return 0;
            if (m1 == null) return -1;
            if (m2 == null) return 1;
            
            int dif = m2.getAtomCount() - m1.getAtomCount();
            if (dif == 0) {
                double dm = m2.getMass() - m1.getMass();
                if (dm > 0) dif = 1;
                else if (dm < 0) dif = -1;
            }
            return dif;
        }
    }

    public static void cleanUp () {
        NEUTRALIZER.remove();
        NORMALIZER.remove();
        DEHYDRATER.remove();
        EXTERNAL_NORMALIZER.remove();
    }
    
    private SMIRKS[] neutralizer = NEUTRALIZER.get();
    private SMIRKS[] normalizer = NORMALIZER.get();
    private SMIRKS[] externalNormalizer = EXTERNAL_NORMALIZER.get();

    private TautomerGenerator taugen;
    private boolean removeSalt = true;
    private Molecule[] frags;

    public LyChIStandardizer () {
        init (new SayleDelanyTautomerGenerator ());
    }

    public LyChIStandardizer (TautomerGenerator taugen) {
        init (taugen);
    }

    public void removeSaltOrSolvent (boolean remove) { removeSalt = remove; }
    public boolean isSaltOrSolventRemoved () { return removeSalt; }

    protected void init (TautomerGenerator taugen) {
        this.taugen = taugen;
    }

    public static boolean containMetals (Molecule mol) {
        MolAtom atoms[] = mol.getAtomArray();
        for (int i = 0; i < atoms.length; ++i) {
            MolAtom a = atoms[i];
            if (ElementData.isMetal(a.getAtno())) {
                return true;
            }
        }
        return false;
    }

    protected static void handleSpecialMetals (Molecule mol) {
        MolBond bonds[] = mol.getBondArray();
        List<MolBond> delete = new ArrayList<MolBond>();
        for (int i = 0; i < bonds.length; ++i) {
            MolAtom a1 = bonds[i].getAtom1();
            MolAtom a2 = bonds[i].getAtom2();
            // handle limited covalent NaO 
            if ((a1.getAtno() == 8) && a2.getAtno() == 11) {
                a2.setCharge(a2.getCharge()+1);
                a1.setCharge(a1.getCharge()-1);         
                delete.add(bonds[i]);
            }
            else if (a1.getAtno() == 11 && (a2.getAtno() == 8)) {
                a1.setCharge(a1.getCharge()+1);
                a2.setCharge(a2.getCharge()-1);
                delete.add(bonds[i]);
            }
        }
        
        for (MolBond b : delete)
            mol.removeEdge(b);
        
        mol.valenceCheck();
    }

    protected static boolean disconnectMetals (Molecule mol) {
        MolAtom atoms[] = mol.getAtomArray();
        Set<MolBond> remove = new HashSet<MolBond>();
        Set<MolAtom> metals = new HashSet<MolAtom>();
        Set<MolAtom> others = new HashSet<MolAtom>();

        for (int i = 0; i < atoms.length; ++i) {
            MolAtom a = atoms[i];
            if (ElementData.isMetal(a.getAtno()) 
                // don't disconnect Hg, Sb, etc....
                //&& !ElementData.isSpecialMetal(a.getAtno())
                ) {
                for (int j = 0; j < a.getBondCount(); ++j) {
                    MolBond b = a.getBond(j);
                    MolAtom xa = b.getOtherAtom(a);
                    // don't break bond if it's F, Cl, Br, I, P
                    switch (xa.getAtno()) {
                    case 9:
                    case 15:
                    case 17:
                    case 35:
                    case 53:
                        // do nothing
                        break;

                    default: 
                        others.add(xa);
                        metals.add(a);
                        remove.add(b);
                    }
                }
            }
        }

        boolean containsMetal = !metals.isEmpty();
        if (containsMetal) {
            for (MolBond b : remove) {
                mol.removeEdge(b);
            }
            mol.valenceCheck();

            for (MolAtom a : metals) {
                a.setCharge(ElementData.getChargeForValence
                            (a.getAtno(), 0, a.getValence() 
                             - (a.getImplicitHcount() 
                                + a.getExplicitHcount())));
            }

            mol.valenceCheck();

            for (MolAtom a : others) {
                /*
                  int r = a.getRadical();
                  if (r > 0) {
                  a.setRadical(r-1);
                  }
                  else {
                  a.setCharge(ElementData.getChargeForValence
                  (a.getAtno(), -1, a.getValence()));
                  }
                */
                if (a.hasValenceError()) {
                    System.err.println
                        ("** atom " + (mol.indexOf(a)+1) 
                         + " has valence error!");
                }
            }
        }

        return containsMetal;
    }

    protected static MolBond[] getRingBonds (Molecule mol, int[] ring) {
        BitSet bs = new BitSet (mol.getAtomCount());
        for (int i = 0; i < ring.length; ++i) {
            bs.set(ring[i]);
        }

        List<MolBond> bonds = new ArrayList<MolBond>();
        for (MolBond b : mol.getBondArray()) {
            int a1 = mol.indexOf(b.getAtom1());
            int a2 = mol.indexOf(b.getAtom2());
            if (bs.get(a1) && bs.get(a2)) {
                bonds.add(b);
            }
        }
        return bonds.toArray(new MolBond[0]);
    }

    protected static void _atom (PrintStream ps, int index, MolAtom a) {
        ps.println("  " + String.format("%1$3d", index)
                   +"[a="+a.getAtno()
                   +",m="+a.getAtomMap()
                   +",c="+a.getCharge()
                   +",h="+a.getImplicitHcount() 
                   +",H="+a.getExplicitHcount()
                   +",r="+a.getRadical()
                   +",q="+a.getQuerystr()
                   +",l="+a.getQueryLabel()
                   +",x="+a.getExtraLabel()
                   +",v="+a.getValence()
                   +",s="+a.getSymbol()
                   +"]");
    }

    public  void deprotonate (Molecule mol) {
        Deprotonator.getInstance().apply(mol);
    }

    public void protonate (Molecule mol) {
        for (int i = 0; i < neutralizer.length; ++i) {
            SMIRKS rule = neutralizer[i];
            if (rule.transform(mol)) {
                if (DEBUG) {
                    System.err.println("++ rule " + i + ": " + rule);
                    System.err.println
                        ("=> " + mol.toFormat("smiles:q"));
                }
            }
        }

        /*
         * once [X-] is converted to XH, then we need to run
         * tautomer to get the preferred tautomeric form, e.g.,
         *  CC(=S)[O-] => CC(=S)O => CC(S)(=O)
         */
        Molecule clone = mol.cloneMolecule();
        try {
            taugen.generate(mol);
            Molecule cantau = taugen.getCanonicalTautomer();
            if (cantau != mol) {
                cantau.clonecopy(mol);
            }
        }
        catch (Exception ex) {
            logger.log(Level.WARNING, 
                       "TautomerGenerator fails: "+ex.getMessage(), ex);
            clone.clonecopy(mol);
        }
    }

    public static Molecule removeSaltOrSolvent (Molecule mol) {
        // make copy of properties...
        Map<String, Object> props = new HashMap<String, Object> ();
        for (int i = 0; i < mol.getPropertyCount(); ++i) {
            String p = mol.getPropertyKey(i);
            props.put(p, mol.getPropertyObject(p));
        }

        SaltIdentifier saltId = SaltIdentifier.getInstance();
        Molecule[] frags = mol.convertToFrags();
        for (int i = 0; i < frags.length; ++i) {
            Molecule f = frags[i];
            if (saltId.isSaltOrSolvent(f)) {
                frags[i] = null;
            }
        }
        Arrays.sort(frags, new MolComparator ());

        StringBuffer fused = new StringBuffer ();
        for (int i = 0; i < frags.length; ++i) {
            if (frags[i] != null) {
                //mol.fuse(frags[i]); // <- doesn't keep stereo
                if (fused.length() > 0) {
                    fused.append('.');
                }
                fused.append(frags[i].toFormat("smiles:q"));
            }
        }

        mol.clear();
        if (fused.length() > 0) {
            try {
                if (!MolImporter.importMol(fused.toString(), mol)) {
                    System.err.println
                        ("** fatal error: unable to import molecule: " 
                         + fused);
                }
                //mol.dearomatize();
            }
            catch (chemaxon.formats.MolFormatException ex) {
                logger.log(Level.SEVERE, 
                           "Bogus fused molecule: " + fused, ex);
            }
        }

        for (Map.Entry<String, Object> e : props.entrySet()) {
            mol.setPropertyObject(e.getKey(), e.getValue());
        }

        return mol;
    }

    public Molecule standardize (String smiles) throws Exception {
        MolHandler mh = new MolHandler (smiles);
        Molecule mol = mh.getMolecule();
        standardize (mol);
        return mol;
    }

    public static boolean containsMetal (Molecule m) {
        for (MolAtom a : m.getAtomArray()) {
            if (ElementData.isMetal(a.getAtno())) {
                return true;
            }
        }
        return false;
    }

    public static boolean isInorganic (Molecule m) {
        for (MolAtom a : m.getAtomArray())
            if (a.getAtno() == 6 && a.getImplicitHcount() > 0)
                return false;
        return true;
    }

    public static boolean isQueryMol (Molecule mol) {
        mol.ungroupSgroups();
        for (MolAtom atom : mol.getAtomArray()) {
            switch (atom.getAtno()) {
            case MolAtom.LIST:
            case MolAtom.NOTLIST:
            case MolAtom.ANY:
            case MolAtom.RGROUP:
                return true;
            }
        }

        for (MolBond bond : mol.getBondArray()) {
            switch (bond.getType()) {
            case 1: case 2: case 3: case MolBond.AROMATIC:
                break;

            default:
                return true;
            }
        }
        return false;
    }


    protected static void depthSearchEdge 
        (Vector<Integer> path, int start, BitSet visited, 
         BitSet subgraph, Molecule mol) {
        if (visited.get(start)) {
            return;
        }
        visited.set(start);
        path.add(start);
        
        MolBond b = mol.getBond(start);
        MolAtom a = b.getAtom1();
        for (int i = 0; i < a.getBondCount(); ++i) {
            int nb = mol.indexOf(a.getBond(i));
            if (subgraph.get(nb)) {
                depthSearchEdge (path, nb, visited, subgraph, mol);
            }
        }

        a = b.getAtom2();
        for (int i = 0; i < a.getBondCount(); ++i) {
            int nb = mol.indexOf(a.getBond(i));
            if (subgraph.get(nb)) {
                depthSearchEdge (path, nb, visited, subgraph, mol);
            }
        }
    }

    protected static int[] depthSearchEdge 
        (int start, BitSet subgraph, Molecule mol) {
        Vector<Integer> path = new Vector<Integer>();
        BitSet visited = new BitSet ();

        depthSearchEdge (path, start, visited, subgraph, mol);
        int[] p = new int[path.size()];
        for (int i = 0; i < p.length; ++i) {
            p[i] = path.get(i);
        }
        return p;
    }

    public static void kekulize (Molecule mol) {
        mol.dearomatize();

        MolBond[] bary = mol.getBondArray();
        int[] aroma = new int[mol.getAtomCount()];
        BitSet bonds = new BitSet (bary.length);

        for (int i = 0; i < bary.length; ++i) {
            MolBond b = bary[i];
            int type = b.getType();
            if (type == MolBond.AROMATIC) {
                MolAtom a1 = b.getAtom1();
                MolAtom a2 = b.getAtom2();
                ++aroma[mol.indexOf(a1)];
                ++aroma[mol.indexOf(a2)];
                b.setType(1); // temporary set it to single
                //b.setFlags(1, MolBond.TYPE_MASK);
                bonds.set(i);
            }
        }

        if (bonds.isEmpty()) {
            return ;
        }

        if (DEBUG) {
            System.err.println
                ("** " + mol.getName() + " contains " + bonds.cardinality() 
                 + " bogus aromatic bond(s)!");
        }

        // now attempt to kekulize the structure; do a first pass to figure
        //   which bonds are candidate for modification
        for (int i = bonds.nextSetBit(0); i >= 0; i = bonds.nextSetBit(i+1)) {
            MolBond b = mol.getBond(i);
            MolAtom a1 = b.getAtom1();
            MolAtom a2 = b.getAtom2();

            if (ElementData.checkValence
                (a1.getAtno(), a1.getCharge(),
                 a1.getValence()+1, a1.getImplicitHcount())
                && ElementData.checkValence
                (a2.getAtno(), a2.getCharge(),
                 a2.getValence()+1, a2.getImplicitHcount())) {
                if (DEBUG) {
                    System.err.println("Bond "+(i+1)+":"+(mol.indexOf(a1)+1)
                                       +"-"+(mol.indexOf(a2)+1)
                                       +" requires adjustment");
                    _atom (System.err, mol.indexOf(a1)+1, a1);
                    _atom (System.err, mol.indexOf(a2)+1, a2);
                }
            }
            else {
                bonds.clear(i);
            }
        }

        // now do second pass to try all combinations; the combination
        //  that yields the largest number of aromatized bonds is 
        //  then selected
        int best = -1;
        int[] order = {};

        for (int i = bonds.nextSetBit(0); i >= 0; i = bonds.nextSetBit(i+1)) {

            bary[i].setType(2); // seed bond
            int[] path = depthSearchEdge (i, bonds, mol);

            for (int j = 0; j < path.length; ++j) {
                if (path[j] == i) {
                    continue;
                }

                MolBond b = bary[path[j]];
                MolAtom a1 = b.getAtom1();
                MolAtom a2 = b.getAtom2();

                // now make sure both neighbors have single bonds
                boolean hasDouble = false;
                for (int k = 0; k < a1.getBondCount(); ++k) {
                    MolBond ab = a1.getBond(k);
                    MolAtom a = ab.getOtherAtom(a1);
                    if (a != a2 && aroma[mol.indexOf(a)] > 0
                        && ab.getType() == 2) {
                        hasDouble = true;
                        break;
                    }
                }

                if (!hasDouble) {
                    // now do the other atom
                    for (int k = 0; k < a2.getBondCount(); ++k) {
                        MolBond ab = a2.getBond(k);
                        MolAtom a = ab.getOtherAtom(a2);
                        if (a != a1 && aroma[mol.indexOf(a)] > 0
                            && ab.getType() == 2) {
                            hasDouble = true;
                            break;
                        }
                    }
                }

                b.setType(hasDouble ? 1 : 2);
            }

            int[] current = new int[bary.length];
            for (int j = 0; j < path.length; ++j) {
                current[path[j]] = bary[path[j]].getType();
            }

            if (DEBUG) {
                System.err.print("** path");
                for (int j = 0; j < path.length; ++j) {
                    MolBond b = bary[path[j]];
                    System.err.print(" " + (mol.indexOf(b.getAtom1())+1)
                                     +"-" + (mol.indexOf(b.getAtom2())+1)
                                     +":" + b.getType());
                }
                System.err.println();
            }

            // now aromatize
            mol.aromatize();

            // now count the number of aromatic bonds
            int naro = 0;
            for (MolBond b : bary) {
                if (b.getType() == MolBond.AROMATIC) {
                    //b.setFlags(1, MolBond.TYPE_MASK);
                    ++naro;
                }
            }

            // reset the bond type
            mol.dearomatize();
            for (MolBond b : bary) {
                if (b.getType() == MolBond.AROMATIC) {
                    b.setType(1);
                }
            }

            if (DEBUG) {
                MolAtom a1 = bary[i].getAtom1();
                MolAtom a2 = bary[i].getAtom2();
                int i1 = mol.indexOf(a1), i2 = mol.indexOf(a2);
                System.err.println
                    ("Start bond "+ (i1+1)+"-" + (i2+1) + " " + naro
                     + (naro > best ? " (best)":""));
            }

            if (naro > best) {
                order = current;
                /*
                  if (naro == arobond) {
                  // bail out early if 
                  if (DEBUG) {
                  System.out.println("** truncating search!");
                  }
                  break;
                  }
                */
                best = naro;
            }
        }
        
        // update the final bond order
        if (DEBUG) {
            System.err.print("final order:");
        }
        for (int i = 0; i < order.length; ++i) {
            if (order[i] != 0) {
                bary[i].setType(order[i]);
                //bary[i].setFlags(order[i], MolBond.TYPE_MASK);
                if (DEBUG) {
                    System.err.print
                        (" " + (mol.indexOf(bary[i].getAtom1())+1)
                         +"-" + (mol.indexOf(bary[i].getAtom2())+1)
                         +":" + bary[i].getType());
                }
            }
        }

        if (DEBUG) {
            System.err.println();
            System.out.println(mol.toFormat("smiles:q"));
        }

        ChemUtil.canonicalSMILES (mol, mol, true);
    }

    // fix a bug in jchem whereby radicals are incorrectly inferred; e.g.,
    //  OC1=[N](=O)c2ccccc2[N](=O)=C1c3ccc(Cl)cc3
    public static void fixRadical (Molecule m) {
        MolAtom[] atoms = m.getAtomArray();
        for (int i = 0; i < atoms.length; ++i) {
            MolAtom a = atoms[i];
            int r = a.getRadical();
            int c = a.getCharge();
            if (r > 0) {
                /*
                  if (ElementData.checkValence
                  (a.getAtno(), c, a.getValence(), 
                  a.getImplicitHcount() + a.getExplicitHcount())) {
                  }
                  else {
                  if (DEBUG) {
                  System.err.println
                  ("fixing charge at atom " + (i+1) +": " + (r+c));
                  }
                  a.setCharge(r+c);
                  }
                */
                a.setRadical(0);
                a.valenceCheck();
            }
        }
    }

    public static void fixValence (Molecule m) {
        MolAtom[] atoms = m.getAtomArray();

        BitSet aromaticAtoms = new BitSet (atoms.length);
        for (int i = 0; i < atoms.length; ++i) {
            MolAtom a = atoms[i];
            if (a.hasValenceError()) {
                // try adjust the charge on N
                if (a.getAtno() == 7 && 
                    ElementData.checkValence
                    (a.getAtno(), a.getCharge()+1, 
                     a.getValence(), a.getImplicitHcount())) {
                    a.setCharge(a.getCharge()+1);
                }

                if (!ElementData.checkValence
                    (a.getAtno(), a.getCharge(), 
                     a.getValence(), a.getImplicitHcount())) {
                    // ok, now try  to see if this is can be fixed with
                    //  aromatization...
                    if (DEBUG) {
                        System.err.println("Adjust bonds connecting to atom");
                        _atom (System.err, i+1, a);
                    }
                    
                    for (int j = 0; j < a.getBondCount(); ++j) {
                        MolBond b = a.getBond(j);
                        if (b.getType() == MolBond.AROMATIC) {
                            // do nothing...
                        }
                        else if (m.isRingBond(m.indexOf(b))) {
                            if (DEBUG) {
                                System.err.println
                                    ("Adjusting bond " 
                                     + (m.indexOf(b)+1)
                                     + ":"+(m.indexOf(b.getAtom1())+1)
                                     + "-"+(m.indexOf(b.getAtom2())+1)
                                     + " from type="+b.getType()
                                     + " to aromatic");
                            }
                            
                            // force all ring bonds connected to this atom
                            //  to be aromatic and let the kekulize method
                            //  fix these bond type
                            b.setType(MolBond.AROMATIC);
                            aromaticAtoms.set(m.indexOf(b.getAtom1()));
                            aromaticAtoms.set(m.indexOf(b.getAtom2()));
                        }
                    }
                }
            }
        }

        // now adjust the remainder of the bonds
        if (!aromaticAtoms.isEmpty()) {

            int[][] rings = m.getSSSR();
            for (int i = 0; i < rings.length; ++i) {
                // if any atom in this ring is aromatic, then assume
                //   the reset of the atoms are also aromatic
                boolean aromatic = false;
                for (int j = 0; !aromatic && j < rings[i].length; ++j) {
                    aromatic = aromaticAtoms.get(rings[i][j]);
                }

                if (aromatic) {
                    // force all members of this ring to be aromatic
                    for (int j = 0;j < rings[i].length; ++j) {
                        aromaticAtoms.set(rings[i][j]);
                    }
                }
            }

            // now force all ring bonds to be aromatic
            for (MolBond b : m.getBondArray()) {
                int a1 = m.indexOf(b.getAtom1());
                int a2 = m.indexOf(b.getAtom2());
                if (aromaticAtoms.get(a1) && aromaticAtoms.get(a2)) {
                    if (DEBUG) {
                        if (b.getType() != MolBond.AROMATIC) {
                            System.err.println
                                ("Adjusting bond " 
                                 + (m.indexOf(b)+1)
                                 + ":"+(m.indexOf(b.getAtom1())+1)
                                 + "-"+(m.indexOf(b.getAtom2())+1)
                                 + " from type="+b.getType()
                                 + " to aromatic");
                        }
                    }
                    b.setType(MolBond.AROMATIC);
                }
            }
            //ChemUtil.canonicalSMILES (m, m, true);
        }
    }

    /**
     * This method performs rudimentary cleanup of bond parities
     * and implicit/explicit hydrogens. It's probably unnecessarily
     * complicated!
     */
    public static void preprocessing (final Molecule m) {
        ChemUtil.sanitize(m);

        // suppress all hydrogens except for isotopes
        m.hydrogenize(false);
        m.implicitizeHydrogens(MolAtom.ALL_H & ~MolAtom.ISOTOPE_H);

        List<MolAtom> addH = new ArrayList<MolAtom>();

        int[] rank = new int[m.getAtomCount()];
        m.getGrinv(rank);

        Map<MolAtom, Integer> chirality = new LinkedHashMap<MolAtom, Integer>();
        Map<MolAtom, MolBond[]> saturated = new HashMap<MolAtom, MolBond[]>();
        for (int i = 0; i < m.getAtomCount(); ++i) {
            int chiral = m.getChirality(i);
            MolAtom a = m.getAtom(i);
            if (chiral == MolAtom.CHIRALITY_R 
                || chiral == MolAtom.CHIRALITY_S) {
                if (a.getImplicitHcount() > 0)  {
                    addH.add(a);
                }
                else {
                    // check that the parity bonds around this stereocenter
                    // if there are more than one, keep the one associate
                    // with the lowest rank. if there is another bond that
                    // is lower rank and for which no parity is assigned, 
                    // then we move the parity to that bond
                    List<MolBond> rbonds = new ArrayList<MolBond>();
                    for (int k = 0; k < a.getBondCount(); ++k) {
                        MolBond b = a.getBond(k);
                        MolAtom xa = b.getOtherAtom(a);

                        // S(=O)(*)(*) 
                        if (rbonds.isEmpty() || b.getType() != 1)
                            rbonds.add(b);
                        else if (xa.getAtno() == 1) {
                            rbonds.clear();
                            break;
                        }
                        else {
                            int r = rank[m.indexOf(xa)];
                            int j = 0;
                            for (; j < rbonds.size(); ++j) {
                                MolAtom ra = rbonds.get(j).getOtherAtom(a);
                                if (r < rank[m.indexOf(ra)]) {
                                    rbonds.add(j, b);
                                    break;
                                }
                            }

                            if (j == rbonds.size())
                                rbonds.add(b); // append
                        }
                    }

                    if (!rbonds.isEmpty()) {
                        saturated.put(a, rbonds.toArray(new MolBond[0]));
                    }
                }
                chirality.put(a, chiral);
            } 
            else if (chiral != 0) {
                // only override the chirality if there are no parity
                int parity = 0;
                
                MolAtom atom  = m.getAtom(i);
                for (int k = 0; k < atom.getBondCount(); ++k) {
                    MolBond b = atom.getBond(k);
                    int flags = b.getFlags() & MolBond.STEREO1_MASK;
                    if (flags == MolBond.UP || flags == MolBond.DOWN)
                        ++parity;
                }
                logger.warning("Stereo flag "+chiral+ " on atom "+(i+1)+
                               " (#"+m.getAtom(i).getAtno()+") with "+parity
                               +" parity defined!");
                /*              
                if (parity == 0)
                m.setChirality(i, 0);*/
            }

            int ch = a.getCharge();
            if (ch < 0) {
                boolean chargeSeparated = false;
                for (int k = 0; k < a.getBondCount(); ++k) {
                    int xch = a.getBond(k).getOtherAtom(a).getCharge();
                    if (xch > 0) {
                        chargeSeparated = true;
                        break;
                    }
                }

                // only neutralize this atom if it's not in charge 
                // separated form
                if (!chargeSeparated) 
                    a.setCharge(ch+1);
            }
        } // endfor each atom

        /*
         * if there are no coordinates in the input structure, the
         * parities on the bonds aren't perceived, so we force them
         * here.
         */
        if (m.getDim() < 2) {
            m.clean(2, null);
            //logger.info(">> "+ChemUtil.canonicalSMILES(m));
        }

        /*
         * now handle stereocenters with explicit Hs
         */
        List<MolAtom> newH = new ArrayList<MolAtom>();
        for (MolAtom a : addH) {
            MolAtom h = new MolAtom (1);
            MolBond bond = new MolBond (a, h);

            // arbitrary set the parity and we'll clean it up later
            bond.setFlags(MolBond.UP, MolBond.STEREO1_MASK);
            for (int k = 0; k < a.getBondCount(); ++k) {
                MolBond b = a.getBond(k);
                int parity = b.getFlags() & MolBond.STEREO1_MASK;
                if (parity == MolBond.UP || parity == MolBond.DOWN) {
                    bond.setFlags(parity == MolBond.UP 
                                  ? MolBond.DOWN : MolBond.DOWN, 
                                  MolBond.STEREO1_MASK);
                    break;
                }
            }

            // now reset all existing bond parities (if any)
            for (int k = 0; k < a.getBondCount(); ++k) {
                MolBond b = a.getBond(k);
                b.setFlags(0, MolBond.STEREO1_MASK);
            }

            int parity = bond.getFlags() & MolBond.STEREO1_MASK;
            if (parity == MolBond.UP || parity == MolBond.DOWN) {
                if (DEBUG) {
                    logger.info("++ adding explicit H ("+getParityFlag (parity)
                                +") at stereocenter "+(m.indexOf(a)+1));
                }
                m.add(h);
                m.add(bond);
                newH.add(h);
            }
            else {
                logger.warning("Stereocenter "+(m.indexOf(a)+1)+" has "
                               +"no parity bonds; likely due "
                               +"to lack of coordinates!");
            }
        }

        if (!newH.isEmpty())
            ChemUtil.localLayout(newH.toArray(new MolAtom[0]));

        if (DEBUG) {
            logger.info("-- "+newH.size()+" chiral Hs added\n"
                        +m.toFormat("mol"));
        }

        /*
         * now process the saturated stereocenters
         */
        for (Map.Entry<MolAtom, MolBond[]> me : saturated.entrySet()) {
            MolAtom a = me.getKey();
            MolBond[] rbonds = me.getValue();
            MolBond b = rbonds[0]; // sorted by ascending rank
            int parity = b.getFlags() & MolBond.STEREO1_MASK;
            
            if (DEBUG) {
                logger.info("Stereocenter "+(m.indexOf(a)+1)
                            +": lowest ranked neighbor "
                            +(m.indexOf(b.getOtherAtom(a))+1)
                            +" with "+getParityFlag (parity)
                            +" parity; number of neighbors is "
                            +rbonds.length);
            }
            
            boolean reset = false;
            if (parity == MolBond.UP || parity == MolBond.DOWN) {
                // sweet.. not much left to do 
                reset = true;
            }
            else {
                MolBond pbond = null;
                for (int k = 1; k < rbonds.length; ++k) {
                    MolBond rb = rbonds[k];
                    int p = rb.getFlags() & MolBond.STEREO1_MASK;
                    if (p == MolBond.UP || p == MolBond.DOWN) {
                        pbond = rb;
                        break;
                    }
                }
                
                if (pbond == null) {
                    logger.warning("There must be at least "
                                   +"one parity at stereocenter!");
                }
                else {
                    int np = countBondParity (b.getOtherAtom(a), 0);
                    // the lowest ranked neighbor doesn't have any
                    // parity so we check to see if it's possible 
                    // to set it
                    if (np == 0) { 
                        if (a == b.getAtom2()) 
                            // make sure it's pointing in the right 
                            // direction
                            b.swap();
                        b.setFlags(pbond.getFlags() & MolBond.STEREO1_MASK, 
                                   MolBond.STEREO1_MASK);
                        reset = true;
                    }
                    else { 
                        
                    }
                }
            }
            
            if (reset) {
                for (int k = 1; k < rbonds.length; ++k) {
                    b = rbonds[k];
                    // only reset the parity of bonds that point toward
                    // the stereocenter
                    if (b.getAtom1() == a) 
                        b.setFlags(0, MolBond.STEREO1_MASK);
                }
            }
        }
        
        if(DEEP_SYMMETRY){
	        try{
	               
	               Map<MolAtom,MolBond> nonChiralStereo = new LinkedHashMap<>();
	               
	               for(int k=0;k<m.getBondCount();k++){
	                   MolBond b = m.getBond(k);
	                   int parity = b.getFlags() & MolBond.STEREO1_MASK;
	                   MolAtom ma1=b.getAtom1();
	                   if(chirality.get(ma1)==null){
	                       if(parity!=0){
	                               //some other parity assigned here
	//                             if(parity==MolBond.UP)System.out.println("UP"); 
	//                             if(parity==MolBond.DOWN)System.out.println("DOWN");
	                               nonChiralStereo.put(ma1, b);
	                       }
	                   }
	               }
	               
	               if(!nonChiralStereo.isEmpty()){
	                       String igprop=m.getProperty("IGNORE_COMPLEX");
	                       
	                       if(!"true".equals(igprop)){
	                               m.setProperty("IGNORE_COMPLEX", "true");
	                               
	                               Set<int[]> rings = new HashSet<int[]>();
	                               
	                               int[][] sssr=m.getSSSR();
	                               for(MolAtom ma:nonChiralStereo.keySet()){
	                                       //need to find all atoms in the ring
	                                       int im=m.indexOf(ma);
	                                       for(int[] ir:sssr){                                             
	                                               for(int i=0;i<ir.length;i++){
	                                                       if(ir[i]==im){
	                                                               rings.add(ir);
	                                                       }
	                                               }
	                                               
	                                       }
	                               }
	                               
	                               for(int[] rr:rings){
	                                       Set<MolAtom> ratoms=Arrays.stream(rr)
	                                                                                          .mapToObj(i->m.getAtom(i))
	                                                                                          .collect(Collectors.toSet());
	                                       
	                                       MolBond[] bonds=ratoms.stream()
	                                                 .filter(a->!chirality.containsKey(a))
	                                             .flatMap(a->IntStream.range(0, a.getEdgeCount()).mapToObj(i->a.getEdge(i)))
	                                             .filter(e->!ratoms.contains(e.getNode1()) || !ratoms.contains(e.getNode2()))
	                                             .map(b->(MolBond)b)
	                                             .filter(b->b.getType()==1)
	                                             .peek(b->{
	                                                 if(ratoms.contains(b.getAtom1()))b.swap();  
	                                             })
	                                             .toArray(i->new MolBond[i]);
	                                                                               
	                                       BitSet bs = new BitSet(bonds.length*2);
	                                       for(int i=0;i<bonds.length;i++){
	                                                 MolBond b=bonds[i];
	                                                 int parity = b.getFlags() & MolBond.STEREO1_MASK;
	                                                 if(parity==MolBond.UP){
	                                                         bs.set(i*2);
	                                                 }else if(parity==MolBond.DOWN){
	                                                         bs.set(i*2+1);
	                                                 }else{
	                                                         bs.set(i*2);
	                                                         bs.set(i*2+1);
	                                                 }
	                                       }
	                                       
	                                       Set<String> allPossible = new HashSet<String>();
	                                       Set<String> currentPossible = new HashSet<String>();
	                                       
	                                       for(int i=0;i<Math.pow(2, bonds.length);i++){
	                                               BitSet onOff = new BitSet(bonds.length*2);
	                                               for(int j=0;j<bonds.length;j++){
	                                                       if((i>>j&1)==1){
	                                                               onOff.set(j*2);
	                                                               bonds[j].setFlags(MolBond.UP, MolBond.STEREO1_MASK);
	                                                       }else{
	                                                               onOff.set(j*2+1);
	                                                               bonds[j].setFlags(MolBond.DOWN, MolBond.STEREO1_MASK);
	                                                       }
	                                               }
	                                               Molecule mclone=m.cloneMolecule();
	                                               //(new LyChIStandardizer()).standardize(mclone);
	                                               String hash1=LyChIStandardizer.hashKey(mclone);
	                                               allPossible.add(hash1);
	                                               onOff.or(bs);
	                                               
	                                               if(onOff.cardinality() == bs.cardinality()){
	                                                       currentPossible.add(hash1);
	                                               }
	                                       }
	                                       if(allPossible.size()==currentPossible.size()){
	                                               for(int j=0;j<bonds.length;j++){
	                                                       bonds[j].setFlags(0, MolBond.STEREO1_MASK);
	                                               }
	                                       }else{
	                                               for(int j=0;j<bonds.length;j++){
	                                                       boolean isUp=bs.get(j*2);
	                                                       boolean isDown=bs.get(j*2+1);
	                                                       
	                                                       if(isUp && ! isDown){
	                                                               bonds[j].setFlags(MolBond.UP, MolBond.STEREO1_MASK);
	                                                       }else if(!isUp && isDown){
	                                                               bonds[j].setFlags(MolBond.DOWN, MolBond.STEREO1_MASK);
	                                                       }else{
	                                                               bonds[j].setFlags(0, MolBond.STEREO1_MASK);
	                                                       }
	                                               }
	                                       }                                       
	                               }
	                               
	                               
	                               
	                               m.setProperty("IGNORE_COMPLEX", null);
	                       }
	               }
	               
	        }catch(Exception e){
	        	logger.warning("Processing symmetry threw an error:" + e.getMessage());
	        }
        }

        /*
         * and finally adjust the wedge/hash based on the chirality
         */
        for (Map.Entry<MolAtom, Integer> me : chirality.entrySet()) {
            int i = m.indexOf(me.getKey());
            m.setChirality(i, me.getValue());
        }
        
        /*
         * Actually, do it again since sometimes there are chiral values
         * which will not be fully set by this round
         */
        for (Map.Entry<MolAtom, Integer> me : chirality.entrySet()) {
            int i = m.indexOf(me.getKey());
            m.setChirality(i, me.getValue());
        }
        
        
        for (Map.Entry<MolAtom, Integer> me : chirality.entrySet()) {
            int i = m.indexOf(me.getKey());
            int oi=me.getValue();
            int ni=m.getChirality(i);
            if(oi!=ni){
            	logger.warning("The stereo on atom " + i + " may have changed chirality during standardization!");
            }
        }

        if (DEBUG) {
            logger.info(m.toFormat("mol"));
            logger.info("++ preprocessing: "+m.toFormat("smiles:q"));
        }
    } // preprocessing

    static void cleanMolecule (int dim, DPoint3[] coords, Molecule mol) {
        MolBond[] bonds = mol.getBondArray();
        int[] flags = new int[bonds.length];
        for (int i = 0; i < bonds.length; ++i) {
            flags[i] = bonds[i].getFlags();
        }

        MolAtom[] atoms = mol.getAtomArray();
        int[] fixed = new int[atoms.length];
        for (int i = 0; i < atoms.length; ++i) {
            MolAtom a = atoms[i];
            fixed[i] = i;
            int amap = a.getAtomMap();
            if (amap > 0) {
                DPoint3 pt = coords[amap-1];
                a.setAtomMap(0);
                a.setXYZ(pt.x, pt.y, pt.z);
                //System.out.println(pt.x+" "+pt.y);
            }
        }

        mol.clean(dim, null);
        for (int i = 0; i < bonds.length; ++i) {
            int stereo = flags[i] & MolBond.STEREO_MASK;
            if (stereo != (bonds[i].getFlags()&MolBond.STEREO_MASK)) {
                System.err.println("Adjusting stereo for bond "+(i+1)+
                                   ": "+(bonds[i].getFlags()&MolBond.STEREO_MASK) + " -> " + stereo);
            }
            bonds[i].setFlags(flags[i]);
            bonds[i].setType(flags[i] & MolBond.TYPE_MASK);
        }

        //mol.valenceCheck();

        // reset E/Z stereo
        //ChemUtil.adjustEZStereo(mol);
    }

    static boolean transform (SMIRKS[] rules, Molecule mol) {
        Molecule copy = mol.cloneMolecule();
        boolean valenceError = mol.hasValenceError();
        for (SMIRKS rule : rules) {
            if (rule.transform(mol)) {
                mol.valenceCheck();
                if (mol.hasValenceError() && !valenceError) {
                    logger.warning
                        (" ++ " + mol.getName() 
                         + ": bogus structure generated by rule " 
                         + ": " + rule + " ==> "
                         + ChemUtil.canonicalSMILES (mol));
                    System.err.println("** Atom with valence error:");
                    MolAtom[] atoms = mol.getAtomArray();
                    for (int k = 0; k < atoms.length; ++k) {
                        if (atoms[k].hasValenceError()) {
                            _atom (System.err, k+1, atoms[k]);
                        }
                    }
                    copy.clonecopy(mol); // restore the original molecule
                    return false;
                }
            }
        }
        return true;
    }

    public static boolean neutralize (Molecule mol) {
        return transform (NEUTRALIZER.get(), mol);
    }
    
    public static boolean normalize (Molecule mol) {
        return transform (NORMALIZER.get(), mol);
    }

    public static boolean dehydrate (Molecule mol) {
        return transform (DEHYDRATER.get(), mol);
    }

    public boolean standardize (Molecule mol) throws Exception {
        String name = mol.getName();

        if (isQueryMol (mol)) {
            logger.log(Level.WARNING, name + " is a query molecule; "
                       +"no standardization performed!");
            return false;
        }

        boolean inputHasValenceError = mol.hasValenceError();
        if (inputHasValenceError) {
            logger.log(Level.WARNING,"Input structure " + name 
                       + " has valence error at the following atom(s):");

            for (int i = 0; i < mol.getAtomCount(); ++i) {
                MolAtom a = mol.getAtom(i);
                if (a.hasValenceError()) {
                    _atom (System.err, i+1, a);
                }
            }
            //return mol;
        }

        Molecule origmol = mol.cloneMolecule();

        // make copy of properties...
        Map<String, Object> props = new HashMap<String, Object> ();
        for (int i = 0; i < mol.getPropertyCount(); ++i) {
            String p = mol.getPropertyKey(i);
            props.put(p, mol.getPropertyObject(p));
        }

        mol.expandSgroups();

        if (DEBUG) {
            logger.info(">> INPUT: "+ChemUtil.canonicalSMILES(mol));
        }

        /*
         * standardize on explicit Hs and bond parities
         */
        preprocessing (mol);

        MolAtom[] atoms = mol.getAtomArray();

        int metals = 0, stereocenters = 0;

        int[] aflags = new int[atoms.length];
        int[] chiral = new int[atoms.length];
        int[][] neighbors = new int[atoms.length][];
        for (int i = 0; i < atoms.length; ++i) {
            /* 
             * make sure no attachments; they can turn off chiral flag
             */
            atoms[i].setAttach(0);
            aflags[i] = atoms[i].getFlags();
            chiral[i] = mol.getChirality(i);

            if (atoms[i].getAtno() == 7
                && (chiral[i] == MolAtom.CHIRALITY_R 
                    || chiral[i] == MolAtom.CHIRALITY_S)) {
                int heavy = 0;
                for (int j = 0; j < atoms[i].getBondCount(); ++j) {
                    MolAtom xa = atoms[i].getBond(j).getOtherAtom(atoms[i]);
                    if (xa.getAtno() != 1)
                        ++heavy;
                }
                
                if (heavy < 4) {
                    logger.warning("** Removing "+getChiralFlag (chiral[i])
                                   +" stereocenter defined "
                                   +"for N atom at "+(i+1));
                    
                    // for N, we make sure it doesn't have a stereocenter
                    //  because of inversion?
                    mol.setChirality(i, 0);
                }
            }

            if (DEBUG) {
                logger.info
                    ("Atom "+(i+1)+": "+atoms[i].getSymbol()+" chiral="
                     +chiral[i] + " flags="+aflags[i]
                     +" parity="+(aflags[i] & MolAtom.PARITY_MASK)
                     +" chirality="+(aflags[i] & MolAtom.CHIRALITY_MASK)
                     +" stereo="+(aflags[i] & MolAtom.ATOMSTEREO_MASK)
                     );
            }

            if (chiral[i] == MolAtom.CHIRALITY_R
                || chiral[i] == MolAtom.CHIRALITY_S) {
                ++stereocenters;
            }

            int[] nb = new int[mol.getNeighborCount(i)];
            for (int j = 0; j < nb.length; ++j) {
                nb[j] = mol.getNeighbor(i, j);
            }
            neighbors[i] = nb;

            if (ElementData.isMetal(atoms[i].getAtno())) {
                ++metals;
            }

            atoms[i].setAtomMap(i+1);
        }

        int dim = mol.getDim();
        if (DEBUG) {
            logger.info("++ Molecule has dimension "+dim);
        }

        DPoint3[] coords = new DPoint3[atoms.length];
        for (int i = 0; i < atoms.length; ++i) {
            coords[i] = atoms[i].getLocation();
        }

        if (DEBUG) {
            logger.info("++ Molecule has "+stereocenters+" stereocenter(s)!");
        }

        if (stereocenters > 0) {
            // count bond parity
            int parity = 0;
            for (MolBond b : mol.getBondArray()) {
                int p = b.getFlags() & MolBond.STEREO1_MASK;
                if (p == MolBond.UP || p == MolBond.DOWN) {
                    ++parity;
                }
            }
            
            if (DEBUG) {
                logger.info("++ Molecule has "+parity+" parity bond(s)!");
            }
            
            if (parity == 0 && dim >= 2) {
                // no bond parity defined, so we generate them
                mol.stereoClean();
            }
        }

        int[][] bflags = new int[atoms.length][atoms.length];
        for (MolBond b : mol.getBondArray()) {
            int ai = mol.indexOf(b.getAtom1());
            int aj = mol.indexOf(b.getAtom2());
            bflags[ai][aj] = b.getFlags();
            bflags[aj][ai] = -1;
        }

        /*disconnectMetals (mol);*/
        handleSpecialMetals (mol);
        
        if (DEBUG) {
            System.err.println("Metals? " + metals 
                               + " " + ChemUtil.canonicalSMILES (mol) 
                               //+ " "+mol.toFormat("smiles:q"));
                               +"\n"+mol.toFormat("mol"));
        }

        // remove bogus radicals that jchem introduces
        fixRadical (mol);
        if (DEBUG) {
            System.err.println("fixRadical: " + ChemUtil.canonicalSMILES (mol)
                               //+" "+mol.toFormat("smiles:q"));
                               +"\n"+mol.toFormat("mol"));
        }

        // fix any bad valence not due to aromaticity
        fixValence (mol);
        if (DEBUG) {
            System.err.println("fixValence: " + ChemUtil.canonicalSMILES (mol) 
                               //+ " " + mol.toFormat("smiles:q"));
                               +"\n"+mol.toFormat("mol"));
        }

        // attempt to fix any bogus aromatic bonds
        kekulize (mol);
        if (DEBUG) {
            System.err.println("kekulize: " + ChemUtil.canonicalSMILES (mol) 
                               //+ " " + mol.toFormat("smiles:q"));
                               +"\n"+mol.toFormat("mol"));
        }

        frags = mol.convertToFrags();

        // vector containing salts
        List<Molecule> salts = new ArrayList<Molecule>();

        SaltIdentifier saltId = SaltIdentifier.getInstance();

        for (int i = 0; i < frags.length; ++i) {
            Molecule f = frags[i];
            f.setName(name);

            if (DEBUG) {
                System.err.println("Normalizing fragment " +i+": " 
                                   + ChemUtil.canonicalSMILES (f) 
                                   //+ " "+f.toFormat("smiles:q"));
                                   +"\n"+f.toFormat("mol"));
            }

            // don't remove anything if there is a metal in the compound
            if (false && metals > 0) {
                // only apply standard normalization rules
                for (int j = 0; j < normalizer.length; ++j) {
                    SMIRKS rule = normalizer[j];
                    if (rule.transform(f)) {
                        if (DEBUG) {
                            System.err.println("++ rule " + j + ": " + rule);
                            System.err.println
                                ("=> " + ChemUtil.canonicalSMILES (f));
                        }
                    }
                }

                if (f.getAtomCount() > 3) {
                    Molecule copy = f.cloneMolecule();
                    try {
                        taugen.generate(f);
                        f = taugen.getCanonicalTautomer();
                    }
                    catch (Exception ex) {
                        logger.log(Level.WARNING, "Can't generate tautomer "
                                   +"for fragment "+i
                                   +"; reason being...", ex);
                        f = copy;
                    }
                }
            }
            else if (f.getAtomCount() > 0) {

                for (int j = 0; j < normalizer.length; ++j) {
                    SMIRKS rule = normalizer[j];

                    boolean valenceError = f.hasValenceError();
                    if (rule.transform(f)) {
                        if (DEBUG) {
                            System.err.println
                                ("++ rule " + j + ": " + rule);
                            System.err.println
                                ("=> " + ChemUtil.canonicalSMILES (f));
                        }
                    }
                    f.valenceCheck();

                    if (f.hasValenceError() && !valenceError) {
                        System.err.println
                            (" ++ " + name + " (frag " + i 
                             + "): bogus structure generated by rule " + j
                             + ": " + rule
                             + " ==> " + ChemUtil.canonicalSMILES (f));
                        System.err.println("** Atom with valence error:");
                        atoms = f.getAtomArray();
                        for (int k = 0; k < atoms.length; ++k) {
                            if (atoms[k].hasValenceError()) {
                                _atom (System.err, k+1, atoms[k]);
                            }
                        }

                        break;
                    }
                    /*
                      else {
                      System.err.println
                      (" ++ " + j + ": " + rules[j] 
                      + " => " + f.toFormat("smiles:a_basu"));
                      }
                    */
                }
                
                // apply external normalization rules
                if (UNMEX) {
                    for (int j = 0; j < externalNormalizer.length; ++j) {
                        SMIRKS rule = externalNormalizer[j];
                        
                        boolean valenceError = f.hasValenceError();
                        if (rule.transform(f)) {
                            if (DEBUG) {
                                System.err.println
                                    ("++ rule " + j + ": " + rule);
                                System.err.println
                                    ("=> " + ChemUtil.canonicalSMILES (f));
                            }
                        }
                        f.valenceCheck();
                        
                        if (f.hasValenceError() && !valenceError) {
                            System.err.println
                                (" ++ " + name + " (frag " + i 
                                 + "): bogus structure generated by rule " + j
                                 + ": " + rule
                                 + " ==> " + ChemUtil.canonicalSMILES (f));
                            System.err.println("** Atom with valence error:");
                            atoms = f.getAtomArray();
                            for (int k = 0; k < atoms.length; ++k) {
                                if (atoms[k].hasValenceError()) {
                                    _atom (System.err, k+1, atoms[k]);
                                }
                            }
                            
                            break;
                        }
                    }
                }

                if (DEBUG) {
                    System.err.println("** Normalized fragment "+i+": " 
                                       + ChemUtil.canonicalSMILES (f) 
                                       + " "+f.toFormat("smiles:q"));
                }

                if (f.getAtomCount() > 3) {
                    Molecule copy = f.cloneMolecule();
                    try {
                        taugen.generate(f);
                        f = taugen.getCanonicalTautomer();
                    }
                    catch (Exception ex) {
                        logger.log(Level.WARNING, "Can't generate tautomer "
                                   +"for fragment "+i
                                   +"; reason being...", ex);
                        f = copy;
                    }

                    if (DEBUG) {
                        int cnt = 0;
                        System.err.println
                            (taugen.getTautomerCount() + " tautomers");
                        for (Enumeration<Molecule> tau = taugen.tautomers();
                             tau.hasMoreElements(); ++cnt) {
                            Molecule t = tau.nextElement();
                            int score = ChemUtil.calcMolScore(t);
                            System.err.println(ChemUtil.canonicalSMILES (t) + " "+cnt
                                               + " " + score);
                        }

                        System.err.println("Canonical tautomer: " 
                                           + ChemUtil.canonicalSMILES (f) 
                                           //+ " " + f.toFormat("smiles:q"));
                                           +"\n"+f.toFormat("mol"));
                    }
                }

                if (removeSalt 
                    && (saltId.isSaltOrSolvent(f) || 
                        (f.getAtomCount() == 1 && 
                         !ElementData.isMetal(f.getAtom(0).getAtno())))) {
                    /*
                      System.err.println("** warning: " + name 
                      + ": fragment " 
                      + f.toFormat("smiles:u") 
                      + " is a salt after transformations!");
                    */
                    if (DEBUG) {
                        System.err.println("## Fragment "+i 
                                           + " is salt/solvent: "
                                           + saltId.identifySaltOrSolvent(f));
                    }
                    salts.add(f);
                    f = null;
                }
            }
            else {
                salts.add(f);
                f = null;
            }

            if (DEBUG) {
                if (f != null) {
                    System.err.println("** Canonicalized fragment "+i+": " 
                                       + ChemUtil.canonicalSMILES (f) 
                                       + " "+f.toFormat("smiles:q"));
                }
                else {
                    System.err.println("** Fragment "+i+" removed");
                }
            }

            frags[i] = f;
        }

        if (salts.size() == frags.length) {
            // all fragments are salt forms... so keep them all
            frags = salts.toArray(new Molecule[0]);
        }

        mol.clear();
        if (frags.length > 1) {
            Arrays.sort(frags, new MolComparator ());
            
            StringBuilder fused = new StringBuilder ();
            List<Molecule> newfrags = new ArrayList<Molecule>();
            for (int i = 0; i < frags.length; ++i) {
                if (frags[i] != null && frags[i].getAtomCount() > 0) {
                    //mol.fuse(frags[i]); // <- doesn't keep stereo
                    if (fused.length() > 0) {
                        fused.append('.');
                    }
                    //fused.append(ChemUtil.canonicalSMILES (frags[i]));
                    // can't use ChemUtil.canonicalSMILES here since we'll loose the
                    //  atom mappings
                    fused.append(frags[i].toFormat("smiles:q"));
                    newfrags.add(frags[i]);
                }
            }
            
            if (DEBUG) {
                System.err.println("Fused: "+fused);
            }
            
            if (fused.length() > 0) {
                if (!MolImporter.importMol(fused.toString(), mol)) {
                    logger.log
                        (Level.SEVERE, "Unable to import molecule: " + fused);
                }
            }

            frags = newfrags.toArray(new Molecule[0]);
            if (frags.length == 1) {
                if (mol instanceof RgMolecule) {
                    ((RgMolecule)mol).setRoot(frags[0]);
                }
                else {
                    frags[0].clonecopy(mol);
                }
            }

            // keep only unique fragments
            /*
            Map<String, Molecule> ufrags = new HashMap<String, Molecule>();
            for (Molecule m : frags) {
                if (m != null && m.getAtomCount() > 0) {
                    for (MolAtom a : f.getAtomArray()) {
                        a.setAtomMap(0);
                    }

                    String hash = hashKeyExt (m, "-");
                    ufrags.put(hash, m);
                }
            }

            frags = ufrags.values().toArray(new Molecule[0]);
            Arrays.sort(frags, new MolComparator ());
            */
        }
        else if (frags.length > 0) {
            if (mol instanceof RgMolecule) {
                ((RgMolecule)mol).setRoot(frags[0]);
            }
            else {
                frags[0].clonecopy(mol);
            }
        }

        if (DEBUG) {
            System.err.println("Postprocessing: " 
                               +ChemUtil.canonicalSMILES (mol)
                               + " "+mol.toFormat("smiles:q")
                               +"\n"+mol.toFormat("mol"));
        }

        mol.setName(name);
        postprocessing (mol);

        if (frags.length > 1) {
            for (int i = 0; i < frags.length; ++i) {
                Molecule f = frags[i];
                if (f != null) {
                    f.setName(name);
                    postprocessing (f, false);
                    if (DEBUG) {
                        logger.info("Component "+(i+1)+": "
                                    +f.toFormat("smiles:q"));
                    }
                }
            }
        }
        else if (frags.length > 0) {
            frags[0] = mol;
        }

        if (DEBUG) {
            logger.info("Number of components: "+frags.length);
            for (Molecule f : frags) {
                if (f != null)
                    System.err.println("  "+f.toFormat("smiles:q"));
            }
        }
        
        /*
         * now fix the atom mapping (if any) for those atoms that
         * don't have a mapping
         */
        fixAtomMapping (neighbors, mol);
        atoms = mol.getAtomArray();

        // restore bond flags
        for (MolBond b : mol.getBondArray()) {
            int i = b.getAtom1().getAtomMap() - 1;
            int j = b.getAtom2().getAtomMap() - 1;
            if (i >= 0 && j >= 0 
                && i < bflags.length && j < bflags[i].length) {
                int f = bflags[i][j];
                if (f < 0) {
                    b.swap();
                    f = bflags[j][i];
                }

                int type = f & MolBond.TYPE_MASK; // original type
                if (type == b.getType()) {
                    b.setFlags(f & MolBond.STEREO_MASK,
                               MolBond.STEREO_MASK);
                    b.setFlags(f & MolBond.TOPOLOGY_MASK, 
                               MolBond.TOPOLOGY_MASK);
                    b.setFlags(f & MolBond.REACTING_CENTER_MASK,
                               MolBond.REACTING_CENTER_MASK);
                }
                else if (b.getType() == 2 
                         && !(b.getAtom1().isTerminalAtom() 
                              || b.getAtom2().isTerminalAtom())
                         && !mol.isRingBond(mol.indexOf(b))) {
                    int s = b.calcStereo2();
                    switch (s) {
                    case MolBond.CIS:
                    case MolBond.TRANS:
                        if (type != 2) {
                            // unspecified due to tautomer
                            b.setFlags(MolBond.CIS|MolBond.TRANS, 
                                       MolBond.CTUMASK);
                            logger.info("Double bond "+(mol.indexOf(b)+1)
                                        +" is unspecified E/Z due "
                                        +"to tautomer!");
                        }
                        else {
                            b.setFlags(s, MolBond.CTUMASK);
                        }
                        break;
                        
                    case MolBond.CTUNSPEC:
                        logger.warning
                            (mol.getName()+": bond "+(mol.indexOf(b)+1)
                             +" has unspecified E/Z");
                    default:
                        b.setFlags(MolBond.CIS|MolBond.TRANS, 
                                   MolBond.CTUMASK);
                    }
                }
            }
        }
        bflags = null;

        if (DEBUG) {
            System.err.println("Restore Bond Flags: " 
                               + ChemUtil.canonicalSMILES (mol)
                               +"\n"+mol.toFormat("mol"));
        }

        // only override stereo if the input molecule has coordinates
        /*
          adjustStereo (dim >= 2, chiral, mol);
        if (DEBUG) {
            System.err.println("Stereo Adjustment: " 
                               + ChemUtil.canonicalSMILES (mol)
                               //+ " "+mol.toFormat("smiles:q"));
                               +"\n"+mol.toFormat("mol"));
        }
        */

        // and finally reset atom map (if any)
        //cleanUp (aflags, mol);
        for (int i = 0; i < atoms.length; ++i) {
            int j = atoms[i].getAtomMap() - 1;
            if (j >= 0 && j < aflags.length) { 
                atoms[i].setXYZ(coords[j].x, coords[j].y, coords[j].z);
                atoms[i].setFlags(aflags[j]);
            }
        }
        // this is important as it recalculate the stereocenters
        mol.setDim(dim);

        if (DEBUG) {
            System.err.println("CleanUp: " + ChemUtil.canonicalSMILES (mol)
                               //+ " "+mol.toFormat("smiles:q"));
                               +"\n"+mol.toFormat("mol"));
        }


        boolean ok = true;
        // if we can't standardize, then revert back???
        if (mol.hasValenceError() && !inputHasValenceError) {
            // the transformations are invalid... revert back
            logger.info("** valence error in canonicalization of " + name);
            System.err.print("  atoms:");
            atoms = mol.getAtomArray();
            for (int i = 0; i < atoms.length; ++i) {
                atoms[i].setAtomMap(i+1);
                if (atoms[i].hasValenceError()) {
                    System.err.print(" " + (i+1)+ "[" 
                                     + atoms[i].getSymbol() + " "
                                     + atoms[i].getCharge() + " "
                                     + atoms[i].getImplicitHcount()+" "
                                     + atoms[i].getValence() + " "
                                     + atoms[i].getRadicalCount()
                                     +"]");
                }
            }
            System.err.println();
            System.err.println("  " + mol.toFormat("smiles:q"));
            for (int i = 0; i < atoms.length; ++i) {
                atoms[i].setAtomMap(0);
            }
            mol.clear();
            origmol.clonecopy(mol);
            frags = origmol.convertToFrags();
            ok = false;
        }

        mol.setName(name);
        for (Map.Entry<String, Object> e : props.entrySet()) {
            mol.setPropertyObject(e.getKey(), e.getValue());
        }

        return ok;
    }

    void adjustStereo (int[] chiral, Molecule mol) {
        adjustStereo (true, chiral, mol);
    }

    static String getChiralFlag (int flag) {
        if (flag == MolAtom.CHIRALITY_R) return "R";
        if (flag == MolAtom.CHIRALITY_S) return "S";
        if (flag == 0) return "R/S";
        return "?";
    }

    static String getParityFlag (int parity) {
        if (parity == MolBond.UP) return "UP";
        if (parity == MolBond.DOWN) return "DOWN";
        return "?";
    }

    void adjustStereo (boolean override, int[] chiral, Molecule mol) {
        if (DEBUG) {
            logger.info("** Stereo adjustments...");
        }

        int[] ringsize = mol.getSmallestRingSizeForIdx();
        Map<MolAtom, List<MolBond>> parity = 
            new HashMap<MolAtom, List<MolBond>>();

        for (MolBond b : mol.getBondArray()) {
            int i = mol.indexOf(b.getAtom1());
            int j = mol.indexOf(b.getAtom2());
            int si = mol.getChirality(i);
            int sj = mol.getChirality(j);
            int mi = b.getAtom1().getAtomMap() - 1;
            int mj = b.getAtom2().getAtomMap() - 1;

            int s = b.getFlags() & MolBond.STEREO1_MASK;
            if (s == MolBond.UP) {                
                if (DEBUG) {
                    System.err.print("Up: ");
                    if (si == MolAtom.CHIRALITY_R) {
                        System.err.print("R");
                    }
                    else if (si == MolAtom.CHIRALITY_S) {
                        System.err.print("S");
                    }
                    else {
                        System.err.print("-");
                    }
                    System.err.print(" ");
                    if (sj == MolAtom.CHIRALITY_R) {
                        System.err.print("R");
                    }
                    else if (sj == MolAtom.CHIRALITY_S) {
                        System.err.print("S");
                    }
                    else {
                        System.err.print("-");
                    }
                    System.err.println();
                }
                if ((si == MolAtom.CHIRALITY_R || si == MolAtom.CHIRALITY_S)
                    && (sj == MolAtom.CHIRALITY_R 
                        || sj == MolAtom.CHIRALITY_S)) {
                    //logger.info("UP: stereo1="+stereo1+" stereo2="+stereo2);
                    /*
                     * a stereocenter must have at least one bond parity such
                     * that point to it
                     */
                    if (b.getAtom1().getAtno() == 7 // point away from N
                        // always point toward the ring
                        || (ringsize[j] > 0 && ringsize[i] == 0) 
                        || (chiral[mi] != si && chiral[mj] != sj)) {
                        // this bond is point in the wrong direction.. flip
                        if (DEBUG) {
                            logger.info("Flipping bond "+(mi+1)+"-"+(mj+1));
                        }
                        b.swap();
                    }
                }
                else if (sj == MolAtom.CHIRALITY_R
                         || sj == MolAtom.CHIRALITY_S) {
                    if (DEBUG) {
                        logger.info("Flipping bond "+(mi+1)+"-"+(mj+1));
                    }
                    b.swap();
                }
                else if (si != MolAtom.CHIRALITY_R 
                         && si != MolAtom.CHIRALITY_S) { // bogus stereo
                    logger.warning("Bond "+(mi+1)+"-"+(mj+1)
                                   +" is UP, but no stereocenter defined!");
                    if (b.getAtom1().getAtno() == 7 // point away from N
                        || b.getAtom1().isTerminalAtom()
                        || (ringsize[j] > 0 && ringsize[i] == 0)
                        ) {
                        b.swap();
                    }
                    else {
                        //b.setFlags(0, MolBond.STEREO1_MASK);
                    }
                }
            }
            else if (s == MolBond.DOWN) {
                if (DEBUG) {
                    System.err.print("Down: ");
                    if (si == MolAtom.CHIRALITY_R) {
                        System.err.print("R");
                    }
                    else if (si == MolAtom.CHIRALITY_S) {
                        System.err.print("S");
                    }
                    else {
                        System.err.print("-");
                    }
                    System.err.print(" ");
                    if (sj == MolAtom.CHIRALITY_R) {
                        System.err.print("R");
                    }
                    else if (sj == MolAtom.CHIRALITY_S) {
                        System.err.print("S");
                    }
                    else {
                        System.err.print("-");
                    }
                    System.err.println();
                }
                
                if ((si == MolAtom.CHIRALITY_R || si == MolAtom.CHIRALITY_S)
                    && (sj == MolAtom.CHIRALITY_R 
                        || sj == MolAtom.CHIRALITY_S)) {
                    //int stereo1 = countBondParity (b.getAtom1(), -1);
                    //int stereo2 = countBondParity (b.getAtom2(), -1);

                    if (b.getAtom1().getAtno() == 7 // N inversion
                        || (ringsize[j] > 0 && ringsize[i] == 0)
                        || (chiral[mi] != si && chiral[mj] != sj)
                        ) {
                        if (DEBUG) {
                            logger.info("Flipping bond "+(mi+1)+"-"+(mj+1));
                        }
                        // this bond is point in the wrong direction.. flip
                        b.swap();
                    }
                }
                else if (sj == MolAtom.CHIRALITY_R
                         || sj == MolAtom.CHIRALITY_S) {
                    if (DEBUG) {
                        logger.info("Flipping bond "+(mi+1)+"-"+(mj+1));
                    }
                    // swap...
                    b.swap();
                }
                else if (si != MolAtom.CHIRALITY_R 
                         && si != MolAtom.CHIRALITY_S) { // bogus stereo
                    logger.warning("Bond "+(mi+1)+"-"+(mj+1)
                                   +" is DOWN, but no stereocenter defined!");
                    // fix it anyway.. 
                    if (b.getAtom1().getAtno() == 7 // point away from N
                        || b.getAtom1().isTerminalAtom()
                        || (ringsize[j] > 0 && ringsize[i] == 0)
                        ) {
                        b.swap();
                    }
                    else { // otherwise, reset it
                        //b.setFlags(0, MolBond.STEREO1_MASK);
                    }
                }
            }

            if (s == MolBond.UP || s == MolBond.DOWN) {
                List<MolBond> bonds = parity.get(b.getAtom2());
                if (bonds == null) {
                    parity.put(b.getAtom2(), bonds = new ArrayList<MolBond>());
                }
                bonds.add(b);
            }
        }

        if (DEBUG) {
            for (MolBond b : mol.getBondArray()) {
                int i = mol.indexOf(b.getAtom1());
                int j = mol.indexOf(b.getAtom2());
                int si = mol.getChirality(i);
                int sj = mol.getChirality(j);
                int s = b.getFlags() & MolBond.STEREO1_MASK;
                if (s == MolBond.UP) {
                    if (DEBUG) {
                        System.err.print("Up: ");
                        if (si == MolAtom.CHIRALITY_R) {
                            System.err.print("R");
                        }
                        else if (si == MolAtom.CHIRALITY_S) {
                            System.err.print("S");
                        }
                        else {
                            System.err.print("-");
                        }
                        System.err.print(" ");
                        if (sj == MolAtom.CHIRALITY_R) {
                            System.err.print("R");
                        }
                        else if (sj == MolAtom.CHIRALITY_S) {
                            System.err.print("S");
                        }
                        else {
                            System.err.print("-");
                        }
                        System.err.println();
                    }
                }
                else if (s == MolBond.DOWN) {
                    if (DEBUG) {
                        System.err.print("Down: ");
                        if (si == MolAtom.CHIRALITY_R) {
                            System.err.print("R");
                        }
                        else if (si == MolAtom.CHIRALITY_S) {
                            System.err.print("S");
                        }
                        else {
                            System.err.print("-");
                        }
                        System.err.print(" ");
                        if (sj == MolAtom.CHIRALITY_R) {
                            System.err.print("R");
                        }
                        else if (sj == MolAtom.CHIRALITY_S) {
                            System.err.print("S");
                        }
                        else {
                            System.err.print("-");
                        }
                        System.err.println();
                    }
                }
            }
        }

        // return if we're not overriding the stereo
        if (!override) 
            return;

        /*
         * Now check verify stereocenters. Here we're only override the 
         * atom's chirality if its original flag is set but for which
         * the calculated chirality based on the coordinates didn't yield
         * anything!!!!
         */
        for (int i = 0; i < mol.getAtomCount(); ++i) {
            MolAtom atom = mol.getAtom(i);
            int map = atom.getAtomMap() - 1;

            List<MolBond> bonds = parity.get(atom);
            if (bonds != null) {
                StringBuilder mesg = new StringBuilder ();
                int up = 0, down = 0;
                for (MolBond b : bonds) {
                    int s = b.getFlags() & MolBond.STEREO1_MASK;
                    if (s == MolBond.UP) ++up;
                    else ++down;
                    mesg.append(" "+(s == MolBond.UP ? "UP" :"DOWN"));
                }
                //logger.info("Atom "+(map+1)+" parity:"+mesg);

                if (up > 1 || down > 1) {
                    logger.warning("Bogus parity at atom "+(map+1)
                                   +"; arbitrary removing one...");
                    MolBond b = bonds.iterator().next();
                    //b.setFlags(0, MolBond.STEREO1_MASK);
                    b.swap();
                }
            }

            int chirality =  mol.getChirality(i);
            if (chirality == MolAtom.CHIRALITY_R 
                || chirality == MolAtom.CHIRALITY_S) {
                // keep this
                if (chirality != chiral[map]) {
                    logger.warning("Stereocenter at atom "
                                   +atom.getSymbol()+":"+(map+1)+" changed "
                                   +"from "+getChiralFlag (chiral[map])
                                   +" to "+getChiralFlag (chirality)+"!");
                }
            }
            else if (chiral[map] == MolAtom.CHIRALITY_R
                     || chiral[map] == MolAtom.CHIRALITY_S) {
                logger.warning("No stereocenter defined at atom "
                               +atom.getSymbol()+":"+(map+1)
                               +" but was original defined as "
                               +getChiralFlag (chiral[map]));
                try {
                    mol.setChirality(i, chiral[map]);
                    int count = countBondParity (atom, 0);
                    logger.warning("Overriding chirality for atom "
                                   +atom.getSymbol()+":"+(map+1)
                                   +" ("+count+"/"+atom.getBondCount()
                                   +" parity) "+ getChiralFlag (chiral[map])
                                   +" "+getChiralFlag (mol.getChirality(i)));

                    /*
                    for (int k = 0; k < atom.getBondCount(); ++k) {
                        MolBond b = atom.getBond(k);
                        MolAtom xa = b.getOtherAtom(atom);
                        int s = b.getFlags() & MolBond.STEREO1_MASK;
                        String dir = "-";
                        if (s == MolBond.UP) dir = "UP";
                        else if (s == MolBond.DOWN) dir = "DOWN";
                        logger.info(dir+" neighbor "+xa.getSymbol()+":"
                                    +(xa.getAtomMap()+1));
                    }
                    */
                }
                catch (Exception ex) {
                    // bogus parity, so clear out the bond
                    logger.warning("Can't override chiral atom "
                                   +(map+1)+"; "+ex);
                }
            }
        }

        // not sure what exact this does
        mol.adjustMultiChiralFlag();
    }

    void cleanUp (int[] aflags, Molecule mol) {
        MolAtom[] atoms = mol.getAtomArray();

        int[] grinv = new int[atoms.length];
        mol.getGrinv(grinv);

        for (int j = 0; j < atoms.length; ++j) {
            MolAtom a = atoms[j];

            int i = a.getAtomMap() - 1;
            if (i >= 0 && i < aflags.length) { 
                a.setFlags(aflags[i]);
            }

            /*
             * now we fix stereocenters for which there is no corresponding
             * wedge bond pointing to them
             */
            int chiral = mol.getChirality(j);
            if (chiral != 0 && a.getAtno() != 7) {
                int stereo1 = countBondParity (a, -1);
                int stereo2 = countBondParity (a, 1);
                if (DEBUG) {
                    logger.info("** Atom "+(j+1)+" stereo2="+stereo2
                                +" stereo1="+stereo1);
                }

                if (stereo2 == 1 && stereo1 == 0) {
                    MolBond ref = null;

                    List<MolBond> cands = new ArrayList<MolBond>();
                    for (int n = 0; n < a.getBondCount(); ++n) {
                        MolBond b = a.getBond(n);
                        int p = b.getFlags() & MolBond.STEREO1_MASK;
                        if (p == MolBond.UP || p == MolBond.DOWN) {
                            ref = b;
                        }
                        else if (p == 0) {
                            cands.add(b);
                        }
                    }

                    if (ref == null) { 
                        // this shouldn't happen!
                        logger.warning("stereo2 is "+stereo2
                                       +" but no defined bond parity found!");
                    }
                    else if (cands.isEmpty()) {
                        logger.warning("Stereocenter at atom "
                                       +(j+1)+" is likely bogus!");
                    }
                    else {
                        MolBond selection = null;
                        if (cands.size() == 1) {
                            selection = cands.get(0);
                        }
                        else {
                            // first prefer bond with no attached 
                            //  stereocenter
                            int minGi = atoms.length+1;
                            for (MolBond b : cands) {
                                MolAtom xa = b.getOtherAtom(a);
                                int xi = mol.indexOf(xa);
                                int c = mol.getChirality(xi);
                                // in case there are multiple options, we 
                                //   select one that has the lowest graph
                                //   invariant
                                if (c != MolAtom.CHIRALITY_R
                                    && c != MolAtom.CHIRALITY_S) {

                                    if (grinv[xi] < minGi) {
                                        minGi = grinv[xi];
                                        selection = b;
                                    }
                                }
                            }

                            if (selection == null) {
                                // saturated... just pick one with lowest gi
                                minGi = atoms.length + 1;
                                for (MolBond b : cands) {
                                    MolAtom xa = b.getOtherAtom(a);
                                    int xi = mol.indexOf(xa);
                                    if (grinv[xi] < minGi) {
                                        minGi = grinv[xi];
                                        selection = b;
                                    }
                                }
                            }
                        }

                        if (selection.getAtom1() != a) {
                            // need to flip this bond
                            selection.swap();
                        }

                        // first try down
                        selection.setFlags(MolBond.DOWN, MolBond.STEREO1_MASK);
                        if (chiral == mol.getChirality(j)) {
                            // do nothing
                        }
                        else { // now try up
                            selection.setFlags
                                (MolBond.UP, MolBond.STEREO1_MASK);
                        }

                        if (chiral != mol.getChirality(j)) {
                            logger.warning("Can't figure out what parity "
                                           +"for bond "
                                           +(mol.indexOf(selection)+1)
                                           +"; revert back!");
                            selection.setFlags(0, MolBond.STEREO1_MASK);
                        }
                        else {
                            int p = selection.getFlags() 
                                & MolBond.STEREO1_MASK;
                            logger.info("Using bond "
                                        +(mol.indexOf(selection)+1)
                                        +" with parity "
                                        +(p==MolBond.UP?"DOWN":"UP")
                                        +" to augment stereocenter at atom "
                                        +(j+1));
                        }
                    }
                }
                else if (stereo1 > 1 && stereo2 == 0) { 
                    // two wedge bonds pointing toward the stereocenter
                    // if one of the bonds has stereocenters on both of
                    // its ends, then we flip it
                    for (int n = 0; n < a.getBondCount(); ++n) {
                        MolBond b = a.getBond(n);
                        int parity = b.getFlags() & MolBond.STEREO1_MASK;
                        if (parity == MolBond.UP || parity == MolBond.DOWN) {
                            MolAtom xa = b.getOtherAtom(a);
                            int xi = mol.indexOf(xa);
                            int c = mol.getChirality(xi);
                            if (c == MolAtom.CHIRALITY_R
                                || c == MolAtom.CHIRALITY_S) {
                                if (DEBUG) {
                                    logger.info("Flipping bond "+(j+1) +" and "
                                                +(xi+1));
                                }
                                b.swap();
                            }
                        }
                    }
                }
            }
        }
    } 

    /**
     * here are the rules for fixing bond parity.  the input bond is 
     * assumed to have a parity of either UP or DOWN and that both of 
     * its ends are stereocenters.
     */
    static void fixBondParity (MolBond bond) {
        Molecule mol = (Molecule)bond.getParent();

        int[] gi = new int[mol.getAtomCount()];
        mol.getGrinv(gi);

        int s = bond.getFlags() & MolBond.STEREO1_MASK;

        // find a neighboring bond with the lowest graph invariant
        //   atom and 
        List<Integer> cands = new ArrayList<Integer> ();
        for (int n = 0; n < bond.getAtom2().getBondCount(); ++n) {
            MolBond nb = bond.getAtom2().getBond(n);
            int f = nb.getFlags() & MolBond.STEREO1_MASK;
            if (nb.getType() == 1 && f == 0) {
                int xi = mol.indexOf(nb.getOtherAtom(bond.getAtom2()));
                int chiral = mol.getChirality(xi);
                if (chiral != MolAtom.CHIRALITY_R
                    && chiral != MolAtom.CHIRALITY_S) {
                    cands.add(xi);
                }
            }
        }
        
        if (cands.isEmpty()) {
            // all neighboring atoms are stereocenters
        }
        else {
            // find the one with the lowest graph invariant
            int nb = cands.get(0);
            int min = gi[nb];
            for (int k = 1; k < cands.size(); ++k) {
                int n = cands.get(k);
                if (gi[n] < min) {
                    min = gi[n];
                    nb = n;
                }
            }
            int bi = mol.getBtab()[nb][mol.indexOf(bond.getAtom2())];
            logger.warning("Moving bond parity "+s+" from "
                           +(mol.indexOf(bond)+1)+" to "+(bi+1));
            MolBond cb = mol.getBond(bi);
            // mv the parity mask to this bond
            cb.setFlags(s, MolBond.STEREO1_MASK);
            bond.setFlags(0, MolBond.STEREO1_MASK);
        }
    }

    /*
     * count number of bond parity that points in the 
     *  designated direction:
     *  0 - either
     *  < 0 - atom1
     *  > 0 - atom2
     */
    static int countBondParity (MolAtom atom, int dir) {
        int count = 0;

        Molecule mol = (Molecule)atom.getParent();
        if (DEBUG) {
            logger.info("** Atom "+(mol.indexOf(atom)+1));
        }
        for (int n = 0; n < atom.getBondCount(); ++n) {
            MolBond nb = atom.getBond(n);
            int parity = nb.getFlags() & MolBond.STEREO1_MASK;
            if (parity == MolBond.UP) {
                if (DEBUG) {
                    System.err.println(" UP: "+(mol.indexOf(nb.getAtom1())+1)
                                       +"-"+(mol.indexOf(nb.getAtom2())+1));
                }
                if (dir < 0) {
                    if (nb.getAtom1() == atom) {
                        ++count;
                    }
                }
                else if (dir > 0) {
                    if (nb.getAtom2() == atom) {
                        ++count;
                    }
                }
                else {
                    ++count;
                }
            }
            else if (parity == MolBond.DOWN) {
                if (DEBUG) {
                    System.err.println(" DOWN: "+(mol.indexOf(nb.getAtom1())+1)
                                       +"-"+(mol.indexOf(nb.getAtom2())+1));
                }
                if (dir < 0) {
                    if (nb.getAtom1() == atom) {
                        ++count;
                    }
                }
                else if (dir > 0) {
                    if (nb.getAtom2() == atom) {
                        ++count;
                    }
                }
                else {
                    ++count;
                }
            }
        }
        return count;
    }


    static int getNeighborStereocenterCount (MolAtom atom) {
        Molecule mol = (Molecule)atom.getParent();
        int stereo = 0; // number of neighboring stereocenters

        for (int n = 0; n < atom.getBondCount(); ++n) {
            MolBond nb = atom.getBond(n);

            int xi = mol.indexOf(nb.getOtherAtom(atom));
            int chiral = mol.getChirality(xi);
            if (chiral == MolAtom.CHIRALITY_R
                || chiral == MolAtom.CHIRALITY_S) {
                ++stereo;
            }
        }
        
        return stereo;
    }

    void fixAtomMapping (int[][] neighbors, Molecule mol) {
        /*
         * there is a bug in chemaxon that doesn't properly
         * parse the atom mapping when the atom is chiral
         * e.g., it doesn't pick up the mapping 8 in [N@:8]
         */
        BitSet mapped = new BitSet (neighbors.length+1);
        List<MolAtom> unmapped = new ArrayList<MolAtom>();
        for (int i = 0; i < mol.getAtomCount(); ++i) {
            MolAtom a = mol.getAtom(i);
            int m = a.getAtomMap();
            if (m == 0) {
                if (DEBUG) {
                    logger.warning("Atom "+(i+1)+" has no mapping!");
                }
                unmapped.add(a);
            }
            else {
                mapped.set(m-1);
            }
        }
                    
        for (MolAtom a : unmapped) {
            int map = -1;
            for (int i = 0; i < a.getBondCount(); ++i) {
                MolBond b = a.getBond(i);
                MolAtom xa = b.getOtherAtom(a);
                int m = xa.getAtomMap();
                if (m == 0) {
                    logger.log(Level.SEVERE, 
                               "Multiple unmapped atoms!");
                }
                else {
                    int[] nb = neighbors[m-1];
                    //System.err.print(m+":");
                    for (int j = 0; j < nb.length; ++j) {
                        //System.err.print(" "+(nb[j]+1));
                        if (!mapped.get(nb[j])) {
                            //System.err.print("*");
                            // this must correspond to the unmapped
                            //  atom?
                            if (map < 0) map = nb[j];
                            else if (map != nb[j]) {
                                logger.log(Level.SEVERE, "Atoms "+(map+1)
                                           +" and "+(nb[j]+1)
                                           +" are unmapped!");
                            }
                        }
                    }
                    //System.err.println();
                }
            }
            if (map < 0) {
                logger.log(Level.SEVERE, "Can't unambigously assign "
                           +"mapping to atom " +(mol.indexOf(a) + 1));
            }
            else {
                a.setAtomMap(map+1);
            }
        }
    }

    int chargeBalance (Molecule mol) {
        int pos = 0, neg = 0;
        for (MolAtom a : mol.getAtomArray()) {
            int ch = a.getCharge();
            if (ch < 0) ++neg;
            else if (ch > 0) ++pos;
        }

        if ((pos > 0 || neg > 0) && (mol.getAtomCount() > 501)) {
            logger.warning("Molecule contains "+mol.getAtomCount()
                           +" atom(s); with +"+pos+" -"+neg+" no attempt to "
                           +"(de-) protonate!");
            return 0;
        }
        
        if (pos > neg) {
            deprotonate (mol);
            if (DEBUG) {
                System.err.println("Deprotonate: " 
                                   + ChemUtil.canonicalSMILES(mol) 
                                   + " " + mol.toFormat("smiles:q"));
            }

            int net = 0;
            for (MolAtom a : mol.getAtomArray()) {
                net += a.getCharge();
            }

            if (net > 0) {
                // now see if we can add electrons to make this
                // a zwitterion
                MolAtom[] atoms = mol.getAtomArray();
                int[] rank = new int[atoms.length];
                mol.getGrinv(rank);

                List<Integer> cand = new ArrayList<Integer>();
                for (int i = 0; i < atoms.length; ++i) {
                    int ch = atoms[i].getCharge();
                    int hc = atoms[i].getImplicitHcount();
                    switch (atoms[i].getAtno()) {
                    case 8: // O
                    case 17: // Cl
                    case 35: // Br
                    case 53: // I
                        if (hc > 0 && ch == 0) {
                            cand.add(i);
                        }
                        break;
                    }
                }

                while (net > 0 && !cand.isEmpty()) {
                    int minrank = Integer.MAX_VALUE;
                    Integer atom = null;
                    for (Integer a : cand) {
                        if (rank[a] < minrank
                            || (rank[a] == minrank
                                && atoms[a].getAtno()
                                < atoms[atom].getAtno())) {
                            minrank = rank[a];
                            atom = a;
                        }
                    }

                    if (atom != null) {
                        cand.remove(atom);
                        atoms[atom].setCharge(-1);
                        --net;
                    }
                }

                neg = 0; // leave the negative charges unchanged
            }
        }
        
        if (neg > pos) {
            protonate (mol);
            if (DEBUG) {
                System.err.println("Standardized molecule has -"+neg+" charge "
                                   +"; attempt to neutralize...");
                System.err.println(ChemUtil.canonicalSMILES (mol));
            }
        }

        int net = 0;
        for (MolAtom a : mol.getAtomArray())
            net += a.getCharge();
        
        return net;
    }

    void postprocessing (Molecule mol) {
        postprocessing (mol, true);
    }
    
    void postprocessing (Molecule mol, boolean chargeBalance) {
        if (chargeBalance) chargeBalance (mol);
        if (mol.getAtomCount() > 0) {
            /*
             * this is needed because of tests/standardizer_case10.smi
             */
            try {
                String smiles = mol.toFormat("cxsmiles:u");
                if (DEBUG) {
                    logger.log(Level.INFO, "Before final molecule cleaning: "
                               +smiles);
                }

                Molecule newmol = MolImporter.importMol(smiles);
                newmol.dearomatize();

                boolean hasArom = false;
                for (MolBond b : newmol.getBondArray()) {
                    if (b.getType() == MolBond.AROMATIC) {
                        hasArom = true;
                        break;
                    }
                }

                if (hasArom) {
                    logger.log(Level.WARNING, mol.getName()
                               +": can't kekulize aromatized form: "+smiles);
                    // leave mol as-is
                }
                else {
                    if (DEBUG) {
                        logger.info("Final molecule clean-up\n"
                                    +newmol.toFormat("mol"));
                    }
                    newmol.clonecopy(mol);
                }
            }
            catch (Exception ex) {
                logger.log(Level.SEVERE, 
                           "Can't cleanup molecule properly", ex);
            }
        }
    }

    public Molecule[] getFragments () { return frags; }
    public int getFragmentCount () { 
        return frags != null ? frags.length : 0; 
    }

    public Molecule getLargestFragment () { return frags[0]; }



    public static String hashKey (Molecule mol) {
        return hashKey (mol, "-");
    }

    public static String hashKey (Molecule mol, String sep) {
        String[] keys = hashKeyArray (mol);
        return keys[0]+sep+keys[1]+sep+keys[2]+sep+keys[3];
    }

    
    private static Molecule getLayer3Equivalent(Molecule m){
    	 Molecule m0=m.cloneMolecule();
    	 
         int[] atno = new int[m0.getAtomCount()];
         for (int i = 0; i < atno.length; ++i) {
             MolAtom a = m0.getAtom(i);
             a.setRadical(0);
             a.setCharge(0);
             a.setFlags(0);
             a.setMassno(0);
             a.setAtomMap(i+1);
         }
         for (MolBond b : m0.getBondArray()) {
            b.setStereo2Flags(b.getNode1(), b.getNode2(), 0);
            if(b.isQuery()){ //hack
            	b.setFlags(1);
            }
         }
         Molecule mout = new Molecule();
         ChemUtil.canonicalSMILES(mout,m0,false);
         
         return mout;
    }
    
    private static Molecule getLayer2Equivalent(Molecule m){
   	 Molecule m0=m.cloneMolecule();
   	 
        int[] atno = new int[m0.getAtomCount()];
        for (int i = 0; i < atno.length; ++i) {
            MolAtom a = m0.getAtom(i);
            a.setRadical(0);
            a.setCharge(0);
            a.setFlags(0);
            a.setMassno(0);
            a.setAtomMap(i+1);
        }
        for (MolBond b : m0.getBondArray()) {
           b.setStereo2Flags(b.getNode1(), b.getNode2(), 0);
           if(b.isQuery()){ //hack
           	b.setFlags(1);
           }
           b.setType(1);    //force single
        }
        Molecule mout = new Molecule();
        ChemUtil.canonicalSMILES(mout,m0,false);
        
        return mout;
   }
    
    private static Molecule getLayer1Equivalent(Molecule m){
      	 Molecule m0=m.cloneMolecule();
      	 
           int[] atno = new int[m0.getAtomCount()];
           for (int i = 0; i < atno.length; ++i) {
               MolAtom a = m0.getAtom(i);
               a.setRadical(0);
               a.setCharge(0);
               a.setFlags(0);
               a.setMassno(0);
               a.setAtno(6);
               a.setAtomMap(i+1);
           }
           for (MolBond b : m0.getBondArray()) {
              b.setStereo2Flags(b.getNode1(), b.getNode2(), 0);
              if(b.isQuery()){ //hack
              	b.setFlags(1);
              }
              b.setType(1);    //force single
           }
           Molecule mout = new Molecule();
           ChemUtil.canonicalSMILES(mout,m0,false);
           
           return mout;
      }
    
    

    /**
     * Takes an int array and replaces each value with the position it would have in a sorted
     * unique list.
     * 
     * For example:
     * [4,2,1,10,42234,1,42234]  
     * 
     * would become:
     * 
     * [2,1,0,3,4,0,4]
     * 
     * 
     * @param rank
     */
    public static void normalizeRanks(int[] rank){

        int[] last=new int[]{-1, -1};
        IntStream.range(0, rank.length)
        	.mapToObj(i->new int[]{i,rank[i]})
        	.sorted((a,b)->b[1]-a[1])
        	.forEach(ab->{
        		if(last[1]!=ab[1]){
        			last[0]++;
        			last[1]=ab[1];
        		}
        		rank[ab[0]]=last[0];
        	});
    }
    
    /**
     * Simple morgan's algorithm for graph invariants. This requires k*N operations
     * where k is a constant that is large enough to "absorb" the whole graph (13 here).
     * 
     * @param m
     * @return
     */
    public static long[] morgans(Molecule m){
    	int MAX_ROUND = 13;
    	int[] atno = new int[m.getAtomCount()];
        for (int i = 0; i < atno.length; ++i) {
            MolAtom a = m.getAtom(i);
            atno[i]=a.getAtno();
        }
        long[] rank;
        {
            
            long[][] hash = new long[MAX_ROUND][atno.length];
            for (int i = 0; i < atno.length; ++i)
                hash[0][i] = atno[i];

            int round = 1;
            for (; round < MAX_ROUND; ++round) {
                int p = round - 1;
                for (int i = 0; i < atno.length; ++i) {
                    MolAtom a = m.getAtom(i);
                    long ha = hash[p][i];
                    for (int j = 0; j < a.getBondCount(); ++j) {
                        MolAtom xa = a.getBond(j).getOtherAtom(a);
                        int k = m.indexOf(xa);
                        ha +=  ha += (a.getBond(j).getType() << xa.getImplicitHcount())
                                + hash[p][k];
                    }
                    if (ha < 0) {
                        if (DEBUG) {
                            logger.log(Level.SEVERE,
                                       "OVERFLOW AT ITERATION "+round+"!");
                        }
                        ha = hash[round-1][i];
                    }
                    hash[round][i] = ha;
                }
            }
            rank = hash[round-1];
        }
        return rank;
    }
    
    /**
     * Return the morgan's algorithm values as a sorted list, encoded
     * as a string. 
     * @param m
     * @return
     */
    public static String morgansAsString(Molecule m){
    	long[] order = Arrays.stream(morgans(m))
    			             .sorted()
    			             .toArray();
    	
    	String s=new BigInteger(BitSet.valueOf(order).toByteArray()).toString(64);
    	return s;
    }
    
    /**
     * Extended version of the hash key that includes the topology+label
     *  layer that sits between the first and second layers of previous
     *  key.
     */
    public static String[] hashKeyArray (Molecule input) {
        if (input == null || input.getAtomCount() == 0) {
            return hashChain45 ("", "", "", "");
        }

        Molecule mol = null;
        String molstr = ChemUtil.canonicalSMILES(input);
        try {
            MolHandler mh = new MolHandler (molstr);
            mol = mh.getMolecule();
        }
        catch (Exception ex) {
            logger.log(Level.SEVERE, 
                       "Can't parse canonical SMILES: "+molstr, ex);
            mol = input.cloneMolecule();
        }

        Molecule m0 = mol.cloneMolecule();
        m0.expandSgroups();
        m0.hydrogenize(false);
        // make sure H isotopes are suppressed
        m0.implicitizeHydrogens(MolAtom.ALL_H);

        // remove any leftover H that didn't get removed 
        // above
        for (MolAtom a : m0.getAtomArray()) {
            if (a.getAtno() == 1)
                m0.removeNode(a);
        }
        
        String molstr1=morgansAsString(getLayer1Equivalent(m0));
        String molstr2=morgansAsString(getLayer2Equivalent(m0));
        String molstr3=morgansAsString(getLayer3Equivalent(m0));
        
        
        return hashChain45 (molstr1, molstr2, molstr3, molstr);
    }

    static String[] hashChain45 (String... strs) {
        String[] hkeys = new String[strs.length];
        try {
            MessageDigest md = MessageDigest.getInstance("SHA1");
            char[] chksum = new char[strs.length-1];
            byte[] ph = {};
            for (int i = 0; i < strs.length; ++i) {
                md.reset();
                md.update(ph); // hash chaining
                byte[] digest = md.digest(strs[i].getBytes("UTF-8"));
                String h = Base32.encode45(digest);
                if (i < chksum.length)
                    chksum[i] = h.charAt(h.length()-1);
                StringBuilder sb = new StringBuilder ();
                for (int j = 0; j < i; ++j)
                    sb.append(chksum[j]);
                sb.append(h);
                hkeys[i] = sb.toString();
                ph = digest;
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return hkeys;
    }


    static void process (String tag, java.io.InputStream is, PrintStream os) 
        throws Exception {
        /*
         * java -Dketo-enol=true ... 
         * to turn on keto-enol tautomerism
         */
        TautomerGenerator tg = 
            // new NCGCTautomerGenerator ()
            new SayleDelanyTautomerGenerator (1001);

        if (Boolean.getBoolean("lychi-keto-enol")) {
            // hanlding keto-enol... might be too slow
            ((SayleDelanyTautomerGenerator)tg).set
                (SayleDelanyTautomerGenerator.FLAG_ALL);
            logger.info("## Keto-Enol tautomerism is on");
        }
        else {
            logger.info("## Keto-Enol tautomerism is off");
        }

        /*
         * -Dkeep-salt={true|false}
         * to turn off/on strip salt/solvent
         */
        boolean removeSaltSolvent = !Boolean.getBoolean("lychi-keep-salt");
        logger.info("## Salt/solvent stripping is "
                    +(removeSaltSolvent ? "on":"off"));

        /*
         * -Doutput-format={smi,smiles,sdf}
         * output format
         */
        String format = System.getProperty("lychi-format", "smi");
        logger.info("## Output format: "+format);

        LyChIStandardizer msz = new LyChIStandardizer (tg);
        msz.removeSaltOrSolvent(removeSaltSolvent);

        MolImporter mi = new MolImporter (is);
        for (Molecule mol = new Molecule (); ; ) {
            try {
                if (!mi.read(mol))
                    break;

                String name = mol.getName();
                if (name != null && name.length() > 0) {
                }
                else if (tag == null) {
                    // assume we have smiles input here
                    for (int i = 0; i < 10 
                             && (name == null || name.equals("")); ++i) {
                        name = mol.getProperty("field_"+i);
                    }
                    mol.setName(name);
                }
                else {
                    name = mol.getProperty(tag);
                    mol.setName(name);
                }
                
                mol.valenceCheck();
                if (mol.hasValenceError()) {
                    System.err.println("** warning: " + mol.getName() 
                                       + " has valence error!");
                }
                
                int dim = mol.getDim();
                msz.standardize(mol);
                String smi = ChemUtil.canonicalSMILES (mol, true);
                if (!format.startsWith("smi")) {
                    mol.setProperty("LyChI_SMILES", smi);
                    String hk = hashKey (mol);
                    mol.setProperty("LyChI_HK", hk);

                    String []toks = hk.split("-");
                    mol.setProperty("LyChI_HK1", toks[0]);
                    mol.setProperty("LyChI_HK2", toks[1]);
                    mol.setProperty("LyChI_HK3", toks[2]);
                    mol.setProperty("LyChI_HK4", toks[3]);

                    os.print(mol.toFormat("sdf"));
                }
                else {
                    //os.println("## fragments: "+msz.getFragmentCount());
                    os.println(smi
                               +"\t"+name+"\t"+hashKey (mol));
                    /*
                    if (!removeSaltSolvent && msz.getFragmentCount() > 1) {
                        // if we're not removing salt then output all components
                        for (Molecule frag : msz.getFragments()) {
                            if (frag != null) {
                                os.println(ChemUtil.canonicalSMILES(frag) 
                                           + "\t" + name 
                                           + "\t" + hashKey (frag)); 
                            }
                        }
                    }
                    */
                }
            }
            catch (Exception ex) {
                logger.log(Level.SEVERE, "** can't process molecule "
                           +mol.getName()+": "+ex.getMessage(), ex);
            }
            catch (StackOverflowError err) {
                logger.log(Level.SEVERE,"** can't parse molecule at line: "
                           +mi.getLineCount());
            }
        }
        mi.close();
    }
    
    public static void main(String[] argv) throws Exception {
    
        String tag = null;
        
        List<String> input = new ArrayList<String>();
        for (int i = 0; i < argv.length; ++i) {
            if (argv[i].startsWith("tag=")) {
                tag = argv[i].substring(4);
            }
            else {
                input.add(argv[i]);
            }
        }

        if (input.isEmpty()) {
            logger.info
                ("** Reading input from STDIN & writing to STDOUT");
            process (tag, System.in, System.out);
        }
        else {
            PrintStream os = System.out;
            String out = System.getProperty("output");
            if (out != null) {
                os = new PrintStream (new FileOutputStream (out));
                logger.info("** Writing output to \""+out+"\"...");
            }

            for (String file : input ) {
                process (tag, new java.io.FileInputStream (file), os);
            }

            if (os != System.out)
                os.close();
        }
    }
}
