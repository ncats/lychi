package lychi;

import java.util.BitSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.concurrent.locks.ReentrantLock;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.formats.MolImporter;
import chemaxon.formats.MolFormatException;
import chemaxon.sss.search.SearchException;
import chemaxon.sss.search.MolSearch;
import static chemaxon.sss.SearchConstants.*;
import chemaxon.util.MolHandler;


/*
 * poor's man version of chemaxon.reaction.Reactor; this class is
 * needed to get around the need to have a license.. sigh!
 */
public class SMIRKS {
    static private boolean debug = false;
    static {
        try {
            debug = Boolean.getBoolean("smirks.debug");
        }
        catch (Exception ex) {
        }
    }
    
    private static final Logger logger = 
        Logger.getLogger(SMIRKS.class.getName());
    
    static final int FP_SIZE = 16;
    static final int FP_LENGTH = 6;
    static final int FP_DEPTH = 2;

    private Molecule source;
    private Molecule target;
    private MolSearch ms;
    private int[] map;
    private String smirks;
    private int[] sourcefp;

    private final ReentrantLock lock = new ReentrantLock ();
    
    public SMIRKS (String smirks) throws MolFormatException {
        int pos = smirks.indexOf(">>");
        if (pos < 0) {
            throw new IllegalArgumentException 
                ("Invalid SMIRKS format: " + smirks 
                 + "; no reactant and/or product specified!");
        }
        setTransform (smirks.substring(0, pos), 
                      smirks.substring(pos+2));
    }

    public SMIRKS (String reactant, String product) 
        throws MolFormatException {
        setTransform (reactant, product);
    }

    public Molecule getReactant () { return source; }
    public Molecule getProduct () { return target; }
    public String getTransform () { return smirks; }

    public void setTransform (String source, String target) 
        throws MolFormatException {
        MolHandler mh = new MolHandler ();
        mh.setQueryMode(true);
        mh.setMolecule(source);
        sourcefp = mh.generateFingerprintInInts
            (FP_SIZE, FP_DEPTH, FP_LENGTH);

        this.source = mh.getMolecule();
        if (target == null || target.length() == 0) {
            this.target = new Molecule ();
        }
        else {
            mh.setMolecule(target);
            //mh.aromatize();
            this.target = mh.getMolecule();
        }

        ms = createMolSearch (this.source);
        this.smirks = source + ">>" + target;

        /*
          if (this.source.getAtomCount() > this.target.getAtomCount()) {
          throw new IllegalArgumentException 
          ("Invalid SMIRKS transform specified; source pattern "
          + " contains more atoms than target!");
          }
        */

        MolAtom[] satoms = this.source.getAtomArray();
        MolAtom[] tatoms = this.target.getAtomArray();
        map = new int[satoms.length];
        BitSet ignore = new BitSet (tatoms.length);
        for (int i = 0; i < satoms.length; ++i) {
            map[i] = -1;
            int m = satoms[i].getAtomMap();
            if (m > 0) { // there is atom mapping...
                // find out which of the target atom has the same atom
                //  mapping...
                for (int j = 0; j < tatoms.length; ++j) {
                    int n = tatoms[j].getAtomMap();
                    if (n == m) {
                        if (map[i] != -1) {
                            // an atom from source is mapped to multiple
                            //  target atoms... 
                            throw new IllegalArgumentException
                                ("Atom " + (i+1) + " in " + source + " "
                                 + "is mapped to multiple atoms in "
                                 + target);
                        }
                        map[i] = j;
                    }
                    else if (n > 0) { 
                        // this atom mapping isn't mapped to anything
                        if (!ignore.get(j)) {
                            /*
                              System.err.println
                              ("** warning: atom label " + n + " in "
                              + "product " + target + " doesn't have "
                              + "any mapping to the corresponding "
                              + "reactant!");
                            */
                            ignore.set(j);
                        }
                    }
                }
            }
        }

        if (tatoms.length == 0) {
            // delete all atoms & bonds
            return;
        }

        /*
          System.out.println(smirks);
          for (int i = 0, j = 0; i < map.length; ++i) {
          if (map[i] < 0) {
          }
          else {
          System.out.println((i+1) + " <=> " + (map[i]+1));
          }
          }
        */

        // now go through the unmatched atoms and map them sequentially
        for (int i = 0, j = 0; i < map.length; ++i) {
            if (map[i] < 0) {
                // skip all mapped atoms
                while (tatoms[j].getAtomMap() > 0 || ignore.get(j)) {
                    ++j;
                }
                ignore.set(j);
                map[i] = j++;
            }
        }

        if (debug) {
            System.out.println("mapping for " + source + ">>" + target);
            for (int i = 0; i < map.length; ++i) {
                int m = satoms[i].getAtomMap();
                int n = tatoms[map[i]].getAtomMap();
                // assert (m == n);
                System.out.println((i+1) +": "
                                   + p(satoms[i]) + " >> " 
                                   + (map[i]+1) + p(tatoms[map[i]]));
            }
        }
    }

    public String toString () {
        StringBuffer sb = new StringBuffer (smirks);
        for (int i = 0; i < map.length; ++i) {
            if (map[i] < 0) {
            }
            else {
                int m = source.getAtom(i).getAtomMap();
                int n = target.getAtom(map[i]).getAtomMap();
                sb.append("\n");
                sb.append((i+1) +": "
                          + p(source.getAtom(i)) + " >> " 
                          + (map[i]+1) + p(target.getAtom(map[i])));
            }
        }
        return sb.toString();
    }

    static MolSearch createMolSearch (Molecule query) {
        MolSearch ms = new MolSearch ();
        ms.setOption(OPTION_VALENCE_MATCHING, VALENCE_MATCHING_ON);
        ms.setOption(OPTION_CHARGE_MATCHING, CHARGE_MATCHING_EXACT);
        //ms.setOption(OPTION_ISOTOPE_MATCHING, ISOTOPE_MATCHING_IGNORE);
        ms.setOption(OPTION_ISOTOPE_MATCHING, ISOTOPE_MATCHING_EXACT);
        ms.setExactStereoMatching(true);
        //ms.setOption(OPTION_RADICAL_MATCHING, RADICAL_MATCHING_EXACT);
        //ms.setHCountMatching(HCOUNT_MATCHING_EQUAL);
        ms.setQuery(query);

        return ms;
    }

    static String p (MolAtom a) {
        StringBuffer sb = new StringBuffer
            ("[a="+a.getAtno()
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
             );
        /*
          String[] props = a.getQPropNames();
          for (String p : props) {
          sb.append(","+p+"="+a.getQProp(p));
          }
        */
        sb.append("]");
        return sb.toString();
    }

    protected boolean _transform (Molecule mol)  {
        // use fingerprint to do quick filter
        MolHandler mh = new MolHandler (mol);
        int[] fp = mh.generateFingerprintInInts
            (FP_SIZE, FP_DEPTH, FP_LENGTH);
        for (int i = 0; i < fp.length; ++i) {
            if ((fp[i] & sourcefp[i]) != sourcefp[i]) {
                return false;
            }
        }

        int[] hit = null;

        lock.lock();
        try {
            ms.setTarget(mol);
            hit = ms.findFirst();
        }
        catch (SearchException ex) {
            logger.log(Level.SEVERE, "MolSearch failed", ex);
            return false;
        }
        finally {
            lock.unlock();
        }

        if (hit == null || hit.length == 0) {
            return false;
        }

        transform (mol, hit);
        mol.valenceCheck();

        return true;
    }

    void transform (Molecule mol, int[] hit) {
        if (debug) {
            System.out.println("Hit: "+smirks);
            for(int j=0; j < hit.length; j++) {
                System.out.print(" " + (hit[j]+1)+":"
                                 +source.getAtom(j).getAtomMap());
            }
            System.out.println();
        }

        if (target.getAtomCount() == 0) {
            // delete this subgraph
            for (int j = 0; j < hit.length; ++j) {
                mol.removeNode(hit[j]);
            }
        }
        else {
            MolAtom[] M = new MolAtom[mol.getAtomCount()];

            // update atom
            for (int j = 0; j < hit.length; ++j) {
                MolAtom src = target.getAtom(map[j]);
                MolAtom dst = mol.getAtom(hit[j]);
                //System.out.println("updating atom " + (hit[j]+1));
                copy (dst, src);
                M[hit[j]] = src;
            }
                    
            // now update the bond
            for (MolBond b1 : mol.getBondArray()) {
                MolAtom a1 = b1.getAtom1();
                MolAtom a2 = b1.getAtom2();
                int i1 = mol.indexOf(a1);
                int i2 = mol.indexOf(a2);

                if (M[i1] != null && M[i2] != null) {
                    for (int j = 0; j < M[i1].getBondCount(); ++j) {
                        MolBond b2 = M[i1].getBond(j);
                        if (b2.getOtherAtom(M[i1]) == M[i2]) {
                            int type = b2.getType();
                            switch (type) {
                            case MolBond.SINGLE_OR_AROMATIC:
                                b1.setType(1);
                                break;

                            case MolBond.DOUBLE_OR_AROMATIC:
                                b1.setType(2);
                                break;

                            case 1: case 2: case 3:
                                b1.setType(type);
                                break;
                            }

                            if (debug) {
                                System.out.println
                                    ("updating bond " + (i1+1)
                                     + "-"+(i2+1) + " to " 
                                     + b1.getType());
                            }
                        }
                    }
                }
            }

            M = null;
        }
    }

    public boolean transform (Molecule mol) {
        int txn = 0;
        while (_transform (mol)) {
            ++txn;
            if (debug) {
                logger.info("transform "+txn+": "+mol.toFormat("smiles:q"));
            }
        }
        return txn > 0;
    }

    static void copy (MolAtom dst, MolAtom src) {
        int map = dst.getAtomMap();
        int atno = dst.getAtno();
        int flags = dst.getFlags();
        int[] list = dst.getList();
        dst.set(src);
        // restoring values that were override by set() above
        dst.setAtno(atno);
        dst.setFlags(flags);
        if (list != null) {
            dst.setList(list);
        }
        //dst.setQuerystr(null);
        //dst.setSMARTS(null);
        dst.setAtomMap(map);
        dst.clearQProps();
    }

    public static void main (String[] argv) throws Exception {
        SMIRKS[] rules = {
            new SMIRKS ("[O-;v2:1]>>[O;+0:1]"),
        };

        MolImporter mi = new MolImporter (System.in);
        /*
        chemaxon.reaction.Reactor reactor = new chemaxon.reaction.Reactor();
        reactor.setTransform(true);
        reactor.setUnique(false);
        reactor.setSingle(false);
        */

        for (Molecule mol1; (mol1 = mi.read()) != null; ) {
            mol1.hydrogenize(false);
            mol1.dearomatize();

            Molecule mol2 = mol1.cloneMolecule();

            try {
                String smi = mol1.toFormat("smiles:q");
                java.util.Vector<SMIRKS> matches = 
                    new java.util.Vector<SMIRKS>();
                for (int i = 0; i < rules.length; ++i) {
                    /*
                    reactor.setReactionString(rules[i].getTransform());
                    reactor.setReactant(mol2);
                    reactor.react();
                    */

                    if (rules[i].transform(mol1)) {
                        /*
                          System.out.println("** " + rules[i].getTransform());
                          System.out.println(rules[i]);
                        */
                        matches.add(rules[i]);
                    }
                }
                
                String out1 = mol1.toFormat("smiles:q");
                String out2 = mol2.toFormat("smiles:q");
                if (!out1.equals(out2)) {
                    System.out.println(smi + "\t" + mol1.getName() 
                                       + "\t" + out1 + "\t" + out2);
                    for (SMIRKS smirks : matches) {
                        System.out.println
                            ("SMIRKS: " + smirks.getTransform());
                        System.out.println(smirks);
                    }
                }
                else {
                    System.out.println(smi + "\t" + mol1.getName() 
                                       + "\t" + out1);
                }
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }
}
