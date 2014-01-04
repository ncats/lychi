package lychi.tautomers;

import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
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
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.struc.StereoConstants;
import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.formats.MolImporter;
import chemaxon.util.MolHandler;

import lychi.TautomerGenerator;
import lychi.util.ChemUtil;

/**
 * tautomer code based on Roger Sayle and Jack Delany
 */
public class SayleDelanyTautomerGenerator implements TautomerGenerator {
    private static boolean debug = false;
    static {
	try {
	    debug = Boolean.getBoolean("lychi-tautomer.debug");
	}
	catch (Exception ex) {
	}
    }

    private static final Logger logger = Logger.getLogger
	(SayleDelanyTautomerGenerator.class.getName());

    static final int TAU_TYPE   = 0xf0;
    static final int TAU_OTHER  = 0x10;
    static final int TAU_DON    = 0x20;
    static final int TAU_ACC    = 0x40;
    static final int TAU_HYB    = 0x80;
    static final int TAU_OCSP3  = 0x11;
    static final int TAU_DCSP3  = 0x21;
    static final int TAU_ACSP2  = 0x42;
    static final int TAU_HCSP2  = 0x82;
    static final int TAU_NOXN   = 0x13;
    static final int TAU_OSP2   = 0x44;

    public static final int FLAG_NOXIDE = 1;
    public static final int FLAG_ALL = 7;


    static class IndexPair implements Comparable<IndexPair> {
	public int idx1, idx2;
	public IndexPair (int idx1, int idx2) {
	    this.idx1 = idx1;
	    this.idx2 = idx2;
	}
	public int compareTo (IndexPair ip) {
	    return idx2 - ip.idx2;
	}
    }

    // maximum number of tautomers to generate
    public static final int MAX_SIZE = 201;

    private int zoneCount; // zone count
    private int []dcounts, acounts; // donor/acceptor
    private int count, maxsize;

    private int flags = FLAG_NOXIDE;
    private long timeoutThreshold = 0;
    private long startTime;

    /*
     * input
     */
    private Molecule mol;
    private MolAtom[] atoms;
    private MolBond[] bonds;
    private int[] bflags; // bond flags
    private int[] amaps; // atom mapping

    /*
     * internal state
     */
    private int[] types;
    private int[] zones;
    private int[] avisit;
    private int[] bvisit;

    private List<Molecule> tautomers = new ArrayList<Molecule>();
    private Molecule canonicalTautomer = null;

    /**
     * Note that if the maximum number of tautomers is reached
     * (truncated is true), then the canonical tautomer is the first 
     * one.  Otherwise, the tautomer with the highest score is 
     * returned as the canonical tautomer!
     */
    private boolean truncated = false;

    public SayleDelanyTautomerGenerator () {
	this (MAX_SIZE);
    }

    public SayleDelanyTautomerGenerator (int maxsize) {
	setMaxTautomers (maxsize);
    }

    public SayleDelanyTautomerGenerator (int maxsize, int mask) {
	setMaxTautomers (maxsize);
	set (mask);
    }

    @Override
    public Object clone () {
	SayleDelanyTautomerGenerator taugen = 
            new SayleDelanyTautomerGenerator();
	taugen.flags = this.flags;
	taugen.timeoutThreshold = this.timeoutThreshold;
	taugen.maxsize = this.maxsize;

	return taugen;
    }

    public int getMaxTautomers () { return maxsize; }
    public void setMaxTautomers (int maxsize) { this.maxsize = maxsize; }

    public void setTimeout (long timeout) {
	this.timeoutThreshold = timeout;
    }
    public long getTimeout () { return timeoutThreshold; }

    public void setMolecule (Molecule mol) {
	init (mol);
    }
    public Molecule getMolecule () { return mol; }

    public void set (int mask) {
	flags |= mask;
    }
    public void unset (int mask) {
	flags &=~mask;
    }
    public boolean get (int mask) {
	return (flags & mask) == mask;
    }

    protected static boolean isKekulized (Molecule m) {
	for (MolBond b : m.getBondArray()) {
	    if (b.getType() == MolBond.AROMATIC) {
		return false;
	    }
	}
	return true;
    }

    private static boolean kekulize (Molecule out, Molecule mol) {
	Molecule m = mol.cloneMolecule();
	m.aromatize(Molecule.AROM_BASIC);
	m.dearomatize();

	boolean kekulized = isKekulized (m);
	if (kekulized) {
	    if (out != null) {
		out.clear();
		m.clonecopy(out);
	    }
	}
	return kekulized;
    }

    protected void init (Molecule m) {
	// first make a copy of the input if it's kekulized...
	//this.mol = isKekulized (m) ? m.cloneMolecule() : m;
	Molecule _m = m;
	m = m.cloneMolecule();
	m.clearProperties();

	String name = m.getName();
	m.expandSgroups();
	m.hydrogenize(false);
	m.aromatize();
	m.dearomatize();

	String smiles = null;
	try {
	    /*
	     * For whatever reason, we have to do this because the internal
	     * state of the input molecule can be quite screwed up.
	     * This allows us to start with a clean slate molecule.
	     */
            Map<Integer, Integer> chiral = new HashMap<Integer, Integer>();
            atoms = m.getAtomArray();
            //logger.info("Tautomer: ++ "+m.toFormat("smiles:q"));
            for (int i = 0; i < atoms.length; ++i) {
                int k = atoms[i].getAtomMap();
                int c = m.getChirality(i);
                if (c == MolAtom.CHIRALITY_R || c == MolAtom.CHIRALITY_S) {
                    //logger.info(k+": "+c);
                    chiral.put(k, c);
                }
            }

	    smiles = m.toFormat("mol");
	    if (debug) {
		logger.info("Input molecule: " + smiles);
	    }

	    MolImporter.importMol(smiles, m);
	    m.dearomatize();

            // restore the chirality flags
            atoms = m.getAtomArray();
            for (int i = 0; i < atoms.length; ++i) {
                Integer c = chiral.get(atoms[i].getAtomMap());
                if (c != null) {
                    //logger.info(atoms[i].getAtomMap()+": "+c);
                    m.setChirality(i, c);
                }
            }
            //logger.info("Tautomer: -- "+m.toFormat("smiles:q"));

	    if (!isKekulized (m)) {
		// can't kekulized this molecule, so revert back to the 
		//   original input
		m.clear();
		_m.clonecopy(m);
	    }

	    m.setName(name);
	}
	catch (Exception ex) {
	    logger.log(Level.SEVERE, "Invalid input structure: "+smiles, ex);
	    throw new IllegalArgumentException 
		("Invalid input structure: "+smiles);
	}
	this.mol = m;

	atoms = mol.getAtomArray();
	amaps = new int[atoms.length];
	/*
	 * TODO: need to verify this
	 * for whatever lame reason, there are instances where having atom
	 * mapping messes up the way tautomer are generated!
	 * e.g., CCN1C2=NNC(=C2C3=NN=C(C)N3C4=C1C=CC=C4)C5=CC=CC=C5
	 * CID 383802
	 */
	for (int i = 0; i < atoms.length; ++i) {
	    MolAtom a = atoms[i];
	    amaps[i] = a.getAtomMap();
	    //a.setAtomMap(0); // now reset it
	}

	bonds = mol.getBondArray();
	bflags = new int[bonds.length];
	for (int i = 0; i < bonds.length; ++i) {
	    MolBond bnd = bonds[i];
	    int flags = bnd.getFlags();
	    if (bnd.getType() == MolBond.AROMATIC) {
		int a = flags & MolBond.TYPE_MASK;
		int b = flags & MolBond.STEREO_MASK;
		int c = flags & MolBond.TOPOLOGY_MASK;
		int d = flags & MolBond.REACTING_CENTER_MASK;
		String str = ("bond " + (this.mol.indexOf(bnd.getAtom1())+1)
			      + "-" + (this.mol.indexOf(bnd.getAtom2())+1)
			      + " type="+a + " stereo="+b 
			      + " topo="+c + " chiral="+d);

		throw new IllegalArgumentException 
		    ("Input structure can't be kekulized: "+smiles+" "+str);
	    }

	    if (debug) {
		if (i == 0) {
		    logger.info("## Init bond annotations for "
				+mol.toFormat("smiles:q"));
		}

		int stereo = flags & MolBond.STEREO_MASK;
		int topo = flags & MolBond.TOPOLOGY_MASK;
		int chiral = flags & MolBond.REACTING_CENTER_MASK;
		if (stereo != 0) {
		    String ez = "";
		    if ((stereo & StereoConstants.TRANS) 
			== StereoConstants.TRANS) {
			ez += "E";
		    }
		    if ((stereo & StereoConstants.CIS) 
			== StereoConstants.CIS) {
			if (ez.length() > 0) {
			    ez += "/";
			}
			ez += "Z";
		    }

		    logger.info("Bond "+(i+1)+":"
				+ (mol.indexOf(bnd.getAtom1())+1)
				+ (bnd.getType()==2?"=":"-")
				+ (mol.indexOf(bnd.getAtom2())+1)
				+ " has stereo flags: "+ez + " ("+stereo+")");
		}
		if (topo != 0) {
		    logger.info("Bond "+(i+1)+":"
				+ (mol.indexOf(bnd.getAtom1())+1)
				+ (bnd.getType()==2?"=":"-")
				+ (mol.indexOf(bnd.getAtom2())+1)
				+" has stereo flags: "+topo);
		}
		if (chiral != 0) {
		    logger.info("Bond "+(i+1)+":"
				+ (mol.indexOf(bnd.getAtom1())+1)
				+ (bnd.getType()==2?"=":"-")
				+ (mol.indexOf(bnd.getAtom2())+1)
				+" has stereo flags: "+chiral);
		}
	    }
	    bflags[i] = flags;
	}
	//System.out.print(this.mol.toFormat("sdf:-a"));

	types = new int[atoms.length];
	zones = new int[atoms.length];
	avisit = new int[atoms.length];
	bvisit = new int[bonds.length];

	tautomers.clear();
	canonicalTautomer = null;
    }

    public int generate (Molecule molecule) {
        init (molecule);

	assignAtomTypes ();
	if (debug) {
	    System.out.println("assignAtomTypes()");
	    displayAtomTypes ();
	}

	if (get (FLAG_NOXIDE)) {
	    assignNOxideAtomTypes ();
	    if (debug) {
		System.out.println("assignNOxideAtomTypes()");
		displayAtomTypes ();
	    }
	}

	if (get (FLAG_ALL)) {
	    assignAdditionalAtomTypes ();
	    if (debug) {
		System.out.println("assignAdditionalAtomTypes()");
		displayAtomTypes ();
	    }
	}

	initializeVisitFlags ();
	initializeTautomer ();

	/*
	  System.err.println("initializeTautomer()");
	  displayAtomTypes ();
	*/

	zoneCount = labelZones ();

	// don't bother if we don't have any h- donor and/or acceptor
	if (0 == countZoneTypes ()) {
	    return 0; 
	}

	fixZones ();
	count = countTotalDonors ();
	if (count == 0) {
	    return 0;
	}

	//displayZoneCounts ();
	normalization ();

	/*
	  System.err.println("normalization()");
	  displayAtomTypes ();
	*/

	//recurseTautomers (2, null);
	startTime = System.currentTimeMillis();
	recurseTautomers (2, graphInvariantOrder (mol));

	if (debug) {
	    logger.info("## Final bond annotations");

	    for (int i = 0; i < bflags.length; ++i) {
		int flags = bflags[i];
		MolBond bnd = mol.getBond(i);

		int stereo = flags & MolBond.STEREO_MASK;
		int topo = flags & MolBond.TOPOLOGY_MASK;
		int chiral = flags & MolBond.REACTING_CENTER_MASK;
		if (stereo != 0) {
		    String ez = "";
		    if ((stereo & StereoConstants.TRANS) 
			== StereoConstants.TRANS) {
			ez += "E";
		    }
		    if ((stereo & StereoConstants.CIS) 
			== StereoConstants.CIS) {
			if (ez.length() > 0) {
			    ez += "/";
			}
			ez += "Z";
		    }

		    logger.info("Bond "+bnd.getAtom1().getSymbol()
				+ (bnd.getType()==2?"=":"-")
				+ bnd.getAtom2().getSymbol()
				+ " has stereo flags: "+ez + " ("+stereo+")");
		}
		if (topo != 0) {
		    logger.info("Bond "+bnd.getAtom1().getSymbol()
				+ (bnd.getType()==2?"=":"-")
				+ bnd.getAtom2().getSymbol()
				+" has stereo flags: "+topo);
		}
		if (chiral != 0) {
		    logger.info("Bond "+bnd.getAtom1().getSymbol()
				+ (bnd.getType()==2?"=":"-")
				+ bnd.getAtom2().getSymbol()
				+" has stereo flags: "+chiral);
		}
	    }
	}

        return tautomers.size();
    }

    public Enumeration<Molecule> tautomers () {
	return Collections.enumeration(tautomers);
    }

    public int getTautomerCount () { return tautomers.size(); }

    public String getCanonicalTautomerAsString () { 
	return ChemUtil.canonicalSMILES(getCanonicalTautomer());
    }

    /*
     * since the tautomers are always generated in the same order
     * (per the graph invariant), the canonical tautomer can simply
     * be the first tautomer found.  however, here we apply a set of
     * scoring heuristics to find the best "looking" one instead.
     */
    public synchronized Molecule getCanonicalTautomer () { 
	if (canonicalTautomer == null) {
	    if (tautomers.isEmpty()) {
		canonicalTautomer = mol;
	    }
	    else if (/*truncated ||*/ tautomers.size() == 1) {
		// no point evaluating the score if we only have on tautomer
		canonicalTautomer = tautomers.iterator().next();
	    }
	    else {
		int bestScore = -1;
		List<Molecule> candidates = new ArrayList<Molecule>();
		for (Molecule tau : tautomers) {
		    int score = scoreTautomer (tau);
		    if (score > bestScore || canonicalTautomer == null) {
			canonicalTautomer = tau;
			candidates.add(tau);
			bestScore = score;
		    }
		    else if (score == bestScore) {
			candidates.add(tau);
		    }
		}

		/*
		if (candidates.size() > 1) {
		    if (debug) {
			logger.info(candidates.size() 
				    + " tautomers with same score: "
				    +bestScore);
		    }

		    // sort all candidates with the same score based on
		    //   (reverse) alphabetical order
		    Collections.sort
			(candidates, new Comparator<Molecule>() {
			    public int compare (Molecule m1, Molecule m2) {
				String s1=MolStandardizer.canonicalSMILES(m1);
				String s2=MolStandardizer.canonicalSMILES(m2);
				return s2.compareTo(s1);
			    }
			});
		    canonicalTautomer = candidates.iterator().next();
		}
		*/
	    }
	}

	return canonicalTautomer;	
    }

    protected static int[] graphInvariantOrder (Molecule m) {
	int[] gi = new int[m.getAtomCount()];
	m.getGrinv(gi, Molecule.GRINV_STEREO|Molecule.GRINV_NOHYDROGEN);

	IndexPair ip[] = new IndexPair[gi.length];
	for (int i = 0; i < gi.length; ++i) {
	    ip[i] = new IndexPair (i, gi[i]);
	}
	Arrays.sort(ip);
	for (int i = 0; i < gi.length; ++i) {
	    gi[i] = ip[i].idx1;
	}
	ip = null;
	return gi;
    }

    protected static void debugStructure (Molecule m) {
	MolAtom[] atoms = m.getAtomArray();
	for (int i = 0; i < atoms.length; ++i) {
	    atoms[i].setAtomMap(i+1);
	}
	System.out.println("  " + m.toFormat("smiles:ua_bas"));
	MolBond[] bonds = m.getBondArray();
	for (int i = 0; i < bonds.length; ++i) {
	    System.out.printf("   %1$2d: %2$2d %3$2d %4$d %5$d\n",
			      i, m.indexOf(bonds[i].getAtom1())+1,
			      m.indexOf(bonds[i].getAtom2())+1,
			      bonds[i].getType(), 
			      bonds[i].getFlags());
	}
    }

    static public int scoreTautomer (Molecule mol) {
	return ChemUtil.calcMolScore(mol);
    }

    protected void normalization () {
	for (int i = 0; i < bonds.length; ++i) {
	    if (bvisit[i] == 0) {
		setBondOrder (bonds[i], 1);
	    }
	}

	for (int i = 0; i < atoms.length; ++i) {
	    MolAtom a = atoms[i];
	    if ((types[i] & TAU_DON) != 0) {
		a.setImplicitHcount(a.getImplicitHcount()+3);
	    }
	    else if ((types[i] & TAU_OTHER) == 0) {
		a.setImplicitHcount(a.getImplicitHcount()+4);
	    }
	}
    }

    protected boolean assignAtomTypes () {
	int acc = 0, don = 0; // count of donor & acceptor

	for (int i = 0; i < atoms.length; ++i) {
	    MolAtom atom = atoms[i];

	    int charge = atom.getCharge();
	    int hcount = atom.getImplicitHcount();
	    int sb = hcount, db = 0, tb = 0;

	    for (int j = 0; j < atom.getBondCount(); ++j) {
		MolBond bond = atom.getBond(j);
		switch (bond.getType()) {
		case 3: ++tb; break;
		case 2: ++db; break;
		case 1: ++sb; break;
		}
	    }

	    int type = TAU_OTHER;
	    switch (atom.getAtno()) {
	    case 6: /* Carbon */
		if (charge == 0) {
		    /* "[#6X3v4+0]=*" */
		    if ((sb==2) && (db==1) && (tb==0)) {
			type = TAU_HCSP2;
		    } 
		    else if ((sb==4) && (db==0) && (tb==0)) {
			if (hcount > 0) 
			    type = TAU_OCSP3;
		    }
		} 
		else if (charge == -1) {
		    /* "[#6X2v3-1]=*" */
		    if ((sb==1) && (db==1) && (tb==0)) {
			type = TAU_HCSP2;
		    }
		} 
		else if (charge == 1) {
		    /* "[#6X2v3+1]=*" */
		    if ((sb==1) && (db==1) && (tb==0)) {
			type = TAU_HCSP2;
		    }
		}
		break;

	    case 7: /* Nitrogen */
		if (charge == 0) {
		    /* "[#7X2v3+0]=*" */
		    if ((sb==1) && (db==1) && (tb==0)) {
			type = TAU_ACC;
			++acc;
			/* "[#7X3v3!H0+0]" */
		    } 
		    else if ((sb==3) && (db==0) && (tb==0)) {
			if (hcount >0) {
			    type = TAU_DON;
			    ++don;
			}
			/* "[#7X3v5+0](=*)=*" */
		    } 
		    else if ((sb==1) && (db==2) && (tb==0)) {
			type = TAU_NOXN;
		    }
		} 
		else if (charge == 1) {
		    /* "[#7X3v4+0]=*" */
		    if ((sb==2) && (db==1) && (tb==0)) {
			type = TAU_HYB;
		    }
		}
		else if (charge == -1) {
		    if ((sb == 2) && (db == 0) && (tb == 0)) {
			type = TAU_DON;
			atom.setImplicitHcount(1);
			atom.setCharge(0);
			++don;
		    }
		}
		break;

	    case 8:  /* Oxygen   */
	    case 16: /* Sulphur  */
	    case 34: /* Selenium */
	    case 52: /* Telurium */
		/* "[#?? X1v2+0]" */
		if (charge == 0) {
		    if ((sb==0) && (db==1) && (tb==0)) {
			type = TAU_OSP2;
			/* "[#?? X2v2!H0+0]" */
		    } 
		    else if ((sb==2) && (db==0) && (tb==0)) {
			if (hcount > 0) {
			    type = TAU_DON;
			    ++don;
			}
		    }
		}
		break;
	    }
	    types[i] = type;
	}

	return acc + don > 0;
    }

    protected void assignNOxideAtomTypes () {
	for (int i = 0; i < bonds.length; ++i) {
	    MolBond bond = bonds[i];
	    
	    MolAtom src = bond.getAtom1();
	    MolAtom dst = bond.getAtom2();

	    int sidx = mol.indexOf(src);
	    int didx = mol.indexOf(dst);
	    
	    int stype = types[sidx];
	    int dtype = types[didx];
	    
	    if ((stype==TAU_NOXN) && (dtype==TAU_OSP2)) {
		setBondOrder (bond, 1);
		types[sidx] = TAU_HYB;
		src.setCharge(1);
		types[didx] = TAU_OTHER;
		dst.setCharge(-1);
	    } 
	    else if ((stype==TAU_OSP2) && (dtype==TAU_NOXN)) {
		setBondOrder (bond, 1);
		types[didx] = TAU_HYB;
		dst.setCharge(1);
		types[sidx] = TAU_OTHER;
		src.setCharge(-1);
	    }
	}
    }


    protected void assignAdditionalAtomTypes () {
	for (int i = 0; i < atoms.length; ++i) {
	    if (types[i] == TAU_HCSP2
		/*&& atoms[i].getImplicitHcount() == 0*/) {
		types[i] = TAU_ACSP2;
	    }
	}

	int type;
	for (int i = 0; i < bonds.length; ++i) {
	    if (mol.isRingBond(i)) {
		MolAtom src = bonds[i].getAtom1();
		MolAtom dst = bonds[i].getAtom2();
		int six = mol.indexOf(src), dix = mol.indexOf(dst);
		if (types[six] == TAU_OCSP3) {
		    type = types[dix];
		    if ((type & TAU_OTHER) == 0 && type != TAU_DCSP3)
			types[six] = TAU_DCSP3;
		}
		if (types[dix] == TAU_OCSP3) {
		    type = types[six];
		    if ((type & TAU_OTHER) == 0 && type != TAU_DCSP3)
			types[dix] = TAU_DCSP3;
		}
	    }
	}
    }

    protected void initializeVisitFlags () {
	for (int i = 0; i < atoms.length; ++i) {
	    if ((types[i] & (TAU_ACC | TAU_DON)) != 0) {
		avisit[i] = 0;
	    }
	    else {
		avisit[i] = 1;
	    }
	}

	for (int i = 0; i < bonds.length; ++i) {
	    MolBond bond = bonds[i];

	    MolAtom src = bond.getAtom1();
	    MolAtom dst = bond.getAtom2();
	    if (((types[mol.indexOf(src)] & TAU_OTHER) == 0)
		&& ((types[mol.indexOf(dst)] & TAU_OTHER) == 0)) {
		bvisit[i] = 0;
	    }
	    else {
		bvisit[i] = 1;
	    }
	}
    }

    protected void initializeTautomer () {
	for (int i = 0; i < atoms.length; ++i) {
	    if ((types[i] & TAU_OTHER) == 0) {
		initializeTautomerAtom (i, atoms[i]);
	    }
	}
    }

    protected void initializeTautomerAtom (int ai, MolAtom atom) {
	if (ai < 0) {
	    ai = mol.indexOf(atom);
	}

	//System.out.print("InitializeTautomerAtom: Atom " + ai + "...");

	if (mustHaveSingle (ai)) {
	    //System.out.println("single");
	    types[ai] = TAU_OTHER;
	    avisit[ai] = 1;
	    for (int i = 0; i < atom.getBondCount(); ++i) {
		MolBond bond = atom.getBond(i);
		int bi = mol.indexOf(bond);
		if (bvisit[bi] == 0) {
		    setBondOrder (bond, 1);
		    bvisit[bi] = 1;
		    initializeTautomerAtom (-1, bond.getOtherAtom(atom));
		}
	    }
	}
	else if (mustHaveDouble (ai)) {
	    //System.out.println("double");
	    types[ai]  = TAU_OTHER;
	    avisit[ai] = 1;
	    for (int i = 0; i < atom.getBondCount(); ++i) {
		MolBond bond = atom.getBond(i);
		int bi = mol.indexOf(bond);
		if (bvisit[bi] == 0) {
		    setBondOrder (bond, 2);
		    bvisit[bi] = 1;
		    initializeTautomerAtom (-1, bond.getOtherAtom(atom));
		}
	    }
	}
	else {
	    //System.out.println("other");
	    int avail = 0;
	    for (int i = 0; i < atom.getBondCount(); ++i) {
		MolBond bond = atom.getBond(i);
		int bi = mol.indexOf(bond);
		if (bvisit[bi] == 0) {
		    ++avail;
		}
	    }

	    if (avail == 0) {
		types[ai] =  TAU_OTHER;
		avisit[ai] = 1;
	    }
	}
	//System.out.println("Atom " + ai);
	//displayAtomTypes ();
    }

    protected int labelZones () {
	int zone = 1;
	for (int i = 0; i < atoms.length; ++i) {
	    if ((types[i] & TAU_OTHER) == 0 && (zones[i] == 0)) {
		dfsZone (i, zone++);
	    }
	}
	return zone - 1;
    }

    protected int countZoneTypes () {
	acounts = new int[zoneCount+1];
	dcounts = new int[zoneCount+1];

	for (int i = 0; i < atoms.length; ++i) {
	    int z = zones[i];
	    switch (types[i] & TAU_TYPE) {
	    case TAU_ACC: acounts[z]++; break;
	    case TAU_DON: dcounts[z]++; break;
	    }
	}

	int total = 0;
	for (int z = 1; z <= zoneCount; ++z) {
	    total += acounts[z] + dcounts[z];
	}
	return total;
    }

    protected int countTotalDonors () {
	int total = 0;
	for (int i = 1; i <= zoneCount; ++i) {
	    total += dcounts[i];
	}
	return total;
    }

    protected int countTotalAcceptors () {
	int total = 0;
	for (int i = 1; i <= zoneCount; ++i) {
	    total += acounts[i];
	}
	return total;
    }	

    protected void fixZones () {
	for (int i = 1; i <= zoneCount; ++i) {
	    if (maxsize > 0 && (acounts[i] + dcounts[i]) > maxsize) {
		fixZone (i);
	    }
	    else if (dcounts[i] == 0) {
		fixZone (i);
	    }
	}
    }

    protected void fixZone (int zone) {
	for (int i = 0; i < atoms.length; ++i) {
	    if (zones[i] == zone) {
		avisit[i] = 1;
	    }
	}
	acounts[zone] = 0;
	dcounts[zone] = 0;
    }

    protected void dfsZone (int ai, int zone) {
	MolAtom atom = atoms[ai];
	zones[ai] = zone;
	for (int i = 0; i < atom.getBondCount(); ++i) {
	    MolAtom xatom = atom.getBond(i).getOtherAtom(atom);
	    int ix = mol.indexOf(xatom);
	    if ((types[ix] & TAU_OTHER) == 0 && (zones[ix] == 0)) {
		dfsZone (ix, zone);
	    }
	}
    }

    protected boolean mustHaveSingle (int ai) {
	if (avisit[ai] != 0) {
	    if ((types[ai] & TAU_DON) != 0)
		return true;
	}

	MolAtom atom = atoms[ai];
	for (int i = 0; i < atom.getBondCount(); ++i) {
	    MolBond bond = atom.getBond(i);
	    if (bond.getType() == 2 && bvisit[mol.indexOf(bond)] != 0) {
		return true;
	    }
	}

	return false;
    }

    protected boolean mustHaveDouble (int ai) {
	if ((types[ai] & TAU_DON) == 0 && avisit[ai] != 0) {
	    int doub = 0, avail = 0;
	    MolAtom atom = atoms[ai];
	    for (int i = 0; i < atom.getBondCount(); ++i) {
		MolBond bond = atom.getBond(i);
		if (bvisit[mol.indexOf(bond)] != 0) {
		    if (bond.getType() == 2)
			++doub;
		}
		else {
		    ++avail;
		}
	    }
	    return avail== 1 && doub == 0;
	}
	return false;
    }

    protected boolean recurseTautomers (int depth, int[] tauOrder) {
	if (debug) {
	    System.out.println("Recurse(depth="+depth+"): count="+count );
	    displayZoneCounts ();
	    //displayAtomTypes ();
	    String smiles = mol.toFormat("smiles:q");
	    System.out.println(smiles);
	}

	int ix = -1;

	if (tauOrder != null) {
	    for (int i = 0; i < tauOrder.length; ++i) {
		if (avisit[tauOrder[i]] == 0) {
		    ix = tauOrder[i];
		    break;
		}
	    }
	}
	else {
	    for (int i = 0; i < atoms.length; ++i) {
		if (avisit[i] == 0) {
		    ix = i;
		    break;
		}
	    }
	}

	if (ix < 0) {
	    if (count == 0 && kekuleTautomer (depth))
		return generateTautomer ();
	    return false;
	}

	// check for timeout...
	if (timeoutThreshold > 0) {
	    long diff = (System.currentTimeMillis() - startTime);
	    if (diff  > timeoutThreshold) {
		System.err.println("** WARNING: " + getClass().getName() 
				   + ": max timeout ("
				   + String.format("%1$.2fs", diff/1000.)
				   + ") reached; search truncated!");
		return true;
	    }
	}

	int zone = zones[ix];
	avisit[ix] = depth;

	int result, rc;
	if (dcounts[zone] > 0) {
	    if ((types[ix] & TAU_ACC) != 0)
		acceptProton (ix);
	    result = propagateTautomer (atoms[ix], depth + 1);
	    if ((result >= 0) && (result < dcounts[zone])) {
		++result;
		dcounts[zone] -= result;
		count -= result;
		if (dcounts[zone] == 0) {
		    rc = kekuleTautomerZone (depth + 2, zone);
		}
		else 
		    rc = depth + 2;
		if (rc != 0) {
		    rc = recurseTautomers (rc, tauOrder) ? 1 : 0;
		}
		dcounts[zone] += result;
		count += result;
		if (rc != 0)
		    return true;
	    }

	    donateProton (ix);
	    resetVisitFlags (depth + 1);
	    result = propagateTautomer (atoms[ix], depth + 1);
	    if ((result >= 0) && (result <= dcounts[zone])) {
		dcounts[zone] -= result;
		count -= result;
		if (dcounts[zone] == 0) {
		    rc = kekuleTautomerZone (depth + 2, zone);
		}
		else 
		    rc = depth + 2;
		if (rc != 0) {
		    rc =  recurseTautomers (rc, tauOrder) ? 1: 0;
		}
		dcounts[zone] += result;
		count += result;
		if (rc != 0)
		    return true;
	    }
	}
	else {
	    if ((types[ix] & TAU_DON) != 0) 
		donateProton (ix);
	    result = propagateTautomer (atoms[ix], depth + 1);
	    if (result != 0) return false;
	    rc = kekuleTautomerZone (depth + 2, zone);
	    if (rc == 0) return false;
	    return recurseTautomers (rc, tauOrder);
	}

	return false;
    }

    protected void acceptProton (int ix) {
	MolAtom atom = atoms[ix];
	atom.setImplicitHcount(atom.getImplicitHcount()+1);
	types[ix] = TAU_DON;
    }

    protected void donateProton (int ix) {
	MolAtom atom = atoms[ix];
	int h = atom.getImplicitHcount();
	if (h > 0) {
	    atom.setImplicitHcount(h-1);
	}
	types[ix] = TAU_ACC;
    }

    protected void resetVisitFlags (int depth) {
	for (int i = 0; i < atoms.length; ++i) {
	    if (avisit[i] >= depth) {
		avisit[i] = 0;
	    }
	}

	for (int i = 0; i < bonds.length; ++i) {
	    if (bvisit[i] >= depth) {
		bvisit[i] = 0;
	    }
	}
    }

    protected boolean kekuleTautomer (int depth) {
	for (int i = 0; i < bonds.length; ++i) {
	    if (bvisit[i] == 0) {
		MolBond bnd = bonds[i];
		MolAtom src = bnd.getAtom1();
		MolAtom dst = bnd.getAtom2();

		bvisit[i] = depth;
		setBondOrder (bnd, 2);

		if (propagateTautomer (src, depth+1) != 0 ||
		    propagateTautomer (dst, depth+1) != 0) {
		    setBondOrder (bnd, 1);
		    resetVisitFlags (depth + 1);
		    
		    if (propagateTautomer (src, depth + 1) != 0 ||
			propagateTautomer (dst, depth + 1) != 0)
			return false;
		}
		depth += 2;
	    }
	}
	return true;
    }

    protected int kekuleTautomerZone (int depth, int zone) {
	for (int i = 0; i < atoms.length; ++i) {
	    if (avisit[i] == 0 && zones[i] == zone) {
		avisit[i] = depth;
		if ((types[i] & TAU_DON) != 0)
		    donateProton (i);
		if (propagateTautomer (atoms[i], depth+1) != 0)
		    return 0;
		depth += 2;
	    }
	}

	for (int i = 0; i < bonds.length; ++i) {
	    if (bvisit[i] == 0) {
		MolBond bnd = bonds[i];
		MolAtom src = bnd.getAtom1();
		MolAtom dst = bnd.getAtom2();
		
		if (zones[mol.indexOf(src)] == zone) {
		    bvisit[i] = depth;
		    setBondOrder (bnd, 2);
		    if (propagateTautomer (src, depth+1) != 0 ||
			propagateTautomer (dst, depth+1) != 0) {
			setBondOrder (bnd, 1);
			resetVisitFlags (depth+1);
			if (propagateTautomer (src, depth+1) != 0 ||
			    propagateTautomer (dst, depth+1) != 0)
			    return 0;
		    }
		    depth += 2;
		}
	    }
	}
	return depth;
    }

    protected int propagateTautomer (MolAtom atom, int depth) {
	int doub = 0, avail = 0;
	for (int i = 0; i < atom.getBondCount(); ++i) {
	    MolBond bnd = atom.getBond(i);
	    int ix = mol.indexOf(bnd);
	    if (bvisit[ix] != 0) {
		if (bnd.getType() == 2)
		    ++doub;
	    }
	    else {
		++avail;
	    }
	}
	if (doub > 1)
	    return -1;

	int result, stat;
	int ix = mol.indexOf(atom);
	if (avisit[ix] != 0) {
	    /* Donor/Acceptor status confirmed! */
	    if ((types[ix] & TAU_DON) != 0) {
		if (doub > 0) return -1;

		if (avail > 0) {
		    /* MustHaveSingle */
		    result = 0;
		    for (int i = 0; i < atom.getBondCount(); ++i) {
			MolBond bnd = atom.getBond(i);
			int bx = mol.indexOf(bnd);
			if (bvisit[bx] == 0) {
			    MolAtom xa = bnd.getOtherAtom(atom);
			    setBondOrder (bnd, 1);
			    bvisit[bx] = depth;
			    stat = propagateTautomer (xa, depth);
			    if (stat < 0)
				return -1;
			    result += stat;
			}
		    }
		    return result;
		}
	    }
	    else {
		if (avail > 0) {
		    if (doub > 0) {
			// MustHaveSingle 
			result = 0;
			for (int i = 0; i < atom.getBondCount(); ++i) {
			    MolBond bnd = atom.getBond(i);
			    int bx = mol.indexOf(bnd);
			    if (bvisit[bx] == 0) {
				MolAtom xa = bnd.getOtherAtom(atom);
				setBondOrder (bnd, 1);
				bvisit[bx] = depth;
				stat = propagateTautomer (xa, depth);
				if (stat < 0)
				    return -1;
				result += stat;
			    }
			}
			return result;
		    }
		    else if (avail == 1) {
			// MustHaveDouble
			result = 0;
			for (int i = 0; i < atom.getBondCount(); ++i) {
			    MolBond bnd = atom.getBond(i);
			    int bx = mol.indexOf(bnd);
			    if (bvisit[bx] == 0) {
				MolAtom xa = bnd.getOtherAtom(atom);
				setBondOrder (bnd, 2);
				bvisit[bx] = depth;
				stat = propagateTautomer (xa, depth);
				if (stat < 0)
				    return -1;
				result += stat;
			    }
			}
			return result;
		    }
		    else return 0;
		}
		else if (doub == 0)
		    return -1;
	    }
	}
	else {
	    if (doub > 0) {
		// Must be acceptor!
		avisit[ix] = depth;
		if ((types[ix] & TAU_DON) != 0)
		    donateProton (ix);
		result = 0;
		if (avail > 0) {
		    // MustHaveSingle
		    for (int i = 0; i < atom.getBondCount(); ++i) {
			MolBond bnd = atom.getBond(i);
			int bx = mol.indexOf(bnd);
			if (bvisit[bx] == 0) {
			    MolAtom xa = bnd.getOtherAtom(atom);
			    setBondOrder (bnd, 1);
			    bvisit[bx] = depth;
			    stat = propagateTautomer (xa, depth);
			    if (stat < 0)
				return -1;
			    result += stat;
			}
		    }
		}
		return result;
	    }
	    else if (avail == 0) {
		// Must be donor!
		avisit[ix] = depth;
		if ((types[ix] & TAU_ACC) != 0)
		    acceptProton (ix);
		return 1;
	    }
	}

	return 0;
    }

    protected static void setBondOrder (MolBond bnd, int order) {
	/*
	 * Note:  setFlags must be used here because, for whatever reason,
	 * using setType affects the order in which tautomers are generated;
	 * that is, tautomers are no longer generated in a canonical order.
	 */
	bnd.setFlags(order, MolBond.TYPE_MASK);
    }

    protected boolean generateTautomer () {
	if (tautomers.size() < maxsize) {
	    mol.valenceCheck();
	    adjustBondAnnotations ();
	    Molecule tau = mol.cloneMolecule(); 
            //for (int i = 0; i < atoms.length; ++i) {
                //tau.getAtom(i).setAtomMap(amaps[i]);
            //}

	    if (debug) {
		logger.info("Tautomer generated: "+tau.toFormat("smiles:q"));
	    }

	    tautomers.add(tau);
	}
	
	truncated = tautomers.size() >= maxsize;
	if (truncated) {
	    if (debug) {
		logger.log(Level.WARNING, "Max tautomers reached ("+maxsize
			   +") for " + mol.getName() +"; "
			   +"canonical tautomer might not have the best score!");
	    }
	}
	return truncated;
    }

    protected void adjustBondAnnotations () {

	for (int i = 0; i < bonds.length; ++i) {
	    MolBond bnd = bonds[i];
	    int flags = bflags[i];

	    int type = bnd.getFlags() & MolBond.TYPE_MASK;
	    int order = flags & MolBond.TYPE_MASK;
	    int stereo = flags & MolBond.STEREO_MASK;
	    int topo = flags & MolBond.TOPOLOGY_MASK;
	    int chiral = flags & MolBond.REACTING_CENTER_MASK;
	    if (order == type) {
		// restore the flags if the bond type is the same
		bnd.setFlags(stereo, MolBond.STEREO_MASK);
		bnd.setFlags(topo, MolBond.TOPOLOGY_MASK);
		bnd.setFlags(chiral, MolBond.REACTING_CENTER_MASK);
	    }
	    else if (type == 1) {
		/* 
		 * here we reset all stereo flags of a mobile hydrogen
		 */
		if (order != type) {
		    // make sure previously generated tautomers turn off
		    //  the stereo flag for this bond
		    for (Molecule tau : tautomers) {
			tau.getBond(i).setFlags(0, MolBond.STEREO_MASK);
			tau.getBond(i).setFlags(0, MolBond.TOPOLOGY_MASK);
			tau.getBond(i).setFlags
			    (0, MolBond.REACTING_CENTER_MASK);
		    }
		    bflags[i] &= ~MolBond.STEREO_MASK;
		}
		bnd.setFlags(0, MolBond.STEREO_MASK);
		bnd.setFlags(0, MolBond.TOPOLOGY_MASK);
		bnd.setFlags(0, MolBond.REACTING_CENTER_MASK);
	    }

	    if (debug) {
		if (stereo != 0 && order != type) {
		    logger.info
			("** warning: "
			 +"stereo flag (" + stereo 
			 +") is ignored due to bond type changing from " 
			 + order + " to " + bnd.getType() 
			 + " at bond " +(i+1)+":"
			 + (mol.indexOf(bnd.getAtom1())+1) 
			 + "-" + (mol.indexOf(bnd.getAtom2())+1));
		}
		if (topo != 0) {
		    logger.info
			("** warning: "
			 +"topology flag (" + topo 
			 +") is ignored due to bond type changing from " 
			 + order + " to " + bnd.getType() 
			 + " at bond " + (i+1)+":"
			 + (mol.indexOf(bnd.getAtom1())+1) 
			 + "-" + (mol.indexOf(bnd.getAtom2())+1));
		}
		if (chiral != 0) {
		    logger.info
			("** warning: "
			 +"chiral flag (" + chiral
			 +") is ignored due to bond type changing from " 
			 + order + " to " + bnd.getType() 
			 + " at bond " + (i+1)+":"
			 + (mol.indexOf(bnd.getAtom1())+1) 
			 + "-" + (mol.indexOf(bnd.getAtom2())+1));
		}
	    }
	}
    }


    protected void displayZoneCounts () {
	int acount = 0,  dcount = 0;
	int hcount = 0,  ocount = 0;

	for (int i = 0; i < atoms.length; ++i) {
	    switch (types[i] & TAU_TYPE ) {
            case TAU_ACC: acount++;  break;
            case TAU_DON:    dcount++;  break;
            case TAU_HYB:  hcount++;  break;
            default: ocount++;  break;
	    }
	}
	
	System.out.printf("Acceptors ... %1$d\n",acount);
	System.out.printf("Donors ...... %1$d\n",dcount);
	System.out.printf("Hybridized .. %1$d\n",hcount);
	System.out.printf("Other ....... %1$d\n",ocount);

	if (zoneCount > 0) {
	    System.out.printf("Zones ....... %1$d\n", zoneCount);
	    for (int i=1; i <= zoneCount; ++i)
		System.out.printf
		    ("  Zone %1$2d: %2$2d donors  %3$2d acceptors\n",
		     i, dcounts[i], acounts[i]);
	}
    }

    private void displayAtomTypes () {
	for (int i = 0; i < atoms.length; ++i) {
	    System.err.printf("Atom %1$2d: %2$3d %3$02x\n", 
			      i+1, atoms[i].getAtno(), types[i]);
	}
	/*
	for (int i = 0; i < bonds.length; ++i) {
	    System.err.printf("Bond %1$2d: %2$d visited=%3$d (%4$d,%5$d)\n", 
			      i+1, bonds[i].getType(), bvisit[i], 
			      mol.indexOf(bonds[i].getAtom1())+1,
			      mol.indexOf(bonds[i].getAtom2())+1);
	}
	*/
    }


    protected static void outputTautomers (TautomerGenerator taugen, 
					   java.io.PrintStream os, 
					   Molecule mol) 
	throws Exception {
	String name = mol.getName();
	if (name == null || name.equals("")) {
	    name = mol.getProperty("field_0");
	    mol.setName(name);
	}

	taugen.generate(mol.cloneMolecule());
	System.out.println(mol.toFormat("smiles:q") + " " 
			   + name + " has " 
			   + taugen.getTautomerCount() 
			   + " tautomer(s)");
	/*
	  mol.setProperty
	  ("TauCount", String.valueOf(taugen.getTautomerCount()));
	  System.out.print(mol.toFormat("sdf:-a"));
	*/
	int cnt = 1;
	for (Enumeration<Molecule> tau = taugen.tautomers();
	     tau.hasMoreElements(); ++cnt) {
	    Molecule t = tau.nextElement();
	    int score = scoreTautomer (t);
	    /*
	      t.setProperty("TauScore", String.valueOf(score));
	      System.out.print(t.toFormat("sdf:-a"));
	    */
	    System.out.println(t.toFormat("smiles:u-a") + " "+cnt
			       + " " + score);
	}

	Molecule cantau = taugen.getCanonicalTautomer();
	//cantau.setName("CanonicalTautomer");
	//System.out.print(cantau.toFormat("sdf:-a"));
	System.out.println(cantau.toFormat("smiles:q") 
			   + "\tCanonical" + "\t" + scoreTautomer (cantau));

	/*
	System.out.println(MolStandardizer.canonicalSMILES(cantau) 
			   + "\tCanonical" + "\t" + scoreTautomer (cantau));
	*/
    }

    public static void main (String argv[]) throws Exception {

	TautomerGenerator taugen = new SayleDelanyTautomerGenerator 
	    (Integer.getInteger("maxtau",1001));
	//taugen.unset(FLAG_NOXIDE);
	//taugen.set(FLAG_ALL);
	//taugen.setTimeout(10000); // 10 seconds

	System.err.println("** Max tautomers generated: " 
			   + taugen.getMaxTautomers());

	if (argv.length == 0) {
	    System.err.println("** reading from STDIN...");

	    MolImporter mi = new MolImporter (System.in);
	    for (Molecule mol = new Molecule (); mi.read(mol); ) {
		outputTautomers (taugen, System.out, mol);
	    }
	}
	else {
	    for (int i = 0; i < argv.length; ++i) {
		MolImporter mi = new MolImporter (argv[i]);
		for (Molecule mol = new Molecule (); mi.read(mol); ) {
		    //st.standardize(mol);
		    outputTautomers (taugen, System.out, mol);
		}
		mi.close();
	    }
	}
    }
}
