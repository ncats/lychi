package lychi.util;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import chemaxon.struc.MolBond;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import chemaxon.formats.MolImporter;


public class ChemUtil {
    static private Logger logger = Logger.getLogger(ChemUtil.class.getName());
    static boolean DEBUG = false;
    static {
	try {
	    DEBUG = Boolean.getBoolean("lychi-chemutil.debug");
	}
	catch (Exception ex) {
	}
    }

    /**
     * Any change to calcMolScore must update this value.  The 
     * MolStandardizer class is very sensitive to how calMolScore
     * behaves!
     */
    public static final int MOLSCORE_VERSION = 0x1;


    static class IntPair implements Comparable<IntPair> {
	public int ival1, ival2;
	public IntPair (int ival1, int ival2) {
	    this.ival1 = ival1;
	    this.ival2 = ival2;
	}
	public int compareTo (IntPair ip) {
	    return ival2 - ip.ival2;
	}
    }

    public static int complexity (Molecule mol) {
	int index = 0;
	
	if (mol != null) {
	    int[][] sssr = mol.getSSSR();
	    for (int i = 0; i < sssr.length; ++i) {
		index += sssr[i].length * 6;
	    }
	    
	    for (int i = 0; i < mol.getAtomCount(); ++i) {
		MolAtom a = mol.getAtom(i);
		int nb = a.getBondCount();
		switch (nb) {
		case 4: index += 24; break;
		case 3: index += 12; break;
		case 2: index += 6; break;
		case 1: index += 3; break;
		}
		
		index += a.getAtno() == 6 ? 3 : 6;
	    }
	}
	return index;
    }

    public static int[] topologyInvariant (Molecule mol) {
	Molecule m = mol.cloneMolecule();
	m.hydrogenize(false);
	m.expandSgroups();

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

	int[] map = new int[m.getAtomCount()];
	m.getGrinv(map);

	return map;
    }

    public static int[] graphInvariantOrder (Molecule mol) {
	int[] gi = new int[mol.getAtomCount()];
	mol.getGrinv(gi);
	IntPair[] pairs = new IntPair[gi.length];
	for (int i = 0; i < pairs.length; ++i) {
	    pairs[i] = new IntPair (i, gi[i]);
	}
	Arrays.sort(pairs);
	for (int i = 0; i < pairs.length; ++i) {
	    gi[i] = pairs[i].ival1;
	}
	pairs = null;
	return gi;
    }

    public static int[] graphInvariantOrder (Molecule mol, final int[] rank) {
	if (rank.length != mol.getAtomCount()) {
	    throw new IllegalArgumentException
		("Input rank doesn't match number of atoms in molecule");
	}

	int[] gi = new int[mol.getAtomCount()];
	mol.getGrinv(gi, Molecule.GRINV_NOHYDROGEN);

	IntPair[] pairs = new IntPair[gi.length];
	for (int i = 0; i < pairs.length; ++i) {
	    pairs[i] = new IntPair (i, gi[i]);
	}

	Arrays.sort(pairs, new Comparator<IntPair>() {
			public int compare (IntPair ip1, IntPair ip2) {
			    int d = ip1.ival2 - ip2.ival2;
			    if (d == 0) {
				d = rank[ip1.ival1] - rank[ip2.ival1];
			    }
			    return d;
			}
		    });

	for (int i = 0; i < pairs.length; ++i) {
	    //System.out.println(i + ": " + pairs[i].ival1 + " " + pairs[i].ival2 + " " + rank[pairs[i].ival1]);
	    gi[i] = pairs[i].ival1;
	}
	pairs = null;
	return gi;
    }


    public static int[] graphInvariantOrder (Molecule mol, Molecule ref) {
	if (mol.getAtomCount() != ref.getAtomCount()) {
	    throw new IllegalArgumentException 
		("Input and ref molecules don't match");
	}
	int[] rank = new int[ref.getAtomCount()];
	for (int i = 0; i < rank.length; ++i) {
	    rank[i] = ref.getAtom(i).getAtno();
	}
	return graphInvariantOrder (mol, rank);
    }


    /*
     * Unlike SSSR, this function considers fused and spiro rings
     *  as part of the same ring system.
     */
    public static int[][] getRingSystems (Molecule m) {
	return getRingSystems (m.getSSSR());
    }

    public static int[][] getRingSystems (int[][] sssr) {

	BitSet[] rings = new BitSet[sssr.length];
	for (int i = 0; i < sssr.length; ++i) {
	    BitSet ri = new BitSet (sssr[i].length);
	    for (int k = 0; k < sssr[i].length; ++k) {
		ri.set(sssr[i][k]);
	    }
	    rings[i] = ri;
	}

	UnionFind eqv = new UnionFind (sssr.length);
	for (int i = 0; i < rings.length; ++i) {
	    for (int j = i+1; j < rings.length; ++j) {
		// see if rings ri & rj share common atoms
		if (rings[i].intersects(rings[j])) {
		    eqv.union(i, j);
		}
	    }
	}
	rings = null;

	// equivalence classes
	int[][] eqc = eqv.getComponents();
	int[][] ringSystems = new int[eqc.length][];
	BitSet bs = new BitSet ();
	for (int i = 0; i < eqc.length; ++i) {
	    bs.clear();
	    for (int j = 0; j < eqc[i].length; ++j) {
		int[] r = sssr[eqc[i][j]];
		for (int n = 0; n < r.length; ++n)
		    bs.set(r[n]);
	    }
	    int[] ring = new int[bs.cardinality()];
	    for (int j = bs.nextSetBit(0), k = 0; 
		 j >= 0; j = bs.nextSetBit(j+1)) 
		ring[k++] = j;
	    ringSystems[i] = ring;
	}
	sssr = null;

	Arrays.sort(ringSystems, new Comparator<int[]>() {
			public int compare (int[] r1, int[] r2) {
			    int d = r2.length - r1.length;
			    for (int i = 0; d==0 && i < r1.length; ++i) {
				d = r1[i] - r2[i];
			    }
			    return d;
			}
		    });
	
	return ringSystems;
    }


    // the scoring function below is a modified version of the one 
    // described by
    //   Oellien et al. j. chem inf model 2006, 46, 2342-54
    public static int calcMolScore (Molecule m) {
	int score = 0;

	Molecule mol = m.cloneMolecule();
	//mol.aromatize(Molecule.AROM_BASIC);
	mol.aromatize();

	int[][] sssr = mol.getSSSR();
	int[] rsizes = mol.getSmallestRingSizeForIdx();

	/* don't use this method... very expensive!
	   int[][][] rings = tau.getAromaticAndAliphaticRings
	   (Molecule.AROM_BASIC, true, false, 0, 0);
	   score += rings[0].length *100;
	*/

	MolBond[] bonds = mol.getBondArray();
	for (int i = 0; i < bonds.length; ++i) {
	    MolBond b = bonds[i];
	    int a1 = b.getAtom1().getAtno();
	    int a2 = b.getAtom2().getAtno();

	    if (b.getType() == 2) {
		// each double bond between carbon and heteroatoms +1
		if ((a1 == 6 && a2 != 6) || (a1 != 6 && a2 == 6)) {
		    score += 1;
		}
		// discourage carbon = carbon
		else if (a1 == 6 && a2 == 6 && !mol.isRingBond(i)) {
		    score -= 55;
		}
		// each double bond between oxygen and nitrogen +2
		else if ((a1 == 7 && a2 == 8) || (a1 == 8 && a2 == 7)) {
		    score += 2;
		}

		if (b.getAtom1().isTerminalAtom()
		    || b.getAtom2().isTerminalAtom()) {
		    score += 1;
		}

		if (!mol.isRingBond(i)) {
		    int s = b.calcStereo2();
		    switch (s) {
		    case MolBond.CIS:
		    case MolBond.TRANS:
		    case MolBond.CTUNSPEC:
			// penalize extra e/z bonds
			score -= 55;
			break;
		    }
		}

		// discourage (non-terminal) =N coming off rings
		if (((a1 == 7 && rsizes[mol.indexOf(b.getAtom2())] > 0) 
		     || (a2 == 7 && rsizes[mol.indexOf(b.getAtom1())] > 0))
		    /*&& ((a1 == 7 && !b.getAtom1().isTerminalAtom())
		      || (a2 == 7 && !b.getAtom2().isTerminalAtom()))*/) {
		    score -= 4;
		}

		// each C=N[OH] fragment +4
		if ((a1 == 6 && a2 == 7) || (a1 == 7 && a2 == 6)) {
		    MolAtom xa = a1 == 6 && a2 == 7 
			? b.getAtom2() : b.getAtom1();
		    int nb = xa.getBondCount();
		    if (nb == 2) {
			for (int j = 0; j < nb; ++j) {
			    MolBond xb = xa.getBond(j);
			    if (xb.getType() == 1 
				&& xb.getOtherAtom(xa).getAtno() == 8) {
				score += 4;
			    }
			}
		    }
		}

		// each double bond between carbon and oxygen +2
		if ((a1 == 6 && a2 == 8) || (a1 == 8 && a2 == 6)) {
		    score += 4;
		    if (rsizes[mol.indexOf(b.getAtom1())] > 0
			|| rsizes[mol.indexOf(b.getAtom2())] > 0) {
			// favor =O coming off of a ring
			score += 48; 
		    }
		}

		if (a1 == 8 || a2 == 8) {
		    score += 2;
		}

		if (mol.isRingBond(i)) {
		    score += 1;
		}
	    }
	    else if (b.getType() == 1) {
		// each P-H, S-H, Se-H, and Te-H bond -1
		if (b.getAtom1().isTerminalAtom()
		    && (a1 == 15 || a1 == 16 || a1 == 34 || a1 == 52)) {
		    score -= 1;
		}
		else if (b.getAtom2().isTerminalAtom()
			 && (a2 == 15 || a2 == 16 || a2 == 34 || a2 == 52)) {
		    score -= 1;
		}
	    }
	}
	bonds = null;



	if (DEBUG) {
	    logger.info("score 1: "+score);
	}

	//mol.hydrogenize(false);

	// each methyl group (applying a penalty to structures with 
	//   terminal double bonds) +1
	MolAtom[] atoms = mol.getAtomArray(); 
	int chargeCount = 0;
	for (MolAtom a : atoms) {
	    // there shouldn't be any explicit H, but what the hell...
	    if ((a.getImplicitHcount() + a.getExplicitHcount()) == 3) {
		score += 1;
	    }

	    // not in the paper
	    int c = 2*a.getCharge(); // scale the charge values
	    if (c != 0) {
		switch (a.getAtno()) {
		case 7: score += c; break;
		case 8: score -= c; break;
		default: // discourage charged ions
		    score -= Math.abs(c); 
		    break;
		}
		++chargeCount;
	    }
	}
	// penalize molecule with too many charged atoms
	score -= chargeCount*10;

	if (DEBUG) {
	    logger.info("score 2: "+score);
	}


	for (int i = 0; i < sssr.length; ++i) {
	    int naro = 0;
	    for (int j = 0; j < sssr[i].length; ++j) {
		if (atoms[sssr[i][j]].hasAromaticBond()) {
		    ++naro;
		}
	    }

	    if (naro == sssr[i].length) {
		//each aromatic ring system +100
		score += 100;
	    }
	}
	atoms = null;

	if (DEBUG) {
	    logger.info("score 3: "+score);
	}

	/*
	 * we have to do this hack to get the alternating single-double
	 * bond to be in a consistent layout for an aromatic ring
	 */
	try {
	    MolHandler mh = new MolHandler (mol.toFormat("smiles:q"));
	    m = mh.getMolecule();
	    m.dearomatize();
	    for (MolAtom a : m.getAtomArray()) {
		if (a.hasAromaticBond()) {
		    throw new Exception ("bailing out");
		}
	    }
	    mol = m;
	    sssr = mol.getSSSR();
	    rsizes = mol.getSmallestRingSizeForIdx();
	}
	catch (Exception ex) {
	    mol.dearomatize();
	}
	
	// prefer double-bond shared between rings
	Map<Integer, Integer> ringCount = new HashMap<Integer, Integer>();
	
	for (int i = 0; i < sssr.length; ++i) {
	    for (int j = 0; j < sssr[i].length; ++j) {
		Integer c = ringCount.get(sssr[i][j]);
		ringCount.put(sssr[i][j], c == null ? 1 : (c+1));
	    }
	}
	
	for (MolBond b : mol.getBondArray()) {
	    int a1 = mol.indexOf(b.getAtom1());
	    int a2 = mol.indexOf(b.getAtom2());
	    
	    if (b.getType() == 2 && rsizes[a1] > 0 && rsizes[a2] > 0) {
		if (ringCount.get(a1) > 1 && ringCount.get(a2) > 1) {
		    score += 1;
		}
	    }
	}
	ringCount.clear();
	sssr = null;

	if (DEBUG) {
	    logger.info("score 4: "+score);
	}

	return score;
    }

    static public String generateSkeleton (Molecule m) {
	List<MolAtom> remove = new ArrayList<MolAtom>();
	for (MolAtom a : m.getAtomArray()) {
	    if (a.isTerminalAtom()) {
		// remove all terminal atoms
		remove.add(a);
	    }
	    else {
		a.setRadical(0);
		a.setCharge(0);
		a.setFlags(0);
		a.setAtno(6);
	    }
	}
	for (MolBond b : m.getBondArray()) {
	    b.setFlags(0);
	    b.setType(1);
	}

	for (MolAtom a : remove) {
	    m.removeNode(a);
	}

	return m.toFormat("cxsmarts");
    }

    public static void resetEZ (Molecule mol, int flags) {
	MolBond[] bonds = mol.getBondArray();
	for (int i = 0; i < bonds.length; ++i) {
	    MolBond b = bonds[i];
	    if (b.getType() == 2) {
		MolAtom a1 = b.getCTAtom1();
		MolAtom a4 = b.getCTAtom4();
		if (a1 != null && a4 != null) {
		    b.setStereo2Flags(a1, a4, flags);
		}
	    }
	}
    }

    public static void resetEZ (Molecule mol) {
	MolBond[] bonds = mol.getBondArray();
	for (int i = 0; i < bonds.length; ++i) {
	    MolBond b = bonds[i];
	    int type = b.getType();
	    if (type == 2) {
		b.setFlags(0, MolBond.CTUMASK);
	    }
	    else if (type == 1) {
		int flags = b.getFlags() & MolBond.STEREO1_MASK;
		if (flags  == (MolBond.UP | MolBond.DOWN)) {
		    // WAVY flag... 
		    b.setFlags(0, MolBond.STEREO1_MASK);
		}
	    }
	}
    }

    public static void adjustEZStereo (Molecule mol) {
	for (MolBond b : mol.getBondArray()) {
	    if (b.getType() == 2) {
		MolAtom a1 = b.getCTAtom1();
		MolAtom a4 = b.getCTAtom4();
		if (a1 != null && a4 != null 
		    && (a1.getX() != 0. || a1.getY() != 0.)
		    && (a4.getX() != 0. || a4.getY() != 0.)) {
		    b.setStereo2Flags(a1, a4, b.calcStereo2(a1, a4));
		}
	    }
	}
    }

    public static String canonicalSMILES (Molecule m) {
	return canonicalSMILES (m, true);
    }

    public static String canonicalSMILES (Molecule m, boolean stereo) {
	return canonicalSMILES (null, m, stereo);
    }

    public static String canonicalSMILES 
	(Molecule out, Molecule mol, boolean stereo) {
	String smiles = null;
	
	Molecule m = mol.cloneMolecule();
	MolAtom[] atoms = m.getAtomArray();
	int[] amap = new int[atoms.length];
	for (int i = 0; i < atoms.length; ++i) {
	    int j = atoms[i].getAtomMap();
	    if (j > 0) {
		amap[i] = j;
		atoms[i].setAtomMap(0);
	    }
	    else {
		amap[i] = -(i+1); 
	    }
	}

	try {
	    /*
	     * here we force the dimension to be zero to turn
	     * off automatic stereo perceptions when the molecule is
	     * being converted into smiles.  this means that the 
	     * molecule's stereo flags must be set appropriately!
	     */
	    //m.setDim(0); // why was this set before????
	    smiles = m.toFormat("smiles:q" + (stereo?"":"0"));
	}
	catch (Exception ex) {
	    try {
		smiles = m.toFormat("cxsmarts:u");

		String[] toks = smiles.split("[\\s]");
		if (toks.length > 0) {
		    // trim off any extra info and treat it as a normal
		    //  smiles/smarts
		    smiles = toks[0];
		}
	    }
	    catch (Exception e) {
		logger.log(Level.SEVERE, "Can't convert input molecule "
			   +"into smiles or cxsmarts", e);
		System.err.println(m.toFormat("sdf"));
		smiles = null;
	    }
	}

	if (smiles != null) {
	    StringBuilder sb = new StringBuilder ();
	    int br = 0;
	    for (int i = 0; i < smiles.length(); ++i) {
		char ch = smiles.charAt(i);
		if (ch == '[' || ch == ']') {
		    br ^= 1;
		}
		else if (ch == '-' && br == 0) {
		    // skip single bond
		    continue;
		}
		sb.append(ch);
	    }
	    smiles = sb.toString();
	}

	/**
	 * Oh, why there isn't a method that allows one to iterate
	 * over the atoms in a canonical order independent of the 
	 * underlying order????  This is really an ugly hack to do 
	 * just that!!!!  The rest of the code simply transform the
	 * input molecule (mol) into a canonical order (out).  The
	 * output molecule is allowed to be the same instance as the
	 * input molecule.
	 */	
	if (out != null && smiles != null) {
	    // save these in case mol == out
	    String name = mol.getName();
	    int dim = mol.getDim();
	    boolean absStereo = mol.isAbsStereo();

	    // copy properties
	    Map<String, Object> props = new HashMap<String, Object>();
	    for (int i = 0; i < mol.getPropertyCount(); ++i) {
		String prop = mol.getPropertyKey(i);
		props.put(prop, mol.getPropertyObject(prop));
	    }

	    int[] molGi = new int[atoms.length];
	    mol.getGrinv(molGi);

	    try {
		out.clear();
		for (int i = 0; i < atoms.length; ++i) {
		    atoms[i].setAtomMap(amap[i] < 0 ? -amap[i] : amap[i]);
		}

		// generate a version with the original atom mapping
		String s = m.toFormat("smiles:q"+(stereo?"":"0"));
		MolImporter.importMol(s, out);
		if (DEBUG) {
		    logger.info(smiles + " <==> " + out.toFormat("smiles:q"));
		}

		/*
		 * For whatever reason, we can call out.setDim() here, doing
		 * so will reset any stereo flags!! So we have to resort
		 * to this hack... sign.
		 */
		if (false && dim > 0) {
		    int[] fixed = new int[out.getAtomCount()];
		    for (int i = 0; i < fixed.length; ++i) {
			fixed[i] = i;
		    }
		    out.partialClean(dim, fixed, null);
		}
		out.setName(name);
		out.setAbsStereo(absStereo);
		if (DEBUG) {
		    logger.info("Restored 0: " + out.toFormat("smiles:q"));
		}

		// now preserve the coordinate
		int unmapped = 0;
		MolAtom[] outAtoms = out.getAtomArray();
		for (MolAtom a : outAtoms) {
		    int j = a.getAtomMap() - 1;
		    if (j >= 0 && j < atoms.length) {
			//a.setAtomMap(amap[j] > 0 ? amap[j] : 0);
			a.setX(atoms[j].getX());
			a.setY(atoms[j].getY());
			a.setZ(atoms[j].getZ());
			a.setFlags(atoms[j].getFlags());
		    }
		    else {
			a.setAtomMap(1023);
			++unmapped;
		    }
		}

		if (DEBUG) {
		    logger.info("Restored 1: " + out.toFormat("smiles:q"));
		}

		if (unmapped > 0) {
		    // there seem to be a bug when parsing atom maps 
		    //   in smiles when the atom has stereo, e.g., 
		    //   ...[S@:4]... 
		    if (unmapped > 1) {
			logger.log(Level.WARNING, name+": "+unmapped
				   +" unmapped atom(s) in "+s);
		    }

		    int[] outGi = new int[atoms.length];
		    out.getGrinv(outGi);

		    for (int j = 0; j < outGi.length; ++j) {
			MolAtom a = out.getAtom(j);
			if (a.getAtomMap() == 1023) {
			    for (int i = 0; i < molGi.length; ++i) {
				if (outGi[j] == molGi[i]) {
				    a.setAtomMap(i+1);
				    a.setX(atoms[i].getX());
				    a.setY(atoms[i].getY());
				    a.setZ(atoms[i].getZ());
				    a.setFlags(atoms[i].getFlags());
				    if (DEBUG) {
					logger.info(name+": Atom "
						    +(j+1)+" maps to "+(i+1));
				    }

				    --unmapped;
				    break;
				}
			    }
			}
		    }

		    if (unmapped > 0) {
			// hopefully we never get here.... but who knows
			for (int i = 0; i < out.getAtomCount(); ++i) {
			    MolAtom a = out.getAtom(i);
			    if (a.getAtomMap() == 1023) {
				logger.log(Level.WARNING, name+": Atom "+(i+1)+
					   " in output is unmapped, so its "
					   +"stereo flags and coordinates "
					   +"might be wrong!");
				a.setAtomMap(0);
			    }
			}
		    }
		}

		// perserve the bond flags
		{
		    int[][] atab = out.getBtab();
		    int[][] btab = m.getBtab();
		    for (int i = 0; i < outAtoms.length; ++i) {
			for (int j = 0; j < outAtoms.length; ++j) {
			    if (atab[i][j] >= 0) {
				int mi = outAtoms[i].getAtomMap()-1;
				int mj = outAtoms[j].getAtomMap()-1;
				if (mi >= 0 && mj >= 0 
				    && mi < atoms.length
				    && mj < atoms.length 
				    && btab[mi][mj] >= 0) {
				    MolBond b = out.getBond(atab[i][j]);
				    b.setFlags
					(m.getBond(btab[mi][mj]).getFlags());
				}
			    }
			}
		    }
		}

		// reset atom mapping 
		for (MolAtom a : outAtoms) {
		    int i = a.getAtomMap() - 1;
		    if (i >= 0 && i < amap.length) {
			// restore the atom mapping if any, else clear out
			a.setAtomMap(amap[i] > 0 ? amap[i] : 0);
		    }
		}

		//logger.info("#2: " + out.toFormat("smiles:q"));
		
		// only recalculate the 2-d stereo if we have 
		//   proper coordinates
		adjustEZStereo (out);

		for (Map.Entry<String, Object> e : props.entrySet()) {
		    out.setPropertyObject(e.getKey(), e.getValue());
		}
	    }
	    catch (Exception ex) {
		logger.log(Level.SEVERE, 
			   "Can't import canonical smiles: "+smiles
			   +" "+mol.toFormat("cxsmiles:q"), ex);
		if (out != mol) {
		    mol.clonecopy(out);
		}
	    }
	}

	return smiles;
    }

    static Pattern ISOTOPE_PATTERN = Pattern.compile
	("\\[(\\d+)([A-Za-z&&[^hH]])([^]]*)\\]");

    public static String clearIsotopes (String smiles) {
	Matcher mp = ISOTOPE_PATTERN.matcher(smiles);
	StringBuffer sb = new StringBuffer ();
	int start = 0;
	while (mp.find()) {
	    /*
	      for (int i = 0; i <= mp.groupCount(); ++i) {
	      System.out.println(" group " + i + " [" + mp.start(i)
	      + "," + mp.end(i) + "): " 
	      + mp.group(i));
	      }
	    */

	    sb.append(smiles.substring(start, mp.start(0)));
	    // skip isotope value: 
	    //   [15n] => n
	    //   [14cH] => [cH]
	    if (mp.end(3) > mp.start(3)) {
		sb.append("[");
		sb.append(smiles.substring(mp.start(2), mp.end(3)));
		sb.append("]");
	    }
	    else {
		sb.append(smiles.substring(mp.start(2), mp.end(3)));
	    }
	    start = mp.end(0);
	}
	sb.append(smiles.substring(start));

	if (start > 0) {
	    try {
		Molecule m = MolImporter.importMol(sb.toString());
		m.dearomatize();
		smiles = canonicalSMILES (m);
	    }
	    catch (Exception ex) {
		System.err.println("** warning: clearIsotopes() generates "
				   +"bogus SMILES: " + sb 
				   + "; revert changes!");
	    }
	}

	return smiles;
    }

    private ChemUtil () {
    }
}
