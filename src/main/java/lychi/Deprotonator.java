package lychi;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;

import lychi.util.ChemUtil;
import lychi.util.MolPath;

/**
 * Deprotonation... please be careful making changes to this class!!!
 */
public class Deprotonator implements Comparator<MolAtom> {

    private static final Logger logger = 
	Logger.getLogger(Deprotonator.class.getName());

    static private boolean debug = false;

    static private final BitSet ATNO = new BitSet (256);
    static {
	ATNO.set(5);  //B
	ATNO.set(26); //Fe
	ATNO.set(35); //Br

	try {
	    debug = Boolean.getBoolean("deprotonator.debug");
	}
	catch (Exception ex) {
	}
    }


    static final int[] priority = new int[120];
    static {
	for (int i = 0; i < priority.length; ++i) {
	    priority[i] = Integer.MAX_VALUE;
	}
	priority[8] = 1; // O
	priority[16] = 2; // S
	priority[7] = 3; // N
    }

    static Deprotonator instance = new Deprotonator ();
    protected Deprotonator () {}
    
    public static Deprotonator getInstance () {
	return instance;
    }

    /*
     * It's assumed that the input molecule DO NOT contain explicit
     * hydrogens!!!
     */
    public void apply (Molecule mol) {
	ArrayList<MolAtom> acc = new ArrayList<MolAtom>();
	ArrayList<MolAtom> don = new ArrayList<MolAtom>();
	ArrayList<MolAtom> hac = new ArrayList<MolAtom>();

	// now see if there is any charge left
	MolAtom[] atoms = mol.getAtomArray();
	for (MolAtom a : atoms) {
	    int ch = a.getCharge();

	    if (ch <= 0) {
		if (ch == 0) {
		    int sb = 0;
		    for (int i = 0; i < a.getBondCount(); ++i) {
			int type = a.getBond(i).getType();
			if (type == 1) {
			    ++sb;
			}
		    }

		    if (a.getAtno() == 7 && sb == 3) {
			if (acc.indexOf(a) < 0) {
			    acc.add(a);
			}
		    }
		    else if (priority[a.getAtno()] < 5 
			     && a.getImplicitHcount() > 0) {
			// atom capable of accepting hydrogens
			if (hac.indexOf(a) < 0) {
			    hac.add(a);
			}
		    }
		}
		continue;
	    }

	    // see if we can neutralize this charge by donating
	    //  a proton to a neighboring atom 
	    for (int i = 0; i < a.getBondCount(); ++i) {
		MolBond b = a.getBond(i);
		MolAtom xa = b.getOtherAtom(a);
		ch += xa.getCharge();
	    }

	    if (ch == 0) { // preserve [-*]~[*+] 
		continue;
	    }

	    for (int i = 0; i < a.getBondCount(); ++i) {
		MolBond b = a.getBond(i);
		MolAtom xa = b.getOtherAtom(a);

		if (b.getType() == 2) {
		    SortedMap<MolAtom, MolBond> cand = 
			new TreeMap<MolAtom, MolBond>(this);
			
		    for (int j = 0; j < xa.getBondCount(); ++j) {
			MolBond nb = xa.getBond(j);
			if (nb != b && nb.getType() == 1) {
			    MolAtom na = nb.getOtherAtom(xa);
			    if (na.getAtno() != 6 
				&& (na.getImplicitHcount() > 0
				    || (na.getCharge() < 0
					&& !ATNO.get(na.getAtno()))
				    // or -[O+]=
				    /*|| a.getAtno() == 8*/)) {
				cand.put(na, nb);
			    }
			}
		    }
		    
		    if (!cand.isEmpty()) {
			MolAtom na = cand.firstKey();
			MolBond nb = cand.get(na);
			nb.setType(2);
			if (na.getCharge() < 0) {
			    na.setCharge(na.getCharge()+1);
			}
			a.setCharge(ch-1);
			b.setType(1);
			break; 
		    }

		    if (don.indexOf(a) < 0) {
			don.add(a);
		    }
		}
	    }
	}


	if (don.isEmpty()) {
	    return;
	}

	/*
	 * now check for exceptions...
	 */
	ArrayList<MolAtom> ignore = new ArrayList<MolAtom>();
	for (MolAtom a : don) {
	    int sb = 0, db = 0;
	    for (int i = 0; i < a.getBondCount(); ++i) {
		int t = a.getBond(i).getType();
		if (t == 1) ++sb;
		else if (t == 2) ++db;
	    }
	    if (db == 2) {
		// ignore N=[N+]=[N-] or [-*]~[*+]
		ignore.add(a);
	    }
	}

	don.removeAll(ignore);
	if (don.isEmpty()) {
	    return;
	}

	/*
	 * this attempt to do "hard" (de)protonation; e.g., methylene blue
	 */
	/*
	if (don.size() > 1) {
	    // can't handle multiple charges at the moment
	    logger.log(Level.WARNING, mol.getName() + " contains multiple "
		       +"charged groups; no attempts will be made to "
		       +"deprotonate them!");
	    return;
	}
	*/

	int[] order = ChemUtil.topologyInvariant(mol);
	for (Iterator<MolAtom> it = don.iterator(); it.hasNext(); ) {
	    MolAtom ref = it.next();
	    if (debug) {
		StringBuilder sb = new StringBuilder 
		    ("Deprotonating ref "+(mol.indexOf(ref)+1)+" => "
		     +ref.getSymbol());
		for (int b = 0; b < ref.getBondCount(); ++b) {
		    sb.append(ref.getBond(b).getOtherAtom(ref).getSymbol());
		}
		logger.info(sb.toString());
	    }
	    // proton propagation
	    if (propagate1 (ref, acc, mol, order)) {
		if (debug) {
		    logger.info("propagate1 succeeds");
		}
	    }
	    // now try hydrogen propagation
	    else if (propagate2 (ref, hac, mol, order)) {
		if (debug) {
		    logger.info("propagate2 succeeds");
		}
	    }
	}
    }

    static MolPath.Annotator debugPath = new MolPath.Annotator() {
		public String getLabel (MolAtom a) {
		    return a.getSymbol();
		}
	    public String getLabel (MolBond b) {
		int type = b.getType();
		switch (type) {
		case 1: return "";
		case 2: return "=";
		case 3: return "#";
		case MolBond.AROMATIC: return ":";
		default: return "?";
		}
	    }
	};

    static Comparator<MolBond[]> pathSorter = new Comparator<MolBond[]> () {
	public int compare (MolBond[] p1, MolBond[] p2) {
	    return p1.length - p2.length;
	}
    };

    // proton propagation
    boolean propagate1 (MolAtom ref, List<MolAtom> candidates, 
			Molecule mol, int[] topologyInvariant) {
	MolAtom[] ends = new MolAtom[2];
	MolAtom[] endsPicked = null;
	Map<MolBond, Integer> mflags = new HashMap<MolBond, Integer>();

	MolPath mp = new MolPath (mol);

	int i = mol.indexOf(ref);
	for (MolAtom atom : candidates) {
	    int j = mol.indexOf(atom);

	    MolBond[][] paths = mp.getPaths(i, j);
	    Arrays.sort(paths, pathSorter);

	    //System.out.println("Paths between "+(i+1) +" "+(j+1));
	    for (MolBond[] path : paths) {
		int[] flags = new int[path.length];

		boolean ok = propagateCharge 
		    (mol, topologyInvariant, path, flags, ends);

		if (ok) {
		    if (mflags.isEmpty()) {
			endsPicked = new MolAtom[2];
			endsPicked[0]  = ends[0];
			endsPicked[1]  = ends[1];

			for (int k = 0; k < path.length; ++k) {
			    mflags.put(path[k], flags[k]);
			}
		    }
		    else {
			for (int k = 0; k < path.length; ++k) {
			    if (!mflags.containsKey(path[k])) {
				// preserve the bond type for paths that 
				//   weren't choosen
				mflags.put(path[k], 
					   (flags[k] & ~MolBond.TYPE_MASK) 
					   | path[k].getType());
			    }
			}
		    }
		}
	    }
	}

	boolean ok = !mflags.isEmpty();
	if (ok) {
	    // now propagate the bond flags
	    for (Map.Entry<MolBond, Integer> e : mflags.entrySet()) {
		MolBond b = e.getKey();
		b.setFlags(e.getValue());
	    }

	    MolAtom head = endsPicked[0];
	    MolAtom tail = endsPicked[1];
	    if (head.getAtno() == tail.getAtno()) {
		int hi = mol.indexOf(head);
		int ti = mol.indexOf(tail);
		if (topologyInvariant[hi] <= topologyInvariant[ti]) {
		}
		else {
		    // movements of proton from head to tail
		    exchangeCharge (head, tail);
		}
	    }
	    else {
		exchangeCharge (head, tail);
	    }
	}

	return ok;
    }

    void exchangeCharge (MolAtom from, MolAtom to) {
	to.setCharge(from.getCharge());
	from.setCharge(0);
	
	// example: C1=CC=C2C(=C1)C3=[N+](C4=CC=CC=C4N(C3=N2)O)[O-]
	// we also check to see if there is any negative charge
	//  that need propagation too; i.e., [-*]~[+*]
	List<MolAtom> pairs = new ArrayList<MolAtom>();
	for (int i = 0; i < from.getBondCount(); ++i) {
	    MolAtom a = from.getBond(i).getOtherAtom(from);
	    if (a.getCharge() != 0 && a != to) { 
		pairs.add(a);
	    }
	}
	
	// now check to see if any of the neighbors of tail
	//  is the same atom
	if (!pairs.isEmpty()) {
	    List<MolAtom> neutral = new ArrayList<MolAtom>();

	    for (MolAtom xa : pairs) {
		for (int i = 0; i < to.getBondCount(); ++i) {
		    MolAtom a = to.getBond(i).getOtherAtom(to);
		    if (a.getAtno() == xa.getAtno() 
			&& a.getCharge() != xa.getCharge()) {
			a.setCharge(xa.getCharge());
			neutral.add(xa);
			break;
		    }
		}
	    }

	    for (MolAtom a : neutral) {
		a.setCharge(0);
	    }
	}
    }

    // hydrogen propagation
    boolean propagate2 (MolAtom ref, List<MolAtom> candidates, 
			Molecule mol, int[] topologyInvariant) {
	MolAtom[] ends = new MolAtom[2];
	MolAtom[] endsPicked = null;
	Map<MolBond, Integer> mflags = new HashMap<MolBond, Integer>();

	MolPath mp = new MolPath (mol);
	int i = mol.indexOf(ref);

	Collections.sort(candidates, this);
	for (final MolAtom atom : candidates) {
	    int j = mol.indexOf(atom);
	    if (debug) {
		StringBuilder sb = new StringBuilder 
		    ("Candicate "+(j+1)+" => "+atom.getSymbol());
		for (int k = 0; k < atom.getBondCount(); ++k) {
		    sb.append(atom.getBond(k).getOtherAtom(atom).getSymbol());
		}
		logger.info(sb.toString());
	    }

	    MolBond[][] paths = mp.getPaths(i, j);
	    Arrays.sort(paths, pathSorter);

	    for (MolBond[] path : paths) {
		int[] flags = new int[path.length];

		if (debug) {
		    logger.info("Path "+MolPath.toString(path, debugPath));
		}
		boolean ok = propagateCharge 
		    (mol, topologyInvariant, path, flags, ends);
		if (debug) {
		    StringBuilder sb = new StringBuilder
			("..ok? "+ ok + " => ");
		    for (int k = 0; k < flags.length; ++k) {
			sb.append(String.valueOf(flags[k]&MolBond.TYPE_MASK));
		    }
		    logger.info(sb.toString());
		}


		if (ok) {
		    if (mflags.isEmpty()) {
			endsPicked = new MolAtom[2];
			endsPicked[0] = ends[0];
			endsPicked[1] = ends[1];

			for (int k = 0; k < path.length; ++k) {
			    mflags.put(path[k], flags[k]);
			}
		    }
		    else {
			for (int k = 0; k < path.length; ++k) {
			    if (!mflags.containsKey(path[k])) {
				mflags.put(path[k], 
					   (flags[k] & ~MolBond.TYPE_MASK) 
					   | path[k].getType());
			    }
			}
		    }
		}
	    }
	}

	boolean ok = !mflags.isEmpty();
	if (ok) {
	    // only movement of hydrogens
	    MolAtom head = endsPicked[0];
	    MolAtom tail = endsPicked[1];

	    // now propagate the bond flags
	    for (Map.Entry<MolBond, Integer> e : mflags.entrySet()) {
		MolBond b = e.getKey();
		b.setFlags(e.getValue());
	    }

	    if (head.getAtno() == tail.getAtno()) {
		int hi = mol.indexOf(head);
		int ti = mol.indexOf(tail);
		if (topologyInvariant[hi] <= topologyInvariant[ti]) {
		}
		else {
		    head.setCharge(0);
		}
	    }
	    else {
		head.setCharge(0);
	    }
	}

	return ok;
    }

    boolean propagateCharge (Molecule mol, int[] order, 
			     MolBond[] path, int[] flags, 
			     MolAtom[] ends) {
	if (path.length % 2 != 0) {
	    return false;
	}

	if (path.length != flags.length) {
	    throw new IllegalArgumentException ("Flags not same size as path");
	}

	if (ends.length < 2) {
	    throw new IllegalArgumentException
		("Not sufficient length for head/tail");
	}
	
	int type = path[0].getType();
	if (type != 2) {
	    return false;
	}

	for (int i = 0; i < path.length; ++i) {
	    if (path[i].getType() != type) {
		// not alternating double-single 
		return false;
	    }
	    type = 3 - type;
	}

	MolAtom head = path[0].getAtom1();
	if (head.getCharge() == 0) {
	    head = path[0].getAtom2();
	}

	if (head.getCharge() == 0) {
	    // bogus path that doesn't have an atom of type =[A+]
	    return false;
	}

	int[] types = new int[path.length];

	MolAtom tail = head;
	int ch = head.getCharge();

	MolAtom[] apath = new MolAtom[path.length+1];
	// original valence of atoms on this path
	int[] hcount = new int[path.length+1]; 

	apath[0] = head;
	hcount[0] = head.getImplicitHcount();
	type = path[0].getType();
	for (int i = 0; i < path.length; ++i) {
	    tail = path[i].getOtherAtom(tail);
	    if (tail == null) {
		logger.log(Level.WARNING, "Disconnected path");
		return false;
	    }
	    apath[i+1] = tail;
	    hcount[i+1] = tail.getImplicitHcount();

	    types[i] = path[i].getType();
	    flags[i] = path[i].getFlags();

	    type = 3 - type;
	    path[i].setFlags(type, MolBond.TYPE_MASK);
	    //path[i].setType(type);
	}

	ends[0] = head;
	ends[1] = tail;

	head.setCharge(0);
	if (tail.getImplicitHcount() == 0) {
	    tail.setCharge(ch);
	}

	// now check valence error
	for (int i = 0; i < apath.length; ++i) {
	    if (!checkValence (apath[i], hcount[i])) {
		// undo...
		for (int j = 0; j < path.length; ++j) {
		    path[j].setFlags(types[j], MolBond.TYPE_MASK);
		    //path[j].setType(types[j]);
		}
		head.setCharge(ch);
		tail.setCharge(0);
		//System.err.println("==> "+mol.toFormat("smiles:q"));

		return false;
	    }
	}

	// all is good
	if (debug) {
	    StringBuilder sb = new StringBuilder ("success propagation: ");
	    for (int i = 0; i < path.length; ++i) {
		sb.append(path[i].getFlags() & MolBond.TYPE_MASK);
	    }
	    logger.info(sb.toString());
	}

	// make sure stereo flags are 0
	for (int i = 0; i < path.length; ++i) {
	    /*
	      int stereo = b.getFlags() & MolBond.STEREO_MASK;
	      if (stereo != 0) {
	      logger.info("Ignore stereo ("+stereo+") for bond "
	      +(mol.indexOf(b) +1) + " due to deprotonation");
	      }
	    */
	    type = path[i].getFlags() & MolBond.TYPE_MASK;
	    // restore the bond type
	    path[i].setFlags(types[i], MolBond.TYPE_MASK); 

	    // update the flags with the new type and turn off stereo mask
	    flags[i] = ((flags[i] & ~MolBond.TYPE_MASK) | type)
		& ~MolBond.STEREO_MASK;
	}

	if (debug) {
	    StringBuilder sb = new StringBuilder ("new flags: ");
	    for (int i = 0; i < flags.length; ++i) {
		sb.append(String.valueOf(flags[i] & MolBond.TYPE_MASK));
	    }
	    logger.info(sb.toString());
	}

	// restore charge 
	head.setCharge(ch);
	tail.setCharge(0);

	//logger.info(mol.toFormat("smiles:q"));

	// ok, now check to see if the propagation is in 
	//  the right direction.  if not, we revert the
	//  changes, but still keep the stereo flags reset
	//  if any.
	if (head.getAtno() == tail.getAtno()) {
	    int i = mol.indexOf(head);
	    int j = mol.indexOf(tail);
	    if (order[i] <= order[j]) {
		// revert the changes... but still keep the changes to 
		//   stereo flag
		for (int k = 0; k < path.length; ++k) {
		    //path[k].setFlags(types[k], MolBond.TYPE_MASK);
		    flags[k] = (flags[k] & ~MolBond.TYPE_MASK) | types[k];
		}

		if (debug) {
		    StringBuilder sb = new StringBuilder ("restored flags: ");
		    for (int k = 0; k < flags.length; ++k) {
			sb.append(String.valueOf
				  (flags[k] & MolBond.TYPE_MASK));
		    }
		    logger.info(sb.toString());
		}
	    }
	}
	
	return true;
    }

    static boolean checkValence (MolAtom a, int hcount) {
	a.valenceCheck();

	if (a.hasValenceError()) {
	    return false;
	}

	/*
	System.out.println("["+a.getSymbol()+",h="
			   +a.getImplicitHcount()
			   +",v="+a.getValence()
			   +",H="+hcount+"]");
	*/

	// only allow alternating path if the valences of atoms
	//   are satisfied
	int h = a.getImplicitHcount();
	int atno = a.getAtno();

	if (atno == 6) {
	    if (h != hcount) {
		return false;
	    }
	}
	else if (h <= hcount) { // [O,S,N]
	    // we're good
	}
	else {
	    return false;
	}

	return true;
    }

    public int compare (MolAtom a, MolAtom b) {
	int d = priority[a.getAtno()] -priority[b.getAtno()];
	if (d == 0) {
	    // have a preference for the atom with the high (implicit) H
	    //   count
	    d = b.getImplicitHcount() - a.getImplicitHcount();
	}
	return d;
    }

    public static void main (String[] argv) throws Exception {
	Deprotonator dp = getInstance ();
	for (String s : argv) {
	    MolHandler mh = new MolHandler (s);
	    Molecule mol = mh.getMolecule();
	    mol.aromatize(false);
	    dp.apply(mol);
	    System.out.println(s + " => " + mol.toFormat("smiles:q"));
	}
	if (argv.length == 0) {
	    String s= "C[N+](C)=C(\\C=C/N)/C=C\\O";
	    MolHandler mh = new MolHandler (s);
	    Molecule mol = mh.getMolecule();
	    mol.aromatize(false);
	    dp.apply(mol);
	    System.out.println(s + " => " + mol.toFormat("smiles:q"));
	}
    }
}
