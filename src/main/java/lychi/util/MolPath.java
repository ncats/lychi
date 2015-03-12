package lychi.util;

import java.io.InputStream;
import java.io.FileInputStream;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.formats.MolImporter;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import chemaxon.sss.search.MolSearch;

public class MolPath {
    private static final Logger logger = 
	Logger.getLogger(MolPath.class.getName());

    // maximum recursion depth for DFS
    static final int MAX_DEPTH = 32;

    public interface PathVisitor {
	void visit (Molecule mol, List<MolBond> path);
    }

    public interface Annotator {
	String getLabel (MolAtom a);
	String getLabel (MolBond b);
    }

    public static class ShortestPaths extends MolPath 
	implements PathVisitor, Comparator<String> {

	List<MolBond[]> paths = new ArrayList<MolBond[]>();
	int[][] cost;
	int[][] btab;

	public ShortestPaths () {}
	public ShortestPaths (Molecule mol) {
	    super (mol);
	}
	
	public MolBond[][] getPaths (MolAtom a, MolAtom b) {
	    Molecule m = getMolecule ();
	    return getPaths (m.indexOf(a), m.indexOf(b));
	}

	@Override
	public void setMolecule (Molecule mol) {
	    super.setMolecule(mol);

	    btab = mol.getBtab();

	    // init the path cost...
	    int acount = getAtomCount ();
	    cost = new int[acount][acount];

	    for (int i = 0; i < acount; ++i) {
		cost[i][i] = 0;
		for (int j = i+1; j < acount; ++j) {
		    cost[i][j] = cost[j][i] = 
			btab[i][j] < 0 ? acount : 1;
		}
	    }
	    
	    // now perform floyd-warshall's all pairs shortest path
	    for (int k = 0; k < acount; ++k) 
		for (int i = 0; i < acount; ++i) 
		    for (int j = 0; j < acount; ++j) 
			cost[i][j] = Math.min
			    (cost[i][j], cost[i][k]+cost[k][j]);
	}

	Comparator<MolBond[]> sorter = new Comparator<MolBond[]>() {
	    public int compare (MolBond[] p1, MolBond[] p2) {
		return p1.length - p2.length;
	    }
	};

	@Override
	public MolBond[][] getPaths (int a, int b) {
	    paths.clear();

	    // cost is defined in terms of atoms and paths are bonds, 
	    //   so make sure we bound it properly
	    walkPaths (this, a, b);

	    return paths.toArray(new MolBond[0][]);
	}

	@Override
	public Map<String, Integer> getSparseVector () {
	    if (sparseVector != null) {
		return sparseVector;
	    }

	    Map<String, Integer> vec = new TreeMap<String, Integer>(this);
	    int size = getAtomCount();
	    for (int i = 0; i < size; ++i) {
		for (int j = i+1; j < size; ++j) {
		    MolBond[][] paths = getPaths (i, j);

		    ArrayList<String> pset = new ArrayList<String>();
		    for (int k = 0; k < paths.length; ++k) {
			if (paths[k] != null) {
			    pset.add(toString (paths[k]));
			}
		    }

		    Collections.sort(pset, this);
		    StringBuilder sb = new StringBuilder ();
		    for (String s : pset) {
			if (sb.length() > 0) {
			    sb.append('.');
			}
			sb.append(s);
		    }
		    String s = sb.toString();
		    Integer c = vec.get(s);
		    vec.put(s, c== null ? 1 : (c+1));
		}
	    }

	    return sparseVector = vec;
	}


	@Override
	public void walkPaths (PathVisitor visitor, int start, int end) {
	    BitSet visited = new BitSet (atoms.length);
	    Stack<MolBond> path = new Stack<MolBond>();
	    dfs (visitor, path, visited, start, start, end);
	}

	@Override
	protected void dfs (PathVisitor visitor, Stack<MolBond> path, 
			    BitSet visited, int start, int a, int end) {
	    if (a == end) {
		visitor.visit(getMolecule (), path);
		return;
	    }
	    
	    visited.set(a);
	    MolAtom atom = atoms[a];
	    for (int b = 0; b < atom.getBondCount(); ++b) {
		MolBond bond = atom.getBond(b);
		int xa = mol.indexOf(bond.getOtherAtom(atom));

		if ((cost[start][xa] + cost[xa][end] <= cost[start][end])
		    && !visited.get(xa)) {
		    path.push(bond);
		    dfs (visitor, path, visited, start, xa, end);
		    path.pop();
		}
	    }
	    visited.clear(a);
	}

	public int compare (String s1, String s2) {
	    int l1 = s1.length(), l2 = s2.length();
	    int d = l1 - l2;
	    if (d == 0) {
		d = s1.compareTo(s2);
	    }
	    return d;
	}

	public void visit (Molecule mol, List<MolBond> path) {
	    paths.add(path.toArray(new MolBond[0]));
	}
    }


    public static class AllPaths extends MolPath implements PathVisitor {
	ArrayList<MolBond[]> paths = new ArrayList<MolBond[]>();

	public AllPaths () {}
	public AllPaths (Molecule m) { 
	    super (m);
	}

	public void visit (Molecule m, List<MolBond> path) {
	    paths.add(path.toArray(new MolBond[0]));
	}

	@Override
	public MolBond[][] getPaths (int a, int b) {
	    paths.clear();
	    walkPaths (this, a, b);
	    return paths.toArray(new MolBond[0][]);
	}

	public MolBond[][] getPaths (MolAtom a, MolAtom b) {
	    Molecule m = getMolecule ();
	    return getPaths (m.indexOf(a), m.indexOf(b));
	}
    }

    protected Molecule mol;
    protected MolAtom[] atoms;
    protected MolBond[] bonds;
    protected Map<String, Integer> sparseVector;
    protected MolBond[][][][] paths;
    protected Map<String, int[]> pairs;
    protected int maxDepth = MAX_DEPTH;

    public MolPath () {
    }

    public MolPath (Molecule mol) {
	setMolecule (mol);
    }

    public int getMaxDepth () { return maxDepth; }
    public void setMaxDepth (int maxDepth) { this.maxDepth = maxDepth; }

    public void setMolecule (Molecule mol) {
	this.mol = mol;
	atoms = this.mol.getAtomArray();
	bonds = this.mol.getBondArray();
	sparseVector = null;
	paths = new MolBond[atoms.length][atoms.length][][];
	pairs = new HashMap<String, int[]>();
    }

    public Molecule getMolecule () { return mol; }
    public MolAtom[] getAtoms () { return atoms; }
    public MolBond[] getBonds () { return bonds; }
    public int getAtomCount () { return atoms != null ? atoms.length : -1; }
    public int getBondCount () { return bonds != null ? bonds.length : -1; }
    public MolAtom getAtom (int i) { return atoms[i]; }
    public MolBond getBond (int j) { return bonds[j]; }

    /*
     * Given a starting and ending atoms, this method returns
     * all unique paths that connect the two atoms
     */
    public void walkPaths (PathVisitor visitor, 
			   MolAtom startAtom, MolAtom endAtom) {
	if (mol == null) {
	    throw new IllegalArgumentException ("No input molecule defined");
	}
	int start = mol.indexOf(startAtom);
	int end = mol.indexOf(endAtom);

	if (start < 0) {
	    throw new IllegalArgumentException
		("Starting atom is not part of molecule");
	}

	if (end < 0) {
	    throw new IllegalArgumentException
		("Ending atom is not part of molecule");
	}

	walkPaths (visitor, start, end);
    }


    public void walkPaths (PathVisitor visitor, int start, int end) {
	walkPaths (visitor, maxDepth, start, end);
    }

    /*
     * length <= atom count
     */
    public void walkPaths (PathVisitor visitor, int length, 
			   int start, int end) {
	if (visitor == null) {
	    throw new IllegalArgumentException ("Path visitor is null");
	}

	if (start < 0 || start >= atoms.length) {
	    throw new IllegalArgumentException
		("Invalid starting atom specified");
	}

	if (end < 0 || end >= atoms.length) {
	    throw new IllegalArgumentException
                ("Invalid ending atom specified");
	}

        if (maxDepth > 0 && length > maxDepth) {
            logger.warning("Truncating path to "+maxDepth);
            length = maxDepth;
        }

	BitSet visited = new BitSet (atoms.length);

	Stack<MolBond> path = new Stack<MolBond>();
	dfs (visitor, path, visited, length, start, end);
    }

    public synchronized Map<String, Integer> getSparseVector () { 
	if (sparseVector == null) {
	    sparseVector = calcAllPairwisePathVector (this, pairs);
	}
	return sparseVector; 
    }

    /*
     * Generate a fingerprint where numInts is the length of the 
     *  fingerprint in 32-bit units.  That is, if an 1024-bit fingerprint
     *  is desired, then numInts is 32.
     */
    public int[] generateFingerprint (int numInts) {
	return generateFingerprint (new int[numInts], 0, numInts);
    }

    public int[] generateFingerprint (int[] fp, int offset, int length) {
	for (int i = 0; i < length; ++i) {
	    fp[offset+i] = 0;
	}

	Map<String, Integer> sv = getSparseVector ();
	int nbits = length * 32;

	Random rand = new Random ();
	for (String s : sv.keySet()) {
	    long hash = (long)s.hashCode() & 0xffffffffl;
	    int bit =  (int)(hash % nbits); 
	    //System.out.println(s + " => " + bit);
	    fp[offset+bit/32] |= 1 << (31 - bit%32);
	    if (s.indexOf('.') > 0) {
		// multiple path, then turn on additional bit
		rand.setSeed(hash);
		bit = rand.nextInt(nbits);
		fp[offset+bit/32] |= 1 << (31 - bit%32);
	    }
	}
	return fp;
    }


    /*
     * Default to return all paths
     */
    public MolBond[][] getPaths (int a, int b) {
	MolBond[][] p = paths[a][b];
	if (p == null) {
	    p = new AllPaths (getMolecule ()).getPaths(a, b);
	    paths[a][b] = paths[b][a] = p;
	}
	return p;
    }

    /*
     * This method takes as input a canonical path and returns all node
     * pairs (encoded as integers) that match the given path.  Each integer
     * can be extracted to get the individual node index as follows:
     *   int node1 = (int)((i >> 16) & 0xffff);
     *   int node2 = (int)(i & 0xffff);
     * where i is an integer returned.
     */
    public int[] getPairs (String path) {
	int[] p = pairs.get(path);
	if (p == null && sparseVector == null) {
	    getSparseVector ();
	    p = pairs.get(path);
	}
	return p;
    }

    public static MolAtom[] toNodePath (MolBond[] path) {
	if (path.length == 0) {
	    return new MolAtom[0];
	}

	if (path.length == 1) {
	    return new MolAtom[]{path[0].getAtom1(), path[0].getAtom2()};
	}

	MolAtom atom = null; // figuring out direction
	if (path[0].getAtom1() == path[1].getAtom1()
	    || path[0].getAtom1() == path[1].getAtom2()) {
	    atom = path[0].getAtom2();
	}
	else if (path[0].getAtom2() == path[1].getAtom1()
		 || path[0].getAtom2() == path[1].getAtom2()) {
	    atom = path[0].getAtom1();
	}
	else {
	    throw new IllegalArgumentException ("Path is disconnected");
	}

	MolAtom[] atoms = new MolAtom[path.length+1];
	atoms[0] = atom;
	for (int i = 0; i < path.length; ++i) {
	    atom = path[i].getOtherAtom(atom);
	    if (atom == null) {
		throw new IllegalArgumentException ("Path is disconnected");
	    }
	    atoms[i+1] = atom;
	}
	return atoms;
    }

    protected static String getDefaultLabel (MolAtom atom) {
	return atom.getSymbol() + (atom.hasAromaticBond() ? "ar" 
				   : "sp"+atom.getHybridizationState());
    }

    private Annotator pathAnnotator = createDefaultAnnotator ();
    public String toString (MolBond[] path) {
	return toString (path, pathAnnotator);
    }

    public static Annotator createDefaultAnnotator () {
	return new Annotator () {
		public String getLabel (MolAtom a) {
		    return MolPath.getDefaultLabel(a);
		}
		public String getLabel (MolBond b) {
		    return "";
		}
	    };
    }

    public static String toString (MolBond[] path, Annotator anno) {
	if (path == null || path.length == 0) {
	    return "";
	}

	if (path.length == 1) {
	    String s1 = anno.getLabel(path[0].getAtom1());
	    String s2 = anno.getLabel(path[0].getAtom2());
	    return s1.compareTo(s2) <= 0 ? s1+anno.getLabel(path[0])+s2 
		: s2+anno.getLabel(path[0])+s1;
	}

	MolAtom head = null; // figuring out direction
	if (path[0].getAtom1() == path[1].getAtom1()
	    || path[0].getAtom1() == path[1].getAtom2()) {
	    head = path[0].getAtom2();
	}
	else if (path[0].getAtom2() == path[1].getAtom1()
		 || path[0].getAtom2() == path[1].getAtom2()) {
	    head = path[0].getAtom1();
	}
	else {
	    throw new IllegalArgumentException ("Path is disconnected");
	}

	MolAtom tail = null;
	{ int i = path.length - 1;
	    if (path[i].getAtom1() == path[i-1].getAtom1()
		|| path[i].getAtom1() == path[i-1].getAtom2()) {
		tail = path[i].getAtom2();
	    }
	    else if (path[i].getAtom2() == path[i-1].getAtom1()
		     || path[i].getAtom2() == path[i-1].getAtom2()) {
		tail = path[i].getAtom1();
	    }
	    else {
		throw new IllegalArgumentException ("Path is disconnected");
	    }
	}

	StringBuilder forward = new StringBuilder ();
	StringBuilder backward = new StringBuilder ();

	forward.append(anno.getLabel(head));
	backward.append(anno.getLabel(tail));
	int dir = head.getAtno() - tail.getAtno();
	for (int i = 0, j = path.length-1; i < path.length; ++i, --j) {
	    head = path[i].getOtherAtom(head);
	    forward.append(anno.getLabel(path[i]));
	    forward.append(anno.getLabel(head));

	    tail = path[j].getOtherAtom(tail);
	    backward.append(anno.getLabel(path[j]));
	    backward.append(anno.getLabel(tail));

	    if (dir == 0) {
		dir = head.getAtno() - tail.getAtno();
	    }
	}

	return dir <= 0 ? forward.toString() : backward.toString();
    }


    protected void dfs (PathVisitor visitor, Stack<MolBond> path, 
			BitSet visited, int length, int a, int end) {
	/*
	if (visited.get(a)) {
	    return;
	}
	*/

	if (a == end) {
	    visitor.visit(getMolecule (), path);
	    return;
	}

	if (length < 0 || path.size() <= length) {
	    visited.set(a);
	    MolAtom atom = atoms[a];
	    for (int b = 0; b < atom.getBondCount(); ++b) {
		MolBond bond = atom.getBond(b);
		int xa = mol.indexOf(bond.getOtherAtom(atom));
		if (!visited.get(xa)) {
		    path.push(bond);
		    dfs (visitor, path, visited, length, xa, end);
		    path.pop();
		}
	    }
	    visited.clear(a);
	}
    }

    public static Map<String, Integer> calcAllPairwisePathVector (MolPath mp) {
	return calcAllPairwisePathVector (mp, null);
    }

    public static Map<String, Integer> 
	calcAllPairwisePathVector (MolPath mp, Map<String, int[]> pairs) {
	Map<String, Integer> vec = new TreeMap<String, Integer>();
	Map<String, Set<Integer>> pv = new TreeMap<String, Set<Integer>>();

	int size = mp.getAtomCount();
	for (int i = 0; i < size; ++i) {
	    for (int j = i+1; j < size; ++j) {
		MolBond[][] paths = mp.getPaths(i, j);
		//System.out.println("** Paths between "+(i+1)+" and "+(j+1));
		for (int k = 0; k < paths.length; ++k) {
		    if (paths[k] != null) {
			String s = mp.toString(paths[k]);
			Integer c = vec.get(s);
			vec.put(s, c== null ? 1 : (c+1));

			Set<Integer> iv = pv.get(s);
			if (iv == null) {
			    iv = new HashSet<Integer>();
			    pv.put(s, iv);
			}

			iv.add((int)((i << 16) | j));
		    }
		}
	    }
	}

	if (pairs != null) {
	    for (Map.Entry<String, Set<Integer>> e : pv.entrySet()) {
		String s = e.getKey();
		Set<Integer> iv = e.getValue();
		int[] ia = new int[iv.size()];
		{ int i = 0;
		    for (Iterator<Integer> it = iv.iterator(); 
			 it.hasNext(); ) {
			ia[i++] = it.next();
		    }
		}
		pairs.put(s, ia);
	    }
	}

	return vec;
    }

    public static Map<String, Integer> calcTerminalPathVector (MolPath mp) {
	int size = mp.getAtomCount();

	ArrayList<Integer> terminals = new ArrayList<Integer>();
	for (int i = 0; i < mp.getAtomCount(); ++i) {
	    MolAtom a = mp.getAtom(i);
	    if (a.isTerminalAtom()) {
		terminals.add(i);
	    }
	}

	Map<String, Integer> vec = new HashMap<String, Integer>();
	for (int i = 0; i < terminals.size(); ++i) {
	    int a = terminals.get(i);
	    for (int j = i+1; j < terminals.size(); ++j) {
		int b = terminals.get(j);
		MolBond[][] paths = mp.getPaths(a, b);
		for (int k = 0; k < paths.length; ++k) {
		    if (paths[k] == null) {
			continue;
		    }
		    String s = mp.toString(paths[k]);
		    Integer c = vec.get(s);
		    vec.put(s, c== null ? 1 : (c+1));
		}
	    }
	}		

	return vec;
    }


    public static void test2 (String[] argv) throws Exception {
	String s1 = "CCCCCCCCCCCC1=NC(=Cc2[nH]c(cc2OC)-c2cc3ccccc3[nH]2)C=C1";
	String s2 = "CCN1C(SC(=C1c1ccccc1)c1ccccc1)=CC=Cc1sc2c(ccc3ccccc23)[n+]1CC";
	Molecule m1 = new MolHandler (s1).getMolecule();
	Molecule m2 = new MolHandler (s2).getMolecule();
	Map<String, Integer> v1 = calcAllPairwisePathVector (new AllPaths (m1));
	Map<String, Integer> v2 = calcAllPairwisePathVector (new AllPaths (m2));
	System.out.println(v1.size() + " vs "+v2.size());
	//System.out.println("dot = "+ SparseVectorMetricFactory.cosine(v1, v2));
    }

    public static void test1 (String[] argv) throws Exception {
	MolHandler mh = new MolHandler ("CCN1\\C(SC(=C1C2=CC=CC=C2)C3=CC=CC=C3)=C/C=C/C4=[N+](CC)C5=C(S4)C6=CC=CC=C6C=C5");

	Molecule mol = mh.getMolecule();
	mol.aromatize();

	MolAtom[] atoms = mol.getAtomArray();
	MolPath.AllPaths mu = new MolPath.AllPaths (mol);

	PathVisitor visitor = new PathVisitor () {
		public void visit (Molecule mol, List<MolBond> path) {
		    logger.info("** Found path");
		    for (MolBond b : path) {
			System.err.println
			    (" "+(mol.indexOf(b.getAtom1())+1)+"-"
			     +(mol.indexOf(b.getAtom2())+1));
		    }
		}
	    };

	Map<String, Integer> frags = new TreeMap<String, Integer>();
	for (int i = 0; i < atoms.length; ++i) {
	    for (int j = i+1; j < atoms.length; ++j) {
		MolBond[][] paths = mu.getPaths(i, j);
		/*
		logger.info("## Shortest paths between atom "+(i+1)
			    +" and "+(j+1) + ": "+paths.length);
		*/
		for (int k = 0; k < paths.length; ++k) {
		    //System.err.print("  "+(k+1)+":");
		    /*
		    for (int p = 0; p < paths[k].length; ++p) {
			System.err.print(" "+(mol.indexOf(paths[k][p].getAtom1())+1)+"-"+(mol.indexOf(paths[k][p].getAtom2())+1));
		    }
		    */
		    //System.err.println(" "+toString (paths[k]));
		    String s = toString (paths[k], createDefaultAnnotator ());
		    Integer c = frags.get(s);
		    frags.put(s, c== null ? 1 : (c+1));
		}
	    } 
	}
	System.out.println(frags.size() + " unique fragments!");
	for (Map.Entry<String, Integer> me : frags.entrySet()) {
	    System.out.println(me.getKey() +":"+me.getValue());
	}
    }

    static void process (InputStream is) throws Exception {
	MolPath mp = new ShortestPaths ();
	MolImporter mi = new MolImporter (is);

	String q = System.getProperty("query", "c1ccccc1");
	System.out.println("Query: " + q);
	Molecule query = new MolHandler (q).getMolecule();
	query.aromatize();
	query.calcHybridization();

	int fpsize = 16;

	mp.setMolecule(query);
	int[] qFp1 = mp.generateFingerprint(fpsize);

	MolHandler mh = new MolHandler (query);
	int[] qFp2 = mh.generateFingerprintInInts(fpsize, 2, 6);

	MolSearch ms = new MolSearch ();
	ms.setQuery(query);
	{ int n1 = 0, n2 = 0;
	    for (int i = 0; i < fpsize; ++i) {
		n1 += Integer.bitCount(qFp1[i]);
		n2 += Integer.bitCount(qFp2[i]);
	    }
	    System.out.println("Query bit count: "+n1 + " " +n2);
	    System.out.println("Query paths");
	    for (Map.Entry<String, Integer> e : 
		     mp.getSparseVector().entrySet()) {
		System.out.println(e.getKey() + " " + e.getValue());
	    }
	}

	int[] fp1 = new int[fpsize];
	double[] avgBits1 = new double[fpsize];
	double[] avgBits2 = new double[fpsize]; 
	double avg1 = 0., avg2 = 0.;
	int total = 0, nmatches = 0;
	int fpc1 = 0, fpc2 = 0; // false positive count
	for (Molecule mol = new Molecule (); mi.read(mol); ) {
	    mol.calcHybridization();
	    mol.aromatize();

	    ms.setTarget(mol);
	    boolean matched = ms.isMatching();

	    mp.setMolecule(mol);
	    mp.generateFingerprint(fp1, 0, fp1.length);

	    mh.setMolecule(mol);
	    int[] fp2 = mh.generateFingerprintInInts(fpsize, 2, 6);

	    int c1 = 0, c2 = 0;
	    for (int i = 0; i < fpsize; ++i) {
		if ((fp1[i] & qFp1[i]) == qFp1[i]) {
		    ++c1;
		}
		if ((fp2[i] & qFp2[i]) == qFp2[i]) {
		    ++c2;
		}

		int a = Integer.bitCount(fp1[i]);
		int b = Integer.bitCount(fp2[i]);
		avgBits1[i] += a;
		avgBits2[i] += b;
		avg1 += a;
		avg2 += b;
	    }

	    if (c1 == fp1.length) {
		if (!matched /*&& c2 != fp2.length*/) {
		    /*
		    System.err.println
			("** false positive found for " 
			 + mol.toFormat("smiles:q") + " " +mol.getName());
		    Map<String, Integer> target = mp.getSparseVector();
		    mp.setMolecule(query);
		    Map<String, Integer> qq = mp.getSparseVector();
		    target.keySet().retainAll(qq.keySet());
		    System.err.println("overlap paths");
		    for (Map.Entry<String, Integer> e : target.entrySet()) {
			System.err.println("  "+e.getKey() + " " + e.getValue());
		    }
		    */
		    
		    ++fpc1;
		}
	    }
	    else if (matched) { // false negative
		System.err.println
		    ("** fatal error: false negative found for " 
		     + mol.toFormat("smiles:q") + " " +mol.getName());
		Map<String, Integer> target = mp.getSparseVector();
		for (Map.Entry<String, Integer> e : target.entrySet()) {
		    System.err.println("  "+e.getKey() + " " + e.getValue());
		}
		
		mp.setMolecule(query);
		Map<String, Integer> qq = mp.getSparseVector();
		qq.keySet().removeAll(target.keySet());
		System.err.println("missing paths");
		for (Map.Entry<String, Integer> e : qq.entrySet()) {
		    System.err.println("  "+e.getKey() + " " + e.getValue());
		}
		System.exit(1);
	    }

	    if (c2 == fp2.length) {
		if (!matched) ++fpc2;
	    }

	    if (matched) {
		++nmatches;
	    }
	    ++total;
	}
	System.out.println("Average bits per bin for "+total);
	for (int i = 0; i < fp1.length; ++i) {
	    System.out.printf("%1$2d: %2$.3f %3$.3f\n", 
			      i+1, avgBits1[i]/total, avgBits2[i]/total);
	}
	System.out.println("Average bit per molecule: " 
			   + String.format("%1$.3f", avg1/total)
			   + " " + String.format("%1$.3f", avg2/total));
	System.out.println("False positive rate: "
			   +String.format("%1$.3f", (double)fpc1/total)
			   +"/" +fpc1 + " " 
			   + String.format("%1$.3f", (double)fpc2/total)
			   + "/" + fpc2);
	System.out.println("Num matches: "+nmatches+"/"+total);
    }

    static void doMolPath (InputStream is) throws Exception {
	System.out.println("Generating MolPath fingerprint...");
	MolPath mp = new ShortestPaths ();
	MolImporter mi = new MolImporter (is);
	int size = 16, count = 0;
	long start = System.currentTimeMillis();
	for (Molecule mol = new Molecule (); mi.read(mol); ) {
	    mol.calcHybridization();
	    mol.aromatize();

	    mp.setMolecule(mol);
	    int[]fp = mp.generateFingerprint(size);

	    //new gov.nih.ncgc.descriptor.AtomPair(mol);
	    ++count;
	}
	double time = (System.currentTimeMillis() - start)*1e-3;
	System.out.println("Total time to generate fp: "
			   +String.format("%1$.3fs", time));
	System.out.println("Average time per structure: "+
			   String.format("%1$.0fms", time*1e3/count));
    }

    static void doFp (InputStream is) throws Exception {
	System.out.println("Generating path-based fingerprint...");
	MolHandler mh = new MolHandler ();
	MolImporter mi = new MolImporter (is);
	int size = 16, count = 0;
	long start = System.currentTimeMillis();
	for (Molecule mol = new Molecule (); mi.read(mol); ) {
	    mol.aromatize();
	    mh.setMolecule(mol);
	    int[]fp = mh.generateFingerprintInInts(size, 2, 6);
	    ++count;
	}
	double time = (System.currentTimeMillis() - start)*1e-3;
	System.out.println("Total time to generate fp: "
			   +String.format("%1$.3fs", time));
	System.out.println("Average time per structure: "+
			   String.format("%1$.0fms", time*1e3/count));
    }
    
    public static void main (String[] argv) throws Exception {
	
	if (argv.length == 0) {
	    logger.info("Reading from STDIN...");
	    process (System.in);
	}
	else {
	    for (String file : argv) {
		FileInputStream fis = new FileInputStream (file);
		System.out.println("## "+file);
		//doMolPath (fis);
		//doFp (fis);
		process (fis);
		fis.close();
	    }
	}
    }
}
