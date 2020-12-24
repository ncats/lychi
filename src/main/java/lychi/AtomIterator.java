package lychi;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.struc.*;
import chemaxon.formats.*;

/**
 * A utility class to iterate the atoms in the same order (regardless
 * how they are ordered in the parent molecule).  Any changes to this
 * class that result in changes in the underlying order will have major
 * consequences in any function that depends on this ordering, e.g.,
 * hash key generation!  So be careful changing this class!!!!
 */
public class AtomIterator implements ListIterator<MolAtom> {
    private static final Logger logger = 
    Logger.getLogger(AtomIterator.class.getName());

    static class IntPair implements Comparable<IntPair> {
        int index, value;
        IntPair (int index, int value) {
            this.index = index;
            this.value = value;
        }
        public int compareTo (IntPair ip) {
            int d = value - ip.value;
            if (d == 0) {
                d = index - ip.index;
            }
            return d;
        }
        public String toString () {
            return "["+index+","+value+"]";
        }
    }

    protected MolAtom[] atoms;
    protected ListIterator<Integer> iter;
    protected List<Integer> index;

    public AtomIterator (MoleculeGraph mol) {
        setMolecule (mol);
    }

    public AtomIterator (MoleculeGraph mol, int[] rank) {
        setMolecule (mol, rank);
    }
    
    public void setMolecule (MoleculeGraph mol) {
        int[] rank = new int[mol.getAtomCount()];
        mol.getGrinv(rank, MoleculeGraph.GRINV_NOHYDROGEN);
        //mol.calcDehydrogenizedGrinv(rank);
        setMolecule (mol, rank);
    }
    
    public void setMolecule (MoleculeGraph mol, int[] r) {
        atoms = mol.getAtomArray();
        
        IntPair[] ip = new IntPair[atoms.length];
        int[] rank = new int[atoms.length]; // make a copy
        for (int i = 0; i < ip.length; ++i) {
            rank[i] = r[i];
            ip[i] = new IntPair (i, rank[i]);
        }
        Arrays.sort(ip);
        
        /*
          System.out.println(((Molecule)mol).toFormat("smiles:q"));
          for (IntPair i : ip) {
          System.out.println(i+": "+atoms[i.index].getSymbol());
          }
        */
        
        index = new ArrayList<Integer>();
        for (int i = 0; i < atoms.length; ++i) {
            dfs (mol, ip[i].index, index, rank);
        }
        ip = null;
        iter = index.listIterator();
    }
    
    /*
     * Be careful changing this function!!!
     */
    void dfs (MoleculeGraph mol, int pos, List<Integer> index, int[] rank) {
        if (rank[pos] < Integer.MAX_VALUE) {
            index.add(pos);
            rank[pos] = Integer.MAX_VALUE;
            
            // find the lowest non-visited atom
            int k = mol.getNeighborCount(pos);
            if (k > 0) {
                int j = mol.getNeighbor(pos, 0);
                for (int i = 1; i < k; ++i) {
                    int nb = mol.getNeighbor(pos, i);
                    if (rank[nb] == rank[j]) {
                        if (nb < j)
                            j = nb;
                    }
                    else if (rank[nb] < rank[j]) {
                        j = nb;
                    }
                }
                dfs (mol, j, index, rank);
            }
        }
    }

    /*
     * ListIterator interface
     */
    public boolean hasNext () { return iter.hasNext(); }
    public boolean hasPrevious () { return iter.hasPrevious(); }
    public MolAtom next () { return atoms[iter.next()]; }
    public MolAtom previous () { return atoms[iter.previous()]; }
    public int nextIndex () { return index.get(iter.nextIndex()); }
    public int previousIndex () { return index.get(iter.previousIndex()); }

    /*
     * Unsupported operations
     */
    public void add (MolAtom atom) {
        throw new UnsupportedOperationException ("add is not supported");
    }
    public void remove () {
        throw new UnsupportedOperationException ("remove is not supported");
    }
    public void set (MolAtom atom) {
        throw new UnsupportedOperationException ("set is not supported");
    }
    
    public static void main (String[] argv) throws Exception {
        for (String s : argv) {
            MolImporter mi = new MolImporter (s);
            for (Molecule m = new Molecule (); mi.read(m); ) {
                System.out.print(m.getName()+":");
                AtomIterator ai = new AtomIterator (m);
                while (ai.hasNext()) {
                    MolAtom a = ai.next();
                    System.out.print(" "+a.getAtomMap());
                }
                System.out.println();
            }
        }
    }
}

