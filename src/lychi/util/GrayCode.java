// $Id: GrayCode.java 2707 2009-06-22 09:29:24Z nguyenda $

package lychi.util;

import java.util.Observer;
import java.util.Observable;

// rip this off of wikipedia; knuth's vol. 4 algorithm M doesn't require
//   all these auxilary arrays.... don't have it with me right now.
public class GrayCode extends Observable {

    private int maxsize = 0, size = 0;
    private int[] n, g, u, c;

    public GrayCode (int N, int k) {
	init (N, k);
    }

    public GrayCode (int[] N) {
	init (N);
    }

    protected void init (int N, int k) {
	n = new int[k+1];
	g = new int[k+1];
	u = new int[k+1];
	c = new int[k]; // copy of g

	for (int i = 0; i <= k; ++i) {
	    g[i] = 0;
	    u[i] = 1;
	    n[i] = N;
	}
	size = 0;
    }

    protected void init (int[] N) {
	int k = N.length;
	n = new int[k+1];
	g = new int[k+1];
	u = new int[k+1];
	c = new int[k]; // copy of g

	int min = Integer.MAX_VALUE;
	for (int i = 0; i < k; ++i) {
	    g[i] = 0;
	    u[i] = 1;
	    n[i] = N[i];
	    if (N[i] < min) {
		min = N[i];
	    }
	}
	g[k] = 0;
	u[k] = 1;
	n[k] = min;
	size = 0;
    }

    public void generate () {
	for(int i, j; g[c.length] == 0;) {
	    System.arraycopy(g, 0, c, 0, c.length);

	    setChanged ();
	    notifyObservers (c);
	    
	    i = 0; 
	    j = g[0] + u[0];
	    while (((j >= n[i]) || (j < 0)) && (i < c.length)) {
		u[i] = -u[i];
		++i;
		j = g[i] + u[i];
	    }
	    g[i] = j;

	    ++size;
	    if (countObservers() == 0 
		|| (maxsize > 0 && size >= maxsize)) {
		break;
	    }
	}
    }

    public void setMaxSize (int maxsize) {
	this.maxsize = maxsize;
    }
    public int getMaxSize () { return maxsize; }

    public int size () { return size; }

    public static GrayCode createBinaryGrayCode (int size) {
	return new GrayCode (2, size);
    }

    public static class Enum {
        public static void main (final String[] argv) throws Exception {
            // all possible subsets = 2^k
            GrayCode g = createBinaryGrayCode (argv.length);
            g.addObserver(new Observer () {
                    public void update (Observable o, Object arg) {
                        int[] c = (int[])arg;
                        int j = 0;
                        for (int i = 0; i < c.length; ++i) {
                            if (c[i] != 0) {
                                if (j == 0) {
                                    System.out.print("[");
                                }
                                else if (j > 0) {
                                    System.out.print(",");
                                }
                                System.out.print(argv[i]);
                                ++j;
                            }
                        }
                        if (j > 0) {
                            System.out.println("]");
                        }
                    }
                });
            g.generate();
        }
    }

    public static void main (String[] argv) throws Exception {
	if (argv.length == 0) {
	    System.out.println("Usage: GrayCode [N k | n0 n1...]\n");
	    System.exit(1);
	}

	/*
	if (argv.length >= 2) {
	    final int N = Integer.parseInt(argv[0]);
	    int k = Integer.parseInt(argv[1]);

	    System.out.println("Gray Code ("+N+","+k+")");
	    GrayCode g = new GrayCode (N, k);
	    Observer obs = new Observer () {
		    int index = 0;
		    public void update (Observable o, Object arg) {
			int[] c = (int[])arg;
			System.out.printf("%1$5d: ", ++index);
			if (N == 2 && c.length < 64) {
			    long x = 0;
			    for (int i = 0; i < c.length; ++i) {
				System.out.print(c[i]+ " ");
				if (c[i] == 1) {
				    x |= 1 << (c.length-i-1);
				}
			    }
			    System.out.print("=> "+x);
			}
			else {
			    System.out.print(c[0]);
			    for (int i = 1; i < c.length; ++i) {
				System.out.print(" " + c[i]);
			    }
			}
			System.out.println();
		    }
		};
	    g.addObserver(obs);
	    g.generate();
	    System.out.println("--");
	}
	*/

	int[] a = new int[argv.length];
	for (int i = 0; i < a.length; ++i) {
	    a[i] = Integer.parseInt(argv[i]);
	}
	GrayCode g = new GrayCode (a);
	g.addObserver(new Observer () {
		public void update (Observable o, Object arg) {
		    int[] c = (int[])arg;
		    System.out.print(c[0]);
		    for (int i = 1; i < c.length; ++i) {
			System.out.print(" " + c[i]);
		    }
		    System.out.println();
		}
	    });
	g.generate();
	System.out.println("** "+g.size()+" values generated!");
    }
}
