package lychi.util;

import java.util.Random;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.security.SecureRandom;

import java.util.concurrent.*;

public class Base32 {
    private static final Logger logger = 
        Logger.getLogger(Base32.class.getName());

    private static final char[] ALPHA = {
        'A','B','C','D','F','G','H','J','K',
        'L','M','N','P','Q','R','S','T','U',
        'V','W','X','Y','Z','1','2','3','4',
        '5','6','7','8','9',
    };

    /*
     * 2^15 = 32768 = 32^3
     * 
     * A triple in base-32 represents 15 bits.
     */
    /*
    static final String[] CODER;
    static {
        String[] lut = new String[32768];
        char[] triplet = new char[3];
        int p = 0, size = ALPHA.length;
        for (int i = 0; i < size; ++i) {
            triplet[0] = ALPHA[i];
            for (int j = 0; j < size; ++j) {
                triplet[1] = ALPHA[j];
                for (int k = 0; k < size; ++k) {
                    triplet[2] = ALPHA[k];
                    lut[p++] = new String (triplet);
                }
            }
        }
        logger.info("## generated "+p+" trigams!");

        CODER = lut;
    }

    static void testUnrank () {
        for (int i = 0; i < CODER.length; ++i) {
            String unrank = unrank (i, 15);
            long rank = rank (unrank);
            System.out.println(String.format("%1$5d", i)+"\t"+CODER[i]
                               +"\t"+unrank+"\t"+rank);
            if (!unrank.equals(CODER[i])) {
                throw new IllegalStateException ("Unrank fails!");
            }
            if (i != rank) {
                throw new IllegalStateException ("Rank fails!");
            }
        }
    }
    */

    public static String unrank (long rank) {
        return unrank (rank, 60);
    }

    public static String unrank (long rank, int size) {
        return unrank (new StringBuilder (), 
                       rank, size > 60 ? 60 : size).toString();
    }

    public static StringBuilder unrank 
        (StringBuilder sb, long rank, int size) {
        for (int r = 0; r < size; r += 5) {
            int i = (int)(rank & 0x1f);
            sb.insert(0, ALPHA[i]);
            rank >>>= 5;
        }
        return sb;
    }

    public static long rank (String str) {
        str = str.toUpperCase();
        long rank = 0, k = 0;
        for (int i = str.length(); --i >= 0; k += 5) {
            char ch = str.charAt(i);
            int j = 0;
            for (; j < ALPHA.length; ++j) 
                if (ALPHA[j] == ch) break;
            if (j == ALPHA.length)
                return -1;
            rank |= (long)j << k;
        }

        return rank;
    }

    public static String encode15 (byte... data) {
	if (data.length < 2) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode15 (new StringBuilder (3), data).toString();
    }

    public static StringBuilder encode15 (StringBuilder sb, byte... data) {
	int b3 = (data[0]&0xff)|((data[1]&0x7f)<<8); //8+7=15[1]
	return sb.append(unrank (b3, 15));
    }

    public static String encode20 (byte... data) {
        if (data.length < 3) {
	    throw new IllegalArgumentException ("Not enough data specified");
	}
	return encode20 (new StringBuilder (4), data).toString();
    }

    public static StringBuilder encode20 (StringBuilder sb, byte... data) {
        int b4 = ((data[1]&0x80)>>>7) | ((data[2] & 0x0f)<<1); //1+4=20[4]
        sb.append(ALPHA[b4]);
        return encode15 (sb, data);
    }

    public static String encode25 (byte... data) {
        if (data.length < 4) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode25 (new StringBuilder (5), data).toString();
    }

    public static StringBuilder encode25 (StringBuilder sb, byte... data) {
        int b5 = ((data[2]&0xf0)>>4) | ((data[3] & 0x01)<<4); //4+1=25[7]
        sb.append(ALPHA[b5]);
        return encode20 (sb, data);
    }

    public static String encode30 (byte... data) {
        if (data.length < 4) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode30 (new StringBuilder (6), data).toString();        
    }

    public static StringBuilder encode30 (StringBuilder sb, byte... data) {
        int b6 = (data[3] & 0x3f) >>> 1; //0+5=30[2]
        sb.append(ALPHA[b6]);
        return encode25 (sb, data);
    }

    public static String encode35 (byte... data) {
        if (data.length < 5) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode35 (new StringBuilder (7), data).toString();        
    }

    public static StringBuilder encode35 (StringBuilder sb, byte... data) {
        int b7 = ((data[3] & 0xc0)>>>6) | ((data[4]&0x07)<<2); //2+3=35[5]
        sb.append(ALPHA[b7]);
        return encode30 (sb, data);
    }

    public static String encode40 (byte... data) {
        if (data.length < 5) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode40 (new StringBuilder (8), data).toString();        
    }

    public static StringBuilder encode40 (StringBuilder sb, byte... data) {
        int b8 = (data[4] & 0xf8)>>>3; //5+0=40[0]
        sb.append(ALPHA[b8]);
        return encode35 (sb, data);
    }

    public static String encode45 (byte... data) {
        if (data.length < 6) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode45 (new StringBuilder (9), data).toString();
    }

    public static StringBuilder encode45 (StringBuilder sb, byte... data) {
        int b9 = data[5] & 0x1f; //0+5=45[3]
        sb.append(ALPHA[b9]);
        return encode40 (sb, data);
    }

    public static String encode50 (byte... data) {
        if (data.length < 7) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode50 (new StringBuilder (10), data).toString();
    }

    public static StringBuilder encode50 (StringBuilder sb, byte... data) {
        int b10 = ((data[6]&0x03)<<3) | ((data[5] & 0xe0)>>>5); //2+3=50[6]
        sb.append(ALPHA[b10]);
        return encode45 (sb, data);
    }

    public static String encode55 (byte... data) {
        if (data.length < 7) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode55 (new StringBuilder (11), data).toString();
    }

    public static StringBuilder encode55 (StringBuilder sb, byte... data) {
        int b11 = (data[6] & 0x7f)>>>2; //5+0=55[1]
        sb.append(ALPHA[b11]);
        return encode50 (sb, data);
    }

    public static String encode60 (byte... data) {
        if (data.length < 8) {
            throw new IllegalArgumentException ("Not enough data specified");
        }
        return encode60 (new StringBuilder (12), data).toString();
    }

    public static StringBuilder encode60 (StringBuilder sb, byte... data) {
        int b12 = ((data[7] & 0x0f)<<1) | ((data[6] & 0x80)>>>7); //4+1=60[4]
        sb.append(ALPHA[b12]);
        return encode55 (sb, data);
    }

    public static String encode30 (int x) {
        byte[] data = new byte[4];
        data[3] = (byte)((x & 0x3fffffff)>>>24);
        data[2] = (byte)((x & 0x00ffffff)>>>16);
        data[1] = (byte)((x & 0x0000ffff)>>>8);
        data[0] = (byte)(x & 0xff);
        return encode30 (data);
    }

    public static String encode40 (long x) {
        byte[] data = new byte[5];
        data[4] = (byte)((x & 0xffffffffffl)>>>32);
        data[3] = (byte)((x & 0x00ffffffffl)>>>24);
        data[2] = (byte)((x & 0x0000ffffffl)>>>16);
        data[1] = (byte)((x & 0x000000ffffl)>>>8);
        data[0] = (byte)(x & 0xff);
        return encode40 (data);
    }


    /*******************************************************************
     ** TESTING
     ******************************************************************/
    static class Test40 implements Callable<Long> {
        long start, end;

        Test40 (long start, long end) {
            this.start = start;
            this.end = end;
        }

        public Long call () throws Exception {
            logger.info("## "+Thread.currentThread().getName()
                        +": start="+start+" end="+end);
            return test40 (start, end);
        }
    }

    static class TestWorker implements Callable<Long> {
        long start, end;
        int bits, mask;

        TestWorker (long start, long end, int bits, int mask) {
            this.start = start;
            this.end = end;
            this.bits = bits;
            this.mask = mask;
        }

        public Long call () throws Exception {
            logger.info("## "+Thread.currentThread().getName()
                        +": start="+start+" end="+end+" bits="
                        +bits+" mask="+mask);
            return test (start, end, bits, mask);
        }
    }

    static class TestRandomWorker implements Callable<Void> {
        long size;
        int bits, mask;

        TestRandomWorker (long size, int bits, int mask) {
            this.size = size;
            this.bits = bits;
            this.mask = mask;
        }

        public Void call () throws Exception {
            logger.info("## "+Thread.currentThread().getName()
                        +": size="+size+" bits="+bits+" mask="+mask);
            testRandom (size, bits, mask);
            return null;
        }
    }

    static void testRandom (long size, int bits, int mask) throws Exception {
        SecureRandom rand = new SecureRandom ();
        byte[] data = new byte[(bits+7)/8];

        logger.info(Thread.currentThread().getName()+
                    ": Sampling "+size +" collision testing for encode"+bits);
        int index = data.length - 1;
        int pct = 0;
        for (long i = 0; i < size; ++i) {
            rand.nextBytes(data);

            long x = 0l;
            data[index] = (byte)(data[index] & mask);
            for (int j = 0, l = 0; j < data.length; ++j) {
                x |= (data[j] & 0xffl) << l;
                l += 8;
            }

            java.lang.reflect.Method method = Base32.class.getMethod
                ("encode"+bits, byte[].class);
            String value = (String) method.invoke(null, data);

            long k = rank (value);
            if (k != x) {
                /*
                logger.log(Level.SEVERE, "Encoded value "
                           +value+" maps to "+k+" but expected "+x);
                System.exit(1);
                */
                throw new IllegalStateException
                    ("Encoded value "+value+" maps to "+k+" but expected "+x);
            }

            String inv = unrank (k, bits);
            if (!inv.equals(value)) {
                throw new IllegalStateException
                    ("Unrank failed to map "+k+"; expecting "+value
                     +" but got "+inv);
            }

            int p = (int)(i*100./size+.5);
            if (p - pct >= 5) {
                logger.info(Thread.currentThread().getName()
                            +": "+x+" <=> "+value+"..."+p+"%");
                pct = p;
            }
        }
        logger.info(Thread.currentThread().getName()
                    +": Sweet, no collisions found in "+size+" random tests!");
    }


    static long test (long start, long end, int bits, int mask) 
        throws Exception {
        byte[] data = new byte[(bits+7)/8];

        int index = data.length - 1;
        int pct = 0;
        for (long i = start; i < end; ++i) {

            long m = 0xffl;
            for (int j = 0, l = 0; j < data.length; ++j) {
                data[j] = (byte)((i & m) >>> l);
                m <<= 8;
                l += 8;
            }
            data[index] &= mask;

            java.lang.reflect.Method method = Base32.class.getMethod
                ("encode"+bits, byte[].class);
            String value = (String) method.invoke(null, data);

            long k = rank (value);
            if (k != i) {
                throw new IllegalStateException
                    ("Encoded value "+value+" maps to "+k+" but expected "+i);
            }

            String inv = unrank (k, bits);
            if (!inv.equals(value)) {
                throw new IllegalStateException
                    ("Unrank failed to map "+k+"; expecting "+value
                     +" but got "+inv);
            }

            int p = (int)(((double)(i-start)/(end-start))*100.+.5);
            if (p - pct >= 5) {
                logger.info(Thread.currentThread().getName()
                            +": "+i+" <=> "+value+"..."+p+"%");
                pct = p;
            }
        }
        logger.info(Thread.currentThread().getName()
                    +": Sweet, no collisions found in "+(end-start)+" tests!");
        return end - start;
    }

    static void test20 () {
        //logger.info("Exact collision testing...");
        int size = 1<<20; // 2^20
        byte[] data = new byte[4];
        for (int i = 0; i < size; ++i) {
            // big endian
            data[3] = (byte)((i & 0x3fffffff)>>>24);
            data[2] = (byte)((i & 0x00ffffff)>>>16);
            data[1] = (byte)((i & 0x0000ffff)>>>8);
            data[0] = (byte)(i & 0xff);

            int j = (data[0] & 0xff) 
                | ((data[1] & 0xff)<<8) 
                | ((data[2] & 0xff)<<16)
                | ((data[3] & 0x3f)<<24);
            if (i != j) {
                throw new IllegalStateException ("Bogus decoding");
            }

            String value = encode20 (data);
            long k = rank (value);
            System.out.println(value+" "+k);
            if (k != i) {
                logger.log(Level.SEVERE, "Encoded value "
                           +value+" maps to "+k+"; expected "+i);
                System.exit(1);
            }
        }
        //logger.info("Sweet, no collisions found in "+size+" tests!");
        System.exit(0);
    }

    static void test25 () {
        System.out.print("Exact collision testing for encode25");
        long start = 0, end = 1l<<25; // 2^20

        byte[] data = new byte[4];
        int pct = 0;
        for (long i = start; i < end; ++i) {
            // big endian
            data[3] = (byte)((i & 0x01ffffffl)>>>24);
            data[2] = (byte)((i & 0x00ffffffl)>>>16);
            data[1] = (byte)((i & 0x0000ffffl)>>>8);
            data[0] = (byte)( i & 0xff);
            
            String value = encode25 (data);
            long k = rank (value);
            if (k != i) {
                logger.log(Level.SEVERE, "Encoded value "
                           +value+" maps to "+k+" but expected "+i);
                System.exit(1);
            }

            int p = (int)(i*100./end+.5);
            if (p - pct >= 5) {
                System.out.print("..."+p+"%");
                pct = p;
            }
        }
        System.out.println("\nSweet, no collisions found in "+end+" tests!");
    }

    static void test40 () {
        test40 (0l, 1l<<40);
    }

    static long test40 (long start, long end) {
        //System.out.print("Exact collision testing for encode40");
        byte[] data = new byte[5];

        int pct = 0;
        for (long i = start; i < end; ++i) {
            // big endian
            data[4] = (byte)((i & 0xffffffffffl)>>>32);
            data[3] = (byte)((i & 0x00ffffffffl)>>>24);
            data[2] = (byte)((i & 0x0000ffffffl)>>>16);
            data[1] = (byte)((i & 0x000000ffffl)>>>8);
            data[0] = (byte)( i & 0xff);
            
            String value = encode40 (data);
            long k = rank (value);
            if (k != i) {
                /*
                logger.log(Level.SEVERE, "Encoded value "
                           +value+" maps to "+k+" but expected "+i);
                */
                throw new IllegalStateException
                    ("Encoded value "+value+" maps to "+k+" but expected "+i);
            }

            String inv = unrank (k, 40);
            if (!inv.equals(value)) {
                throw new IllegalStateException
                    ("Unrank failed to map "+k+"; expecting "+value
                     +" but got "+inv);
            }

            int p = (int)(((double)(i-start)/(end-start))*100.+.5);
            if (p - pct >= 5) {
                logger.info(Thread.currentThread().getName()
                            +": "+i+" <=> "+value+"..."+p+"%");
                pct = p;
            }
        }
        logger.info(Thread.currentThread().getName()
                    +": Sweet, no collisions found in "+(end-start)+" tests!");
        return end-start;
    }

    static String toHex (byte[] data) {
        StringBuilder sb = new StringBuilder ();
        for (int j = 0; j < data.length; ++j) {
            sb.append(String.format("%1$02x", data[j] & 0xff));
        }
        return sb.toString();
    }

    static void testEndian () {
        byte[] data = new byte[4];
        for (int i=0;i<100;i++){
            data[3] = (byte)((i & 0x3fffffff)>>>24);
            data[2] = (byte)((i & 0x00ffffff)>>>16);
            data[1] = (byte)((i & 0x0000ffff)>>>8);
            data[0] = (byte)(i & 0xff);

            java.nio.ByteBuffer big = java.nio.ByteBuffer.allocate(4);
            big.order(java.nio.ByteOrder.LITTLE_ENDIAN);
            big.putInt(i);

            java.nio.ByteBuffer little = java.nio.ByteBuffer.allocate(4);
            little.order(java.nio.ByteOrder.BIG_ENDIAN);
            little.putInt(i);

            System.out.println("default:"+toHex(data)
                               +" little:"+toHex(little.array())
                               +" big:"+toHex(big.array()));
            
            String hash=encode30(data);
            int j = (data[0]&0xff) | ((data[1]&0xff)<<8) 
                | ((data[2]&0xff)<<16) | ((data[3]&0x3f)<<24);
            System.out.println(hash+"\t"+j);
        }
    }

    static void testAll (long size) throws Exception {
        logger.info("## Testing collisions for "+size+" values...");

        final int[][] codec = new int[][]{
            {15, 0x7f},
            {20, 0x0f},
            {25, 0x01},
            {30, 0x3f},
            {35, 0x07},
            {40, 0xff},
            {45, 0x1f},
            {50, 0x03},
            {55, 0x7f},
            {60, 0x0f}
        };

        for (int i = 0; i < codec.length; ++i) {
            testRandom (size, codec[i][0], codec[i][1]);
        }
    }

    static class Decode {
        public static void main (String[] argv) throws Exception {
            for (String a : argv) {
                System.out.println(String.format("%1$9s",a)+": "+rank (a));
            }
        }
    }

    static class Encode {
        public static void main (String[] argv) throws Exception {
            byte[] data = new byte[4];
            for (String a : argv) {
                int x = Integer.parseInt(a);
                data[3] = (byte)((x & 0x01ffffff)>>>24);
                data[2] = (byte)((x & 0x00ffffff)>>>16);
                data[1] = (byte)((x & 0x0000ffff)>>>8);
                data[0] = (byte)(x & 0xff);
                System.out.println(x+": "+encode25 (data));
            }
        }
    }

    static class TestRandom {
        public static void main (String[] argv) throws Exception {
            if (argv.length < 3) {
                System.out.println
                    ("Usage: Base32$TestRandom SIZE BITS MASK [NTHREADS=1]");
                System.out.println("where SIZE is power of 2 (i.e., 1<<SIZE)");
                System.out.println
                    ("BITS is the number of bits (20, 25, 30, ...)");
                System.out.println("MASK must be in hex");
                System.exit(1);
            }

            long size = 1l<<Integer.parseInt(argv[0]);
            int bits = Integer.parseInt(argv[1]);
            int mask = Integer.parseInt(argv[2], 16);
            int threads = 1;
            if (argv.length > 3) {
                threads = Integer.parseInt(argv[3]);
            }
            logger.info("## size="+size+" bits="+bits
                        +" mask="+mask+" threads="+threads);
            if (threads < 2) {
                testRandom (size, bits, mask);
            }
            else {
                ExecutorService es = Executors.newFixedThreadPool(threads);
                try {
                    Future[] futures = new Future[threads];
                    for (int i = 0; i < threads; ++i) {
                        futures[i] = es.submit
                            (new TestRandomWorker (size, bits, mask));
                    }

                    for (Future f : futures) {
                        f.get();
                    }
                }
                finally {
                    es.shutdown();
                }
            }
        }
    }

    public static void main (String[] argv) throws Exception {
        //test20 ();
        //testUnrank ();
        //test40 (1l<<31);
        //testAll (1<<30);
        //test (1l<<40, 1l<<41, 45, 0x1f);

        if (argv.length < 4) {
            System.out.println
                ("Usage: Base32 START END BITS MASK [NTHREADS=1]");
            System.out.println("where START and END are bit positions!");
            System.out.println("BITS is the number of bits (20, 25, 30, ...)");
            System.out.println("MASK must be in hex");
            System.exit(1);
        }

        long start = Long.parseLong(argv[0]);
        if (start < 64) {
            start = 1l<<start;
        }

        long end = Long.parseLong(argv[1]);
        if (end < 64) {
            end = 1l<<end;
        }

        if (end < start) {
            end += start;
        }

        int bits = Integer.parseInt(argv[2]);
        int mask = Integer.parseInt(argv[3], 16);

        int threads = 1;
        if (argv.length > 4) {
            threads = Integer.parseInt(argv[4]);
        }

        logger.info("## start="+start+" end="+end+" bits="+bits
                    +" mask="+mask+" threads="+threads);
        if (threads < 2) {
            test (start, end, bits, mask);
        }
        else {
            ExecutorService es = Executors.newFixedThreadPool(threads);
            try {
                Future[] futures = new Future[threads];
                long range = (end-start)/threads;
                for (int i = 0; i < threads; ++i) {
                    end = start+range;
                    futures[i] = es.submit
                        (new TestWorker (start, end, bits, mask));
                    start = end;
                }
                
                long total = 0;
                for (Future f : futures) {
                    Long c = (Long)f.get();
                    total += c;
                }
                logger.info("## "+total+" tests performed!");
            }
            finally {
                es.shutdown();
            }
        }
    }
}
