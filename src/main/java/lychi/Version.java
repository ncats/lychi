package lychi;
public final class Version {
   public static final String VERSION = "20180207";    
   public static final String COMMIT = "fe2ea2a";
   public static final String USER = "tyler";
   public static final String TIMESTAMP = "02/07/2018 at 17:57:34 EST";
   private Version () {}
   public static void main (String[] argv) {
      System.out.println("LyChI Version: commit="+COMMIT+" version="+VERSION+" timestamp="+TIMESTAMP+" user="+USER);
   }
}
