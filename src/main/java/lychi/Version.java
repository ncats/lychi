package lychi;
public final class Version {
   public static final String VERSION = "20170222";    
   public static final String COMMIT = "fe2ea2a";
   public static final String USER = "tyler";
   public static final String TIMESTAMP = "12/14/2018 at 14:53:34 EST";
   private Version () {}
   public static void main (String[] argv) {
      System.out.println("LyChI Version: commit="+COMMIT+" version="+VERSION+" timestamp="+TIMESTAMP+" user="+USER);
   }
}
