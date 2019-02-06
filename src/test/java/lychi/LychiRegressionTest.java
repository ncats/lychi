package lychi;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;

@RunWith(Parameterized.class)
public class LychiRegressionTest {
	
	public static class LychiTestInstance{
		String name;
		String input;
		String expectedLychi;
		
		
		public static LychiTestInstance of(String smi, String lychi){
			LychiTestInstance ltest= new LychiTestInstance();
			ltest.input=smi;
			ltest.expectedLychi=lychi;
			ltest.name=smi;
			return ltest;
		}
		
		public LychiTestInstance name(String n){
			this.name=n;
			return this;
		}
		
		
		public Object[] asJunitInput(){
			return new Object[]{this.name, this};
		}
	}
	

	private LychiTestInstance spec;

	public LychiRegressionTest(String ignored, LychiTestInstance spec){
		this.spec = spec;
	}
	
	@Test
	public void correctLychiFirstTime() throws Exception{
		MolHandler mh = new MolHandler();
		mh.setMolecule(spec.input);
		Molecule m =mh.getMolecule();
		LyChIStandardizer std = new LyChIStandardizer();
		std.standardize(m);
		String fullKey=LyChIStandardizer.hashKey(m);
		assertEquals(spec.expectedLychi,fullKey);
	}
	
	@Test
	public void correctLychiAfter10Times() throws Exception{
		for (int i=0;i<10;i++){
			correctLychiFirstTime();
		}
	}
	
	
	@Parameterized.Parameters(name = "{0}")
	public static List<Object[]> data(){
		List<LychiTestInstance> tests = new ArrayList<>();
		
		tests.add(LychiTestInstance.of("CCCCCC","U28WVSD82-2YWKXKS36P-2PQPUKLGUWW-2PWK7BQNJSU6").name("simple carbon chain"));
		tests.add(LychiTestInstance.of("O=C(O[C@H]1C[C@H]2C[C@H]3C[C@@H](C1)N2CC3=O)C4=CNC5=C4C=CC=C5","38C4U16JU-UC5KDUPMVH-UHFJLJL661C-UHCV935QVNKQ").name("simple carbon chain"));

		return tests.stream().map(ls->ls.asJunitInput()).collect(Collectors.toList());
	}
}
