package lychi;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import chemaxon.formats.MolFormatException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
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
		
		public Molecule getMolecule() throws MolFormatException{
			MolHandler mh = new MolHandler();
			mh.setMolecule(input);
			return mh.getMolecule();
			
		}
	}
	

	private LychiTestInstance spec;

	public LychiRegressionTest(String ignored, LychiTestInstance spec){
		this.spec = spec;
	}
	
	public static void basicTest(Molecule m, String expected) throws Exception{
		LyChIStandardizer std = new LyChIStandardizer();
		std.standardize(m);
		String fullKey=LyChIStandardizer.hashKey(m);
		assertEquals(expected,fullKey);
	}
	
	public static Molecule shuffleMolecule(Molecule m, int[] map){
		MolAtom[] mas=m.getAtomArray();
		MolBond[] mbs=m.getBondArray();
		
		Molecule m2=new Molecule();
		
		int[] rmap = new int[map.length];
		
		for(int i=0;i<map.length;i++){
			m2.add((MolAtom)mas[map[i]].clone());
			rmap[map[i]]=i;
		}
		MolAtom[] nmas = m2.getAtomArray();
		for(int j=0;j<mbs.length;j++){
			MolBond mb = mbs[j];
			int oi1=m.indexOf(mb.getAtom1());
			int oi2=m.indexOf(mb.getAtom2());
			int ni1=rmap[oi1];
			int ni2=rmap[oi2];
			MolBond nmb=mb.cloneBond(nmas[ni1], nmas[ni2]);
			m2.add(nmb);
		}
		
		return m2;
		
	}
	
	@Test
	public void correctLychiFirstTime() throws Exception{
		basicTest(spec.getMolecule(),spec.expectedLychi);
	}
	
	@Test
	public void correctLychiAfter10Times() throws Exception{
		for (int i=0;i<10;i++){
			correctLychiFirstTime();
		}
	}
	
	@Test
	public void correctLychiAfterRandomShuffle() throws Exception{
		Molecule m=spec.getMolecule();
		m.clean(2, null);
		Random r= new Random(0xFAB5EED);
		for (int i=0;i<10;i++){
			
			List<Integer> iatoms=IntStream.range(0, m.getAtomCount()).mapToObj(i2->i2).collect(Collectors.toList());
			Collections.shuffle(iatoms, r);
			int[] map =iatoms.stream().mapToInt(i1->i1).toArray();
			Molecule s=shuffleMolecule(m,map);
			
			basicTest(s,spec.expectedLychi);
		}
	}
	
	@Test
	public void daisyChainLychiAfter10Times() throws Exception{
		Molecule m=spec.getMolecule();
		m.clean(2, null);
		for (int i=0;i<10;i++){
			
			List<Integer> iatoms=IntStream.range(0, m.getAtomCount()).mapToObj(i2->i2).collect(Collectors.toList());
			Collections.shuffle(iatoms);
			int[] map =iatoms.stream().mapToInt(i1->i1).toArray();
			Molecule s=shuffleMolecule(m,map);
			basicTest(s,spec.expectedLychi);
			m=s;
		}
	}
	
	
	@Parameterized.Parameters(name = "{0}")
	public static List<Object[]> data(){
		List<LychiTestInstance> tests = new ArrayList<>();
		
		tests.add(LychiTestInstance.of("CCCCCC","U28WVSD82-2YWKXKS36P-2PQPUKLGUWW-2PWK7BQNJSU6").name("simple carbon chain"));
		tests.add(LychiTestInstance.of("O=C(O[C@H]1C[C@H]2C[C@H]3C[C@@H](C1)N2CC3=O)C4=CNC5=C4C=CC=C5","38C4U16JU-UC5KDUPMVH-UHFJLJL661C-UHCRHDK74DXU").name("cage-like structure"));
		tests.add(LychiTestInstance.of("C[C@@H]1CC[C@@H](C)CC1","T75RBW5S8-8D9T563A7Y-8YC8NQXD9W5-8Y5MFVTVS3J3").name("trans across ring"));
		tests.add(LychiTestInstance.of("C[C@H]1CC[C@@H](C)CC1","T75RBW5S8-8D9T563A7Y-8YC8NQXD9W5-8Y5JH5RWXRLR").name("cis across ring"));
		
		
		return tests.stream().map(ls->ls.asJunitInput()).collect(Collectors.toList());
	}
}
