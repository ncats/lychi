package lychi;

import java.util.BitSet;

public class ElementData {
    public static final int MAX_CHARGE =  2;
    public static final int MIN_CHARGE = -2;

    // atomic number, charge, valence
    public static final int[][][] ElData;
    public static final BitSet NonMetals = new BitSet (255);
    public static final BitSet SpecialMetals = new BitSet (255);

    static {
	// charge[index]: -2 [0], -1[1], 0[2], 1[3], 2[4]
	ElData = new int[255][5][];
	int i = 0;

	// "H"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	/*
	// "D"
	ElData[][0] = new int[]{0};
	ElData[][1] = new int[]{0};
	ElData[][2] = new int[]{1};
	ElData[][3] = new int[]{0};
	ElData[][4] = new int[]{0};

	// "T"
	ElData[][0] = new int[]{0};
	ElData[][1] = new int[]{0};
	ElData[][2] = new int[]{1};
	ElData[][3] = new int[]{0};
	ElData[][4] = new int[]{0};
	*/

	// "He"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{0};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	// "Li"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Be"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{1};
	ElData[i][4] = new int[]{0};

	// "B"
	++i;
	ElData[i][0] = new int[]{3};
	ElData[i][1] = new int[]{4};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{2};
	ElData[i][4] = new int[]{1};
	NonMetals.set(i);

	// "C"
	++i;
	ElData[i][0] = new int[]{2};
	ElData[i][1] = new int[]{3};
	ElData[i][2] = new int[]{4,3,2};
	ElData[i][3] = new int[]{3};
	ElData[i][4] = new int[]{2};
	NonMetals.set(i);

	// "N"
	++i;
	ElData[i][0] = new int[]{1};
	ElData[i][1] = new int[]{2};
	ElData[i][2] = new int[]{3,5};
	ElData[i][3] = new int[]{4};
	ElData[i][4] = new int[]{3};
	NonMetals.set(i);

	// "O"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{1};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{3,5};
	ElData[i][4] = new int[]{4};
	NonMetals.set(i);

	// "F"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{2};
	ElData[i][4] = new int[]{3,5};
	NonMetals.set(i);

	// "Ne"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{0};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	// "Na"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Mg"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{1};
	ElData[i][4] = new int[]{0};

	// "Al"
	++i;
	ElData[i][0] = new int[]{3,5};
	ElData[i][1] = new int[]{4};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{2};
	ElData[i][4] = new int[]{1};

	// "Si"
	++i;
	ElData[i][0] = new int[]{2};
	ElData[i][1] = new int[]{3,};
	ElData[i][2] = new int[]{4};
	ElData[i][3] = new int[]{3};
	ElData[i][4] = new int[]{2};
	NonMetals.set(i);

	// "P"
	++i;
	ElData[i][0] = new int[]{1,3,5,7};
	ElData[i][1] = new int[]{2,4,6};
	ElData[i][2] = new int[]{3,5};
	ElData[i][3] = new int[]{4};
	ElData[i][4] = new int[]{3};
	NonMetals.set(i);

	// "S"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{1,3,5,7};
	ElData[i][2] = new int[]{2,4,6};
	ElData[i][3] = new int[]{3,5};
	ElData[i][4] = new int[]{4};
	NonMetals.set(i);

	// "Cl"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1,3,5,7};
	ElData[i][3] = new int[]{2,4,6};
	ElData[i][4] = new int[]{3,5};
	NonMetals.set(i);

	// "Ar"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{0};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	// "K"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ca"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{1};
	ElData[i][4] = new int[]{0};

	// "Sc"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ti"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "V"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,4,5};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Cr"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Mn"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,4,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Fe"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,4,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Co"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ni"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Cu"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1,2};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Zn"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ga"
	++i;
	ElData[i][0] = new int[]{3,5};
	ElData[i][1] = new int[]{4};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{1};

	// "Ge"
	++i;
	ElData[i][0] = new int[]{2,4,6};
	ElData[i][1] = new int[]{3,5};
	ElData[i][2] = new int[]{4};
	ElData[i][3] = new int[]{3};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	// "As"
	++i;
	ElData[i][0] = new int[]{1,3,5,7};
	ElData[i][1] = new int[]{2,4,6};
	ElData[i][2] = new int[]{3,5};
	ElData[i][3] = new int[]{4};
	ElData[i][4] = new int[]{3};
	NonMetals.set(i);

	// "Se"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{1,3,5,7};
	ElData[i][2] = new int[]{2,4,6};
	ElData[i][3] = new int[]{3,5};
	ElData[i][4] = new int[]{4};
	NonMetals.set(i);

	// "Br"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1,3,5,7};
	ElData[i][3] = new int[]{2,4,6};
	ElData[i][4] = new int[]{3,5};
	NonMetals.set(i);

	// "Kr"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{0};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	// "Rb"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Sr"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{1};
	ElData[i][4] = new int[]{0};

	// "Y"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Zr"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Nb"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,5};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Mo"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4,5,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Tc"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{7};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ru"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,4,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Rh"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Pd"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ag"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Cd"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "In"
	++i;
	ElData[i][0] = new int[]{3,5};
	ElData[i][1] = new int[]{2,4};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{1};

	// "Sn"
	++i;
	ElData[i][0] = new int[]{2,4,6};
	ElData[i][1] = new int[]{3,};
	ElData[i][2] = new int[]{2,4};
	ElData[i][3] = new int[]{3};
	ElData[i][4] = new int[]{0};

	// "Sb"
	++i;
	ElData[i][0] = new int[]{1,3,5,7};
	ElData[i][1] = new int[]{2,4,6};
	ElData[i][2] = new int[]{3,5};
	ElData[i][3] = new int[]{2,4};
	ElData[i][4] = new int[]{3};
	SpecialMetals.set(i);

	// "Te"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{1,3,5,7};
	ElData[i][2] = new int[]{2,4,6};
	ElData[i][3] = new int[]{3,5};
	ElData[i][4] = new int[]{2,4};
	NonMetals.set(i);

	// "I"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1,3,5,7};
	ElData[i][3] = new int[]{2,4,};
	ElData[i][4] = new int[]{3,5};
	NonMetals.set(i);

	// "Xe"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{0};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	// "Cs"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ba"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{1};
	ElData[i][4] = new int[]{0};

	// "La"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ce"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Pr"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Nd"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Pm"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Sm"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Eu"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Gd"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Tb"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Dy"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ho"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Er"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Tm"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Yb"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Lu"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Hf"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ta"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{5};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "W"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4,5,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Re"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,4,6,7};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Os"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,4,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ir"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,3,4,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Pt"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Au"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Hg"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1,2};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	SpecialMetals.set(i);

	// "Tl"
	++i;
	ElData[i][0] = new int[]{3,5};
	ElData[i][1] = new int[]{2,4};
	ElData[i][2] = new int[]{1,3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Pb"
	++i;
	ElData[i][0] = new int[]{2,4,6};
	ElData[i][1] = new int[]{3,};
	ElData[i][2] = new int[]{2,4};
	ElData[i][3] = new int[]{3};
	ElData[i][4] = new int[]{0};

	// "Bi"
	++i;
	ElData[i][0] = new int[]{1,3,5,7};
	ElData[i][1] = new int[]{2,4,6};
	ElData[i][2] = new int[]{3,5};
	ElData[i][3] = new int[]{2,4};
	ElData[i][4] = new int[]{3};

	// "Po"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{1,3,5,7};
	ElData[i][2] = new int[]{2,4,6};
	ElData[i][3] = new int[]{3,5};
	ElData[i][4] = new int[]{2,4};

	// "At"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1,3,5,7};
	ElData[i][3] = new int[]{2,4,};
	ElData[i][4] = new int[]{3,5};
	NonMetals.set(i);

	// "Rn"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{0};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
	NonMetals.set(i);

	// "Fr"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Ra"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{2};
	ElData[i][3] = new int[]{1};
	ElData[i][4] = new int[]{0};

	// "Ac"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Th"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Pa"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4,5};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "U"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4,5,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Np"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4,5,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Pu"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4,5,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Am"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4,5,6};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Cm"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Bk"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3,4};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Cf"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Es"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Fm"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Md"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{3};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "No"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Lr"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};

	// "Rf"
	++i;
	ElData[i][0] = new int[]{0};
	ElData[i][1] = new int[]{0};
	ElData[i][2] = new int[]{1};
	ElData[i][3] = new int[]{0};
	ElData[i][4] = new int[]{0};
    }

    public static boolean checkValence (int atno, int charge, 
					int valence, int hcount) {
	int[][] cv = ElData[atno];

	if (cv != null && charge >= MIN_CHARGE && charge <= MAX_CHARGE) {
	    int[] v = cv[charge-MIN_CHARGE];
	    int av = valence - hcount;
	    for (int i = 0; i < v.length; ++i) {
		if (v[i] == av) {
		    return true;
		}
	    }
	}
	return false;
    }

    public static boolean isMetal (int atno) {
	return !NonMetals.get(atno);
    }

    public static boolean isSpecialMetal (int atno) {
	return SpecialMetals.get(atno);
    }

    public static int getChargeForValence 
	(int atno, int minCharge, int valence) {
	int[][]cv = ElData[atno];
	if (minCharge < MIN_CHARGE) {
	    minCharge = MIN_CHARGE;
	}
	for (int c = minCharge; c <= MAX_CHARGE; ++c) {
	    int[] v = cv[c - MIN_CHARGE];
	    for (int i = 0; i < v.length; ++i) {
		if (v[i] == valence) {
		    return c;
		}
	    }
	}
	return Integer.MAX_VALUE;
    }
}
