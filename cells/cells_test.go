package cells

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

// 	testdata/small.vcf Data:
// 	Chr		Pos	Ref	Alt	Qual	Format					Cell 1																Cell 2																Cell 3
//	chr1	1	G	A	50		GT:AD:DP:GQ:PGT:PID:PL	0/0:100,0,0:100:99:.:.:0,120,1800,120,1800,1800						0/0:99,1,0:100:99:.:.:0,120,1800,120,1800,1800						0/0:99,1,0,0,0:100:99:.:.:0,120,1800,120,1800,1800
//	chr1	2	A	C,G	1000	GT:AD:DP:GQ:PL			0/0:100,0,0,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800	0/1:70,30,2,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800	1/2:10,50,40,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800
//	chr1	3	T	C,A	1500	GT:AD:DP:GQ:PL			0/0:98,2,0,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800		1/1:9,90,1,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800		0/2:50,0,50,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800

// Parsed info
// Id			 				Zygosity		Depth			AltReads	RefReads	GenotypeQual
// Allele 1 (Var 1, Alt 1)		WT,WT,WT		100,100,100		0,1,1		100,99,99	99,99,99
// Allele 2 (Var 2, Alt 1)		WT,Het,Het		100,100,100		0,30,50		100,70,10	99,99,99
// Allele 3 (Var 2, Alt 2)		WT,WT,Het		100,100,100		0,2,40		100,70,10	99,99,99
// Allele 4 (Var 3, Alt 1)		WT,Hom,WT		100,100,100		2,90,0		98,9,50		99,99,99
// Allele 5 (Var 3, Alt 2)		WT,WT,Het		100,100,100		0,1,50		98,9,50		99,99,99

var expectedData = Data{
	Cells: []Cell{expectedCell1, expectedCell2, expectedCell3},
	Variants: []Variant{expectedAllele2, expectedAllele3, expectedAllele4, expectedAllele5},
}

// Expected Cells
var expectedCell1 = Cell{
	Id: 0,
	Genotypes: cellVar1,
	GenotypesPresent: 1,
}
var cellVar1 = []CellVar{{
	Vid: 0,
	Genotype: WildType,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 0,
	Af: 0,
}, {
	Vid: 1,
	Genotype: WildType,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 0,
	Af: 0,
}, {
	Vid: 2,
	Genotype: WildType,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 2,
	Af: float64(2)/float64(100),
}, {
	Vid: 3,
	Genotype: WildType,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 0,
	Af: 0,
}}

var expectedCell2 = Cell{
	Id: 1,
	Genotypes: cellVar2,
	GenotypesPresent: 1,
}
var cellVar2 = []CellVar{{
	Vid: 0,
	Genotype: Heterozygous,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 30,
	Af: float64(30)/float64(100),
}, {
	Vid: 1,
	Genotype: WildType,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 2,
	Af: float64(2)/float64(100),
}, {
	Vid: 2,
	Genotype: Homozygous,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 90,
	Af: float64(90)/float64(100),
}, {
	Vid: 3,
	Genotype: WildType,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 1,
	Af: float64(1)/float64(100),
}}

var expectedCell3 = Cell{
	Id: 2,
	Genotypes: cellVar3,
	GenotypesPresent: 1,
}
var cellVar3 = []CellVar{{
	Vid: 0,
	Genotype: Heterozygous,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 50,
	Af: float64(50)/float64(100),
}, {
	Vid: 1,
	Genotype: Heterozygous,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 40,
	Af: float64(40)/float64(100),
}, {
	Vid: 2,
	Genotype: WildType,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 0,
	Af: 0,
}, {
	Vid: 3,
	Genotype: Heterozygous,
	GenotypeQuality: 99,
	ReadDepth: 100,
	AltReads: 50,
	Af: float64(50)/float64(100),
}}

// Expected Variants
// Allele 1 is removed by vcf quality filter
var expectedAllele2 = Variant{
	Id: 0	,
	Chr: "chr1"	,
	Pos: 1	,
	Ref: dna.StringToBases("A"),
	Alt: dna.StringToBases("C"),
	CellsGenotyped:	[]int{0,1,2},
	CellsMutated: []int{1,2},
	GenotypedFrac: 1,
	CellAf: float64(2)/float64(3),
}
var expectedAllele3 = Variant{
	Id: 1	,
	Chr: "chr1"	,
	Pos: 1	,
	Ref: dna.StringToBases("A"),
	Alt: dna.StringToBases("G"),
	CellsGenotyped:	[]int{0,1,2},
	CellsMutated: []int{2},
	GenotypedFrac: 1,
	CellAf: float64(1)/float64(3),
}
var expectedAllele4 = Variant{
	Id: 2	,
	Chr: "chr1"	,
	Pos: 2	,
	Ref: dna.StringToBases("T"),
	Alt: dna.StringToBases("C"),
	CellsGenotyped:	[]int{0,1,2},
	CellsMutated: []int{1},
	GenotypedFrac: 1,
	CellAf: float64(1)/float64(3),
}
var expectedAllele5 = Variant{
	Id: 3	,
	Chr: "chr1"	,
	Pos: 2	,
	Ref: dna.StringToBases("T"),
	Alt: dna.StringToBases("A"),
	CellsGenotyped:	[]int{0,1,2},
	CellsMutated: []int{2},
	GenotypedFrac: 1,
	CellAf: float64(1)/float64(3),
}
var expectedNumVariants = 4


var defaultVcfQual float64 = 100
var defaultCellFilter = CellFilterParam{MinGenotypeQuality: 30, MinGenotypeDepth: 10, MinReadAf: 0.2}
var defaultGlobalFilter = GlobalFilterParam{MinGenotypedFrac: 0.5, MinGenotypesPresent: 0.5, MinCellAf: 0.01}

func TestReadVcf(t *testing.T) {
	data := ReadVcf("testdata/small.vcf", defaultCellFilter, defaultGlobalFilter, defaultVcfQual)
	if !equal(&expectedData, data) {
		t.Errorf("problem with vcf readin")
	}
}

func equal(a *Data, b *Data) bool {
	return equalCells(a.Cells, b.Cells) && equalVariants(a.Variants, b.Variants)
}

func equalCells(a []Cell, b []Cell) bool {
	if len(a) != len(b) {
		return false
	}

	for i := range a {
		switch {
		case a[i].Id != b[i].Id:
			return false
		case a[i].GenotypesPresent != b[i].GenotypesPresent:
			return false
		case !equalCellVar(a[i].Genotypes, b[i].Genotypes):
			 return false
		}
	}
	return true
}

func equalCellVar(a []CellVar, b []CellVar) bool {
	if len(a) != len(b) {
		return false
	}

	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}

	return true
}

func equalVariants(a []Variant, b []Variant) bool {
	if len(a) != len(b) {
		return false
	}

	for i := range a {
		switch {
		case a[i].Id != b[i].Id:
			return false
		case a[i].Chr != b[i].Chr:
			return false
		case a[i].Pos != b[i].Pos:
			return false
		case dna.CompareSeqsCaseSensitive(a[i].Ref, b[i].Ref) != 0:
			return false
		case dna.CompareSeqsCaseSensitive(a[i].Alt, b[i].Alt) != 0:
			return false
		case !equalInt(a[i].CellsGenotyped, b[i].CellsGenotyped):
			return false
		case !equalInt(a[i].CellsMutated, b[i].CellsMutated):
			return false
		case a[i].GenotypedFrac != b[i].GenotypedFrac:
			return false
		case a[i].CellAf != b[i].CellAf:
			return false
		}
	}
	return true
}

func equalInt(a []int, b []int) bool {
	if len(a) != len(b) {
		return false
	}

	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
