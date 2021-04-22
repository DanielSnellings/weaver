package cells

import (
	"fmt"
	"testing"
)

// 	testdata/small.vcf Data:
// 	Chr		Pos	Ref	Alt	Qual	Format					Cell 1																Cell 2																Cell 3
//	chr1	1	G	A	50		GT:AD:DP:GQ:PGT:PID:PL	0/0:227,0,0:100:99:.:.:0,120,1800,120,1800,1800						0/0:72,3,3:100:99:.:.:0,120,1800,120,1800,1800						0/0:163,0,0:100:99:.:.:0,120,1800,120,1800,1800
//	chr1	2	A	C,G	1000	GT:AD:DP:GQ:PL			0/0:227,0,0,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800	0/1:74,1,1,1:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800		1/2:162,1,1,1:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800
//	chr1	3	T	C,A	1500	GT:AD:DP:GQ:PL			0/0:225,2,2,2:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800	1/1:74,1,1,1:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800		0/2:163,0,0,0:100:99:0,120,1800,120,1800,1800,120,1800,1800,1800

// Parsed info
// Id			 				Zygosity		Depth			AltReads	RefReads	GenotypeQual
// Allele 1 (Var 1, Alt 1)		WT,WT,WT		100,100,100		0,1,1		100,99,99
// Allele 2 (Var 2, Alt 1)		WT,Het,Het		100,100,100		0,30,50		100,70,50
// Allele 3 (Var 2, Alt 2)		WT,WT,Het		100,100,100		0,2,40		100,98,60
// Allele 4 (Var 3, Alt 1)		WT,Hom,WT		100,100,100		2,90,0		98,10,100
// Allele 5 (Var 3, Alt 2)		WT,WT,Het		100,100,100		0,1,50		100,99,50

var defaultVcfQual float64 = 100
var defaultCellFilter = CellFilterParam{MinGenotypeQuality: 30, MinGenotypeDepth: 10, MinReadAf: 0.2}
var defaultGlobalFilter = GlobalFilterParam{MinGenotypedFrac: 0.5, MinGenotypesPresent: 0.5, MinCellAf: 0.01}

func TestReadVcf(t *testing.T) {
	data := ReadVcf("testdata/small.vcf", defaultCellFilter, defaultGlobalFilter, defaultVcfQual)
	fmt.Print(data)
}
