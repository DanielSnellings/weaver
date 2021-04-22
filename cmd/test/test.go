package main

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"github.com/vertgenlab/gonomics/dna"
)

var defaultVcfQual float64 = 100
var defaultCellFilter = cells.CellFilterParam{MinGenotypeQuality: 30, MinGenotypeDepth: 10, MinReadAf: 0.2}
var defaultGlobalFilter = cells.GlobalFilterParam{MinGenotypedFrac: 0.5, MinGenotypesPresent: 0.5, MinCellAf: 0.01}

func main() {
	data := cells.ReadVcf("/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz", defaultCellFilter, defaultGlobalFilter, defaultVcfQual)
	fmt.Println(len(data.Variants))
	fmt.Println(len(data.Cells))

	for _, variant := range data.Variants {
		fmt.Println(variant.Chr,"\t", variant.Pos,"\t", dna.BasesToString(variant.Ref),"\t", dna.BasesToString(variant.Alt),"\t", variant.CellAf, "\t", len(variant.CellsMutated), "\t", len(variant.CellsGenotyped))
	}

	//for i := range data.Variants{
	//	if len(data.Variants[i].CellsGenotyped) != len(data.Variants[i].CellsMutated) {
	//		//fmt.Println(len(data.Variants[i].CellsGenotyped))
	//		//fmt.Println(len(data.Variants[i].CellsMutated))
	//		//fmt.Println()
	//	}
	//}
}
