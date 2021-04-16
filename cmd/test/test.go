package main

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
)

var defaultVcfQual float64 = 100

var defaultCellFilter = cells.CellFilterParam{MinGenotypeQuality: 30, MinGenotypeDepth: 10, MinReadAf: 0.2}
var defaultGlobalFilter = cells.GlobalFilterParam{MinGenotypedFrac: 0.5, MinGenotypesPresent: 0.5, MinCellAf: 0.1}

func main() {
	data := cells.ReadVcf("/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz", defaultCellFilter, defaultGlobalFilter, defaultVcfQual)
	fmt.Println(len(data.Variants))
}
