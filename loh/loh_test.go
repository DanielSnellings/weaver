package loh

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"testing"
)

var infile = "/Users/danielsnellings/Desktop/21-10-14_Tapestri_Run/Hu150.cells.vcf.gz"
var v1file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz"
var v2file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/V2_Analysis/CM2001_2.vcf.gz"
var testfile = "../cells/testdata/small.vcf"

func TestFindRoh(t *testing.T) {
	d := cells.ReadVcf(infile, cells.DefaultCellFilter, cells.DefaultGlobalFilter, cells.DefaultVcfQual)

	//roh2 := FindAllRunsOfHomozygosity(d, 2)
	//roh3 := FindAllRunsOfHomozygosity(d, 3)
	//roh4 := FindAllRunsOfHomozygosity(d, 4)
	roh5 := FindAllRunsOfHomozygosity(d, 5)

	//counts2 := CountRohHaplotypes(roh2, d)
	//counts3 := CountRohHaplotypes(roh3, d)
	//counts4 := CountRohHaplotypes(roh4, d)
	counts5 := CountRohHaplotypes(roh5, d)

	fmt.Println("Chr\tStart\tEnd\tCount_min5")
	for key, val := range counts5 {
		for i := range val.Haplotypes {
			if val.HaplotypeCounts[i] < 2 {
				continue
			}
			fmt.Printf("%s\t%d\t%d\t%s\t%d\n", key.Chr, key.Start, key.End, val.Haplotypes[i], val.HaplotypeCounts[i])
		}
	}
	for i := range d.Variants {
		fmt.Println(i, d.Variants[i])
	}

	//var curr []variants.CellVar
	//for i := range d.Cells {
	//	curr = d.Cells[i].Genotypes
	//	if curr[0].Genotype == variants.WildType &&
	//		curr[2].Genotype == variants.WildType &&
	//		curr[3].Genotype == variants.WildType &&
	//		curr[5].Genotype == variants.WildType &&
	//		curr[6].Genotype == variants.WildType &&
	//		curr[7].Genotype == variants.WildType {
	//		fmt.Println(curr)
	//	}
	//}

/*
	out := fileio.EasyCreate("roh.csv")
	_, err := fmt.Fprintln(out, "Cell,Chromosome,StartPos,EndPos,NumVariants")
	if err != nil {
		log.Panic()
	}

	for i := range roh5 {
		if len(roh5[i]) == 0 {
			continue
		}
		for _, run := range roh5[i] {
			_, err = fmt.Fprintf(out, "Cell_%d,%s,%s,%d\n", i, d.Variants[run[0]].Chr, d.Variants[run[0]].Pos+1, d.Variants[run[len(run)-1]].Pos+1, len(run))
			if err != nil {
				log.Panic()
			}
		}
	}

	err = out.Close()
	if err != nil {
		log.Panic()
	}
 */

}
