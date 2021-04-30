package loh

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"testing"
)

var v1file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz"
var v2file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/V2_Analysis/CM2001_2.vcf.gz"
var testfile = "../cells/testdata/small.vcf"

func TestFindRoh(t *testing.T) {
	d := cells.ReadVcf(v1file, cells.DefaultCellFilter, cells.DefaultGlobalFilter, cells.DefaultVcfQual)

	//roh2 := FindAllRunsOfHomozygosity(d, 2)
	//roh3 := FindAllRunsOfHomozygosity(d, 3)
	//roh4 := FindAllRunsOfHomozygosity(d, 4)
	roh5 := FindAllRunsOfHomozygosity(d, 5)

	//var found bool
	//for i := range roh5 {
	//	for _, j := range roh5[i] {
	//		if j.Region(d.Variants).End == 169237430 {
	//			for _, id := range j {
	//				fmt.Println(d.Variants[id])
	//			}
	//			found = true
	//			break
	//		}
	//	}
	//	if found {
	//		break
	//	}
	//}

	//counts2 := CountRohHaplotypes(roh2, d)
	//counts3 := CountRohHaplotypes(roh3, d)
	//counts4 := CountRohHaplotypes(roh4, d)
	counts5 := CountRohHaplotypes(roh5, d)

	fmt.Println("Chr\tStart\tEnd\tCount_min5")
	for key, val := range counts5 {
		for i := range val.Haplotypes {
			fmt.Printf("%s\t%d\t%d\t%s\t%d\n", key.Chr, key.Start, key.End, val.Haplotypes[i], val.HaplotypeCounts[i])
		}
	}

	/*
	out := fileio.EasyCreate("roh.csv")
	_, err := fmt.Fprintln(out, "Cell,Chromosome,StartOrEnd,Pos")
	if err != nil {
		log.Panic()
	}

	for i := range roh5 {
		if len(roh5[i]) == 0 {
			continue
		}
		for _, run := range roh5[i] {
			_, err = fmt.Fprintf(out, "Cell_%d,%s,%s,%d\n", i, d.Variants[run[0]].Chr, "start", d.Variants[run[0]].Pos+1)
			if err != nil {
				log.Panic()
			}
			_, err = fmt.Fprintf(out, "Cell_%d,%s,%s,%d\n", i, d.Variants[run[0]].Chr, "end", d.Variants[run[len(run)-1]].Pos+1)
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
