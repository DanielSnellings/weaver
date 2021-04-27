package loh

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"testing"
)

var v1file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz"
var v2file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/V2_Analysis/CM2001_2.vcf.gz"
var testfile = "../cells/testdata/small.vcf"

func TestFindRoh(t *testing.T) {
	d := cells.ReadVcf(v1file, cells.DefaultCellFilter, cells.DefaultGlobalFilter, cells.DefaultVcfQual)
	answer := FindAllRunsOfHomozygosity(d)

	out := fileio.EasyCreate("roh.csv")
	_, err := fmt.Fprintln(out, "Cell,Chromosome,Start,End")
	if err != nil {
		log.Panic()
	}

	for i := range answer {
		if len(answer[i]) == 0 {
			continue
		}
		for _, run := range answer[i] {
			_, err = fmt.Fprintf(out, "Cell_%d,%s,%d,%d\n", i, d.Variants[run[0]].Chr, d.Variants[run[0]].Pos+1, d.Variants[run[len(run)-1]].Pos+1)
			if err != nil {
				log.Panic()
			}
		}
	}

	err = out.Close()
	if err != nil {
		log.Panic()
	}
}
