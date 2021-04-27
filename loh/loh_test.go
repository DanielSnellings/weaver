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
	roh := FindAllRunsOfHomozygosity(d)

	counts := CountRunsOfHomozygosity(roh, d.Variants)

	fmt.Println("Length\tCount")
	for key, val := range counts {
		fmt.Println(key.End-key.Start, "\t", val)
	}


	out := fileio.EasyCreate("roh.csv")
	_, err := fmt.Fprintln(out, "Cell,Chromosome,StartOrEnd,Pos")
	if err != nil {
		log.Panic()
	}

	for i := range roh {
		if len(roh[i]) == 0 {
			continue
		}
		for _, run := range roh[i] {
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
}
