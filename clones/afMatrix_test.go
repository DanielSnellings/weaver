package clones

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"github.com/vertgenlab/gonomics/fileio"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
	"log"
	"strings"
	"testing"
)

var defaultVcfQual float64 = 100
var defaultCellFilter = cells.CellFilterParam{MinGenotypeQuality: 30, MinGenotypeDepth: 10, MinReadAf: 0.2}
var defaultGlobalFilter = cells.GlobalFilterParam{MinGenotypedFrac: 0.5, MinGenotypesPresent: 0.5, MinCellAf: 0.01}

var v1file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz"
var v2file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/V2_Analysis/CM2001_2.vcf.gz"
var testfile = "../cells/testdata/small.vcf"

func TestPC(t *testing.T) {
	d := cells.ReadVcf(v1file, defaultCellFilter, defaultGlobalFilter, defaultVcfQual)
	afMat := generateAfMatrix(d)
	var pc stat.PC
	pc.PrincipalComponents(afMat, nil)

	fmt.Printf("variances = %.4f\n\n", pc.VarsTo(nil))

	// Project the data onto the first 2 principal components.
	k := 2
	var proj mat.Dense
	var vec mat.Dense
	pc.VectorsTo(&vec)
	proj.Mul(afMat, vec.Slice(0, len(d.Variants), 0, k))

	out := fileio.EasyCreate("pc.csv")
	_, err := fmt.Fprintln(out, getHeader(d))
	if err != nil {
		log.Panic()
	}
	var vals []float64
	rows, _ := proj.Dims()
	for i := 0; i < rows; i++ {
		vals = proj.RawRowView(i)
		_, err = fmt.Fprintf(out, "%.4f,%.4f%s\n", vals[0], vals[1], getAfStringCSV(i, afMat))
		if err != nil {
			log.Panic()
		}
	}
	err = out.Close()
	if err != nil {
		log.Panic()
	}

	//fmt.Printf("proj = %.4f", mat.Formatted(&proj, mat.Prefix("       ")))
}

func getHeader(d *cells.Data) string {
	var answer strings.Builder
	answer.Grow(len(d.Variants) * 20)
	answer.WriteString("pc1,pc2")
	for i := range d.Variants {
		answer.WriteString(fmt.Sprintf(",%s", d.Variants[i]))
	}
	return answer.String()
}

func getAfStringCSV(i int, m *mat.Dense) string {
	var answer strings.Builder
	vals := m.RawRowView(i)
	answer.Grow(len(vals) * 5)
	for _, val := range vals {
		answer.WriteString(fmt.Sprintf(",%.4f", val))
	}
	return answer.String()
}
