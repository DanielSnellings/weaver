package clones

import (
	"github.com/ddsnellings/weaver/cells"
	"github.com/ddsnellings/weaver/variants"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
	"log"
)

func PrincipalComponent(d *cells.Data) stat.PC {
	var answer stat.PC
	ok := answer.PrincipalComponents(generateAfMatrix(d), nil)
	if !ok {
		log.Panic("problem with PC analysis")
	}
	return answer
}

// generateAfMatrix creates a matrix of allele frequencies where each column is variant and each row is a cell.
// The allele frequency for missing genotypes is set to -0.5.
func generateAfMatrix(d *cells.Data) *mat.Dense {
	afSlice := fetchAfFromData(d)
	return mat.NewDense(len(d.Cells), len(d.Variants), afSlice)
}

// fetchAfFromData retrieves a slice of allele frequencies for each variant in each cell. Missing genotypes
// are set to -50. The slice is ordered by cell then by variant
// e.g. {C1_VAF1, C1_VAF2, C1_VAF3, C2_VAF1, C2_VAF2, C2_VAF3, ... Cn_VAF3}
func fetchAfFromData(d *cells.Data) []float64 {
	answer := make([]float64, len(d.Cells)*len(d.Variants))
	for i := range d.Cells {
		for j := range d.Cells[i].Genotypes {
			if d.Cells[i].Genotypes[j].Genotype == variants.NoGenotype {
				answer[i*j] = -0.5
			} else {
				answer[i*j] = d.Cells[i].Genotypes[j].Af
			}
		}
	}
	return answer
}
