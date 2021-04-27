package variants

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
)

// Variant stores information about a DNA change and the Id numbers
// for each cell with a valid genotype and mutation
type Variant struct {
	Id               int // position in variant slice
	Chr              string
	Pos              int // zero base pos
	Ref              []dna.Base
	Alt              []dna.Base
	CellsGenotyped   []int   // all Cell.Id with variant genotyped
	CellsMutated     []int   // all Cell.Id with variant present
	GenotypedFrac    float64 // % of post-filter cells with passing genotype
	CellsMutatedFrac float64 // fraction of genotyped cells mutated
	CellAf           float64 // allele frequency in cells. genotype aware
}

func (v Variant) String() string {
	return fmt.Sprintf("%s:%d:%s:%s", v.Chr, v.Pos, dna.BasesToString(v.Ref), dna.BasesToString(v.Alt))
}

// CellVar stores information about a paritcular variant inside a cell
type CellVar struct {
	Vid             int // variant ID, equivalent to Variant.Id
	Genotype        Zygosity
	GenotypeQuality int     // GQ
	ReadDepth       int     // DP
	AltReads        int     // AD[alleleIdx]
	Af              float64 // allele frequency by read count
}

// Zygosity of a given variant
type Zygosity byte

// String converts type Zygosity to a string. Mainly for debugging purposes
func (z Zygosity) String() string {
	switch z {
	case NoGenotype:
		return "NA"
	case WildType:
		return "WT"
	case Heterozygous:
		return "Het"
	case Homozygous:
		return "Hom"
	case Hemizygous:
		return "Hemi"
	default:
		return "NOT FOUND"
	}
}

// Region defines a genomic span. left-closed, right-open.
type Region struct {
	Chr 	string
	Start 	int // base 0
	End 	int // base 1
}

const (
	NoGenotype Zygosity = iota
	WildType
	Heterozygous
	Homozygous
	Hemizygous
)
