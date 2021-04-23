package cells

import (
	"github.com/vertgenlab/gonomics/dna"
)

// Variant stores information about a DNA change and the Id numbers
// for each cell with a valid genotype and mutation
type Variant struct {
	Id             int // position in variant slice
	Chr            string
	Pos            int // zero base pos
	Ref            []dna.Base
	Alt            []dna.Base
	CellsGenotyped []int   // all Cell.Id with variant genotyped
	CellsMutated   []int   // all Cell.Id with variant present
	GenotypedFrac  float64 // % of post-filter cells with passing genotype
	CellAf         float64 // allele frequency by cell count //TODO split to true allele frequency (het/hom aware) and cells mutated fraction
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

const (
	NoGenotype Zygosity = iota
	WildType
	Heterozygous
	Homozygous
	Hemizygous
)
