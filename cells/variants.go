package cells

import (
	"github.com/vertgenlab/gonomics/dna"
)

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

type CellVar struct {
	Vid             int // variant ID, equivalent to Variant.Id
	Genotype        Zygosity
	GenotypeQuality int     // GQ
	ReadDepth       int     // DP
	AltReads        int     // AD[alleleIdx]
	Af              float64 // allele frequency by read count
}

type Zygosity byte

const (
	WildType Zygosity = iota
	Heterozygous
	Homozygous
	Hemizygous
)
