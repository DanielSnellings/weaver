package variants

import (
	"github.com/vertgenlab/gonomics/dna"
	"sort"
)

// FindHeterozygous finds all variants in the input that are likely to be
// heterozygous in constitutional DNA (>20% AF && < 80% AF).
// Returns a slice of Variant Ids corresponding to heterozygous variants.
func FindHeterozygous(v []Variant) []int {
	return FindVariantsInAfRange(v, 0.2, 0.8)
}

// FindVariantsInAfRange finds all variants with a cell allele frequency
// (pseudobulk) greater than the input minAf and less than the input maxAf.
func FindVariantsInAfRange(v []Variant, minAf float64, maxAf float64) []int {
	answer := make([]int, 0, len(v))
	for i := range v {
		if v[i].CellAf > minAf && v[i].CellAf < maxAf {
			answer = append(answer, v[i].Id)
		}
	}
	return answer
}

// SortByCoord sorts a slice of variants by genomic coordinate
// Ties by position are broken by alt allele length, then
// lexicographically by alt sequence.
func SortByCoord(v []Variant) {
	sort.Slice(v, func(i, j int) bool {
		return lessVariants(v[i], v[j])
	})
}

// SortIdsByCoord sorts an input []VariantIds (ints) by genomic coordinate.
// Does not alter the input []Variant ordering.
func SortIdsByCoord(ids []int, v []Variant) {
	vCopy := make([]Variant, len(v))
	copy(vCopy, v)
	sort.Slice(ids, func(i, j int) bool {
		return lessVariants(v[ids[i]], v[ids[j]])
	})
}

// lessVariants returns whether a should be before b in a sorted slice.
func lessVariants(a, b Variant) bool {
	switch {
	case a.Chr < b.Chr:
		return true
	case a.Chr > b.Chr:
		return false
	case a.Pos < b.Pos:
		return true
	case a.Pos > b.Pos:
		return false
	case len(a.Alt) < len(b.Alt):
		return true
	case len(a.Alt) > len(b.Alt):
		return false
	default:
		return dna.CompareSeqsIgnoreCase(a.Alt, b.Alt) < 0
	}
}
