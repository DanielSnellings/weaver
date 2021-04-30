// Package loh provides tools for identifying cells with loss of heterozygosity (LOH)
// and discriminating LOH from linked allelic dropout (LADO)
//
// The basic idea is to identify informative heterozygous (constitutional) SNPs
// and use them to identify runs of homozygosity (ROH) indicative of a lost alleles.
// Plotting the ROH on a stacked interval plot (like stacking reads in IGV) should
// allow easy visualization and help in discriminating genuine loss of heterozygosity
// from linked allelic dropout (LADO).
package loh

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"github.com/ddsnellings/weaver/variants"
)

// RunOfHomozygosity stores the variant Id for each contiguous homozygous variant
// identified in a single cell. Missing genotypes are not appended to the VariantIds
// slice, however they do not break the run.
type RunOfHomozygosity []int

// Region converts a RunOfHomozygosity into a region struct (able to be key in map).
func (r RunOfHomozygosity) Region(v []variants.Variant) variants.Region {
	return variants.Region{
		Chr: v[r[0]].Chr,
		Start: v[r[0]].Pos,
		End: v[r[len(r)-1]].Pos + 1,
	}
}

// Haplotype stores the variant and zygosity information for a set of variants
type Haplotype struct {
	VariantIds []int
	Genotypes []variants.Zygosity
}

func (h Haplotype) String() string {
	return variantIdsString(h.VariantIds) + ":" + genotypesToString(h.Genotypes)
}

func variantIdsString(ids []int) string {
	if len(ids) == 0 {
		return ""
	}
	answer := fmt.Sprintf("%d", ids[0])
	for i := 1; i < len(ids); i++ {
		answer += fmt.Sprintf("-%d", ids[i])
	}
	return answer
}

func genotypesToString(z []variants.Zygosity) string {
	if len(z) == 0 {
		return ""
	}
	answer := fmt.Sprint(z[0])
	for i := 1; i < len(z); i++ {
		answer += fmt.Sprint("-", z[i])
	}
	return answer
}

// RohHaplotypes stores the unique haplotypes present in a cells with a shared region.
// e.g. Haplotype1: WT-WT-WT Haplotype2: Hom-Hom-Hom Haplotype3: Hom-WT-Hom Haplotype4: Hom-NA-Hom
// A count of cells is stored such that HaplotypeCounts[i] is the number of cells with Haplotype[i].
type RohHaplotypes struct {
	Haplotypes			 	[]Haplotype
	HaplotypeCounts 		[]int
}

// FindAllRunsOfHomozygosity identifies constitutionally heterozygous variants go to
// homozygosity across a contiguous genomic span. Returns a slice of slices of
// RunOfHomozygosity where return[i] is a slice of all homozygous runs for cell
// with Id == i. A putative ROH is only considered after is defined by > minVars
// of variants. i.e. an ROH defined by 2 SNPs is not returned if minVars == 3
func FindAllRunsOfHomozygosity(d *cells.Data, minVars int) [][]RunOfHomozygosity {
	answer := make([][]RunOfHomozygosity, len(d.Cells))
	hetVariantIds := variants.FindHeterozygous(d.Variants)
	variants.SortIdsByCoord(hetVariantIds, d.Variants)

	for i := range d.Cells {
		answer[i] = FindRunsOfHomozygosity(d.Cells[i], hetVariantIds, d.Variants, minVars)
	}
	return answer
}

// FindRunsOfHomozygosity identifies constitutionally heterozygous variants go to
// homozygosity across a contiguous genomic span in the input cell.
// The hetVariantIds input should be a slice of constitutional heterozygous variants
// sorted by genomic coordinate with variants.SortIdsByCoord. A putative ROH is only
// considered after is defined by > minVars of variants.
// i.e. an ROH defined by 2 SNPs is not returned if minVars == 3
func FindRunsOfHomozygosity(c cells.Cell, hetVariantIds []int, vars []variants.Variant, minVars int) []RunOfHomozygosity {
	var answer []RunOfHomozygosity
	var currRun RunOfHomozygosity
	for i := range hetVariantIds {

		// do not let ROH span different chromosomes
		if i > 0 && vars[hetVariantIds[i]].Chr != vars[hetVariantIds[i-1]].Chr {
			if len(currRun) >= minVars {
				answer = append(answer, currRun)
			}
			currRun = nil
		}

		switch c.Genotypes[hetVariantIds[i]].Genotype {
		case variants.Heterozygous:
			if len(currRun) >= minVars {
				answer = append(answer, currRun)
			}
			currRun = nil
			continue

		case variants.Homozygous, variants.WildType:
			currRun = append(currRun, hetVariantIds[i])

		case variants.NoGenotype, variants.Hemizygous:
			continue // do not break start a new run, but don't break an existing one
		}
	}

	if len(currRun) >= minVars { // in case we end in an ongoing run
		answer = append(answer, currRun)
	}
	return answer
}

// CountRohHaplotypes determines the number of cells harboring each unique run of homozygosity.
// The return map is keyed on the genomic region spanned by the ROH and give the number of cells with
// ROH across the keyed region.
func CountRohHaplotypes(r [][]RunOfHomozygosity, d *cells.Data) map[variants.Region]*RohHaplotypes {
	answer := make(map[variants.Region]*RohHaplotypes)
	for i := range r {
		for _, run := range r[i] {
			addHaplotype(i, run, answer, d)
		}
	}
	return answer
}

// addHaplotype determines the haplotype zygosity for the input run of homozygosity and fills the answer map accordingly.
func addHaplotype(cellId int, run RunOfHomozygosity, answer map[variants.Region]*RohHaplotypes, d *cells.Data) {
	var haplotype Haplotype
	haplotype.VariantIds = run
	haplotype.Genotypes = make([]variants.Zygosity, len(run))
	for i := range run {
		haplotype.Genotypes[i] = d.Cells[cellId].Genotypes[run[i]].Genotype
	}

	var rohHap *RohHaplotypes
	var ok bool
	if rohHap, ok = answer[run.Region(d.Variants)]; !ok {
		answer[run.Region(d.Variants)] = &RohHaplotypes{}
		rohHap = answer[run.Region(d.Variants)]
	}

	idx, found := findHaplotypeIdx(haplotype, rohHap.Haplotypes)

	if found {
		rohHap.HaplotypeCounts[idx]++
	} else {
		rohHap.Haplotypes = append(rohHap.Haplotypes, haplotype)
		rohHap.HaplotypeCounts = append(rohHap.HaplotypeCounts, 1)
	}
}

func findHaplotypeIdx(val Haplotype, searchHaps []Haplotype) (idx int, found bool) {
	for i := range searchHaps {
		if matchingHaplotype(val, searchHaps[i]) {
			return i, true
		}
	}
	return -1, false
}

// matchingHaplotype returns true if a and b have identical haplotype zygosity.
func matchingHaplotype(a, b Haplotype) bool {
	if len(a.Genotypes) != len(b.Genotypes) {
		return false
	}
	for i := range a.Genotypes {
		if a.Genotypes[i] != b.Genotypes[i] || a.VariantIds[i] != b.VariantIds[i] {
			return false
		}
	}
	return true
}


