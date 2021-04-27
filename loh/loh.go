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
	"github.com/ddsnellings/weaver/cells"
	"github.com/ddsnellings/weaver/variants"
)

// RunOfHomozygosity stores the variant Id for each contiguous homozygous variant
// identified in a single cell. Missing genotypes are not appended to the VariantIds
// slice, however they do not break the run.
type RunOfHomozygosity []int

// FindAllRunsOfHomozygosity identifies constitutionally heterozygous variants go to
// homozygosity across a contiguous genomic span. Returns a slice of slices of
// RunOfHomozygosity where return[i] is a slice of all homozygous runs for cell
// with Id == i.
func FindAllRunsOfHomozygosity(d *cells.Data) [][]RunOfHomozygosity {
	answer := make([][]RunOfHomozygosity, len(d.Cells))
	hetVariantIds := variants.FindHeterozygous(d.Variants)
	variants.SortIdsByCoord(hetVariantIds, d.Variants)

	for i := range d.Cells {
		answer[i] = FindRunsOfHomozygosity(d.Cells[i], hetVariantIds, d.Variants)
	}
	return answer
}

// FindRunsOfHomozygosity identifies constitutionally heterozygous variants go to
// homozygosity across a contiguous genomic span in the input cell.
// The hetVariantIds input should be a slice of constitutional heterozygous variants
// sorted by genomic coordinate with variants.SortIdsByCoord.
func FindRunsOfHomozygosity(c cells.Cell, hetVariantIds []int, vars []variants.Variant) []RunOfHomozygosity {
	var answer []RunOfHomozygosity
	var currRun RunOfHomozygosity
	for i := range hetVariantIds {

		// do not let ROH span different chromosomes
		if i > 0 && vars[hetVariantIds[i]].Chr != vars[hetVariantIds[i-1]].Chr {
			if len(currRun) > 1 {
				answer = append(answer, currRun)
			}
			currRun = nil
		}

		switch c.Genotypes[hetVariantIds[i]].Genotype {
		case variants.Heterozygous:
			if len(currRun) > 1 {
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

	if len(currRun) > 1 { // in case we end in an ongoing run
		answer = append(answer, currRun)
	}
	return answer
}
