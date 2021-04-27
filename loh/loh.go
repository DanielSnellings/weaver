// Package loh provides tools for identifying cells with loss of heterozygosity (LOH)
// and discriminating LOH from linked allelic dropout (LADO)
//
// The basic idea is to identify informative heterozygous (constitutional) SNPs
// and use them to identify runs of homozygosity (ROH) indicative of a lost alleles.
// Plotting the ROH on a stacked interval plot (like stacking reads in IGV) should
// allow easy visualization and help in discriminating genuine loss of heterozygosity
// from linked allelic dropout (LADO).
package loh

type ROH struct {
	VariantIDs []int
}
