package cells

// CellFilterParam defines the minimum values used for filtering cell.
// These filters are applied on a per-cell basis
type CellFilterParam struct {
	MinGenotypeQuality int     // remove genotype in cell with quality < MinGenotypeQuality // Default 30
	MinGenotypeDepth   int     // remove genotype in cell with read depth < MinGenotypeDepth // Default 10
	MinReadAf          float64 // remove genotype in cell with alternate allele frequency < MinCellAf // Default 0.2
}

func (f CellFilterParam) Passes(c CellVar) bool {
	return c.GenotypeQuality > f.MinGenotypeQuality &&
		c.ReadDepth > f.MinGenotypeDepth &&
		c.Af > f.MinReadAf
}

// GlobalFilterParam defines minimum values used for filtering cells and variants.
// These filters are applied on a global basis.
type GlobalFilterParam struct {
	MinGenotypedFrac    float64 // remove variants genotyped in < MinGenotypedFrac of cells // Default 0.5
	MinGenotypesPresent float64 // remove cells with < MinGenotypesPresent of valid genotypes present // Default 0.5
	MinCellAf           float64 // remove variants mutated in < MinCellAf of cells // Default 0.1
}

func (f GlobalFilterParam) Apply(d *Data) {
	//TODO
}
