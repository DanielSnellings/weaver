package main

import (
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type row struct {
	Chr       string
	Pos       int
	Ref       []dna.Base
	Alt       []dna.Base
	Genotypes []cells.Zygosity
}

func generateTable(filename string, outfile string, cellFilter cells.CellFilterParam, globalFilter cells.GlobalFilterParam, minVcfQual float64, delim string) {
	data := cells.ReadVcf(filename, cellFilter, globalFilter, minVcfQual)
	rows := getRows(data)
	writeTable(outfile, nil, rows, delim)
}

func getRows(d *cells.Data) []row {
	rows := make([]row, len(d.Variants))
	for i := range d.Variants {
		rows[i].Chr = d.Variants[i].Chr
		rows[i].Pos = d.Variants[i].Pos
		rows[i].Ref = d.Variants[i].Ref
		rows[i].Alt = d.Variants[i].Alt
		rows[i].Genotypes = make([]cells.Zygosity, len(d.Cells))
		for _, cellId := range d.Variants[i].CellsGenotyped {
			rows[i].Genotypes[cellId] = d.Cells[cellId].Genotypes[i].Genotype
		}
	}
	return rows
}

func writeTable(outfile string, colNames []string, rows []row, delim string) {
	out := fileio.EasyCreate(outfile)
	_, err := fmt.Fprintln(out, strings.Join(colNames, delim))
	if err != nil {
		log.Panic(err)
	}

	for i := range rows {
		_, err := fmt.Fprintln(out, rowToString(rows[i], delim))
		if err != nil {
			log.Panic(err)
		}
	}

	err = out.Close()
	if err != nil {
		log.Panic(err)
	}
}

func rowToString(r row, delim string) string {
	var answer strings.Builder
	_, err := answer.WriteString(fmt.Sprintf("%s%s%d%s%s%s%s", r.Chr, delim, r.Pos, delim, dna.BasesToString(r.Ref), delim, dna.BasesToString(r.Alt)))
	if err != nil {
		log.Panic(err)
	}
	for _, i := range r.Genotypes {
		answer.WriteByte('\t')
		answer.WriteString(fmt.Sprint(genotypeToInt(i)))
	}

	return answer.String()
}

func genotypeToInt(z cells.Zygosity) int {
	switch z {
	case cells.NoGenotype:
		return -1
	case cells.WildType:
		return 0
	case cells.Heterozygous:
		return 1
	case cells.Homozygous:
		return 2
	case cells.Hemizygous:
		return 3
	default:
		return -1
	}
}

var defaultVcfQual float64 = 100
var defaultCellFilter = cells.CellFilterParam{MinGenotypeQuality: 30, MinGenotypeDepth: 10, MinReadAf: 0.2}
var defaultGlobalFilter = cells.GlobalFilterParam{MinGenotypedFrac: 0.5, MinGenotypesPresent: 0.5, MinCellAf: 0.01}

func main() {
	generateTable("/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz", "CM2001.tsv", defaultCellFilter, defaultGlobalFilter, defaultVcfQual, "\t")
}
