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

func generateTable(filename string, outfile string, cellFilter cells.CellFilterParam, globalFilter cells.GlobalFilterParam, minVcfQual float64, delim string, genotypeAsString bool, minCellAf float64, maxCellAf float64) {
	data := cells.ReadVcf(filename, cellFilter, globalFilter, minVcfQual)
	rows := getRows(data, minCellAf, maxCellAf)
	writeTable(outfile, generateColNames(data, delim), rows, delim, genotypeAsString)
}

func generateColNames(d *cells.Data, delim string) string {
	var s strings.Builder
	s.WriteString("Chromosome" + delim + "Position" + delim + "Ref" + delim + "Alt")
	s.Grow(len(d.Cells) * (len(delim) + 5 + 5)) // 5 bytes for "Cell_" and 5 bytes for up to 99999 cells. More will be added dynamically if needed.
	for i := range d.Cells {
		s.WriteString(fmt.Sprintf("%sCell_%d", delim, d.Cells[i].Id)) // could just be i, but using the Id to be safe
	}
	return s.String()
}

func getRows(d *cells.Data, minCellAf float64, maxCellAf float64) []row {
	rows := make([]row, len(d.Variants))

	for i := range d.Variants {
		if d.Variants[i].CellAf > maxCellAf || d.Variants[i].CellAf < minCellAf {
			continue
		}
		rows[i].Chr = d.Variants[i].Chr
		rows[i].Pos = d.Variants[i].Pos+1 // back to 1-base for user
		rows[i].Ref = d.Variants[i].Ref
		rows[i].Alt = d.Variants[i].Alt
		rows[i].Genotypes = make([]cells.Zygosity, len(d.Cells))
		for _, cellId := range d.Variants[i].CellsGenotyped {
			rows[i].Genotypes[cellId] = d.Cells[cellId].Genotypes[i].Genotype
		}
	}
	return rows
}

func writeTable(outfile string, colNames string, rows []row, delim string, genotypeAsString bool) {
	out := fileio.EasyCreate(outfile)
	var err error
	_, err = fmt.Fprintln(out, colNames)
	if err != nil {
		log.Panic(err)
	}

	for i := range rows {
		if rows[i].Chr == "" { // row was filtered out
			continue
		}
		rowToString(rows[i], delim, genotypeAsString)
		_, err := fmt.Fprintln(out, rowToString(rows[i], delim, genotypeAsString))
		if err != nil {
			log.Panic(err)
		}
	}

	err = out.Close()
	if err != nil {
		log.Panic(err)
	}
}

func rowToString(r row, delim string, genotypeAsString bool) string {
	var answer strings.Builder
	_, err := answer.WriteString(fmt.Sprintf("%s%s%d%s%s%s%s", r.Chr, delim, r.Pos, delim, dna.BasesToString(r.Ref), delim, dna.BasesToString(r.Alt)))
	if err != nil {
		log.Panic(err)
	}

	if genotypeAsString {
		for _, genotype := range r.Genotypes {
			answer.WriteString(delim + genotype.String())
		}
	} else {
		for _, i := range r.Genotypes {
			if i == cells.NoGenotype {
				answer.WriteString(delim + "NA")
			} else {
				answer.WriteString(delim + fmt.Sprint(genotypeToInt(i)))
			}
		}
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

var v1file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/CM2001_R/CM2001.vcf.gz"
var v2file = "/Users/danielsnellings/Desktop/Data/21-04-06_Tapestri_Run/V2_Analysis/CM2001_2.vcf.gz"

func main() {
	generateTable(v1file, "CM2001.csv", defaultCellFilter, defaultGlobalFilter, defaultVcfQual, ",", false, .2, .8)
}
