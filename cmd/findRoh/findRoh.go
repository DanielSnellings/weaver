package main

import (
	"flag"
	"fmt"
	"github.com/ddsnellings/weaver/cells"
	"github.com/ddsnellings/weaver/loh"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"strings"
)

func usage() {
	fmt.Print(
		"findRoh - Find runs of homozygosity.\n\n" +
			"Usage:\n" +
			"  findRoh -i infile.vcf.gz\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func findRoh(infile string, outfile string, varfile string, minRunLength int, minCounts int) {
	d := cells.ReadVcf(infile, cells.DefaultCellFilter, cells.DefaultGlobalFilter, cells.DefaultVcfQual)
	roh := loh.FindAllRunsOfHomozygosity(d, minRunLength)
	counts := loh.CountRohHaplotypes(roh, d)

	var err error
	outRoh, err := os.Create(outfile)
	exception.PanicOnErr(err)
	defer outRoh.Close()
	outVar, err := os.Create(varfile)
	exception.PanicOnErr(err)
	defer outRoh.Close()

	_, err = fmt.Fprintln(outRoh, "Chr,Start,End,Length,Variants,Zygosity,Count")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintln(outVar, "Id,Chr,Pos,Ref,Alt")
	exception.PanicOnErr(err)

	for key, val := range counts {
		for i := range val.Haplotypes {
			if val.HaplotypeCounts[i] < minCounts {
				continue
			}
			_, err = fmt.Fprintf(outRoh, "%s,%d,%d,%d,%s,%d\n", key.Chr, key.Start, key.End, key.End - key.Start, val.Haplotypes[i], val.HaplotypeCounts[i])
			exception.PanicOnErr(err)
		}
	}
	for i := range d.Variants {
		_, err = fmt.Fprintf(outVar, "%d,%s,%d,%s,%s\n", i,
			d.Variants[i].Chr,
			d.Variants[i].Pos,
			getBaseString(d.Variants[i].Ref),
			getBaseString(d.Variants[i].Alt))
		exception.PanicOnErr(err)
	}
}

func getBaseString(b []dna.Base) string {
	s := dna.BasesToString(b)
	if s == "" {
		s = "-"
	}
	return s
}

func main() {
	var minRunLength *int = flag.Int("minRunLength", 5, "Minimum number of adjacent homozygous SNPs for output")
	var minCounts *int = flag.Int("minCounts", 2, "Minimum number of cells with run to qualify for output")
	var infile *string = flag.String("i", "", "Input vcf file (may be vcf.gz)")
	var outfile *string = flag.String("o", "infile.roh.csv", "Output roh file")
	var varfile *string = flag.String("v", "infile.var.csv", "Output variant ID file")
	flag.Parse()

	if *infile == "" {
		usage()
		return
	}

	if *outfile == "infile.roh.csv" {
		*outfile = strings.TrimSuffix(strings.TrimSuffix(*infile, ".gz"), ".vcf") + ".roh.csv"
	}

	if *varfile == "infile.var.csv" {
		*varfile = strings.TrimSuffix(strings.TrimSuffix(*infile, ".gz"), ".vcf") + ".var.csv"
	}

	findRoh(*infile, *outfile, *varfile, *minRunLength, *minCounts)
}
