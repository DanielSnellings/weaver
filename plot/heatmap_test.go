package plot

import (
	"os"
	"os/exec"
	"testing"
)

// TRY OUT gonum/plot FOR PLOTTING AS WELL AS GGPLOT2
func TestHeatmap(t *testing.T) {
	cmd := exec.Command("./heatmap.R")

	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	err := cmd.Run()
	if err != nil {
		t.Errorf(err.Error())
	}
}
