package plotting

import (
	"os"
	"os/exec"
	"testing"
)

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
