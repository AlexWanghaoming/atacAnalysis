#!/bin/bash
## bamCoverage --normalizeUsing RPGC was just used for visualization in IGV
bamCoverage --bam a.bam -o a.bw --binSize 1 --normalizeUsing RPGC\
    --effectiveGenomeSize 36696207 --extendReads
