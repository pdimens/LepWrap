rule validation:
  input:
  output:
  message:
  threads: 30
  shell:
    """
    awk -f LA/scripts/liftover.awk chr1.agp order1.input | sort -V | grep CHR > order1.liftover
    """
