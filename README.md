ngstools
=============

# fastq parser
Eficent fastq and other parsers for NGS analysis.  
Short example: 

```python 
inputfile="example.fastq"
barcodes = {
	"TGAACG" : 1,
	"ACCAGA" : 2,
	"TATGCC" : 3,
	"TCTGGA" : 4
}

outputdir = results 

five_const = "GGGCAACTCC"
three_const = "AAAATGGCTA"

biotools.MakeDirs(outputdir)

df1, quality = biotools.fastryderlist(inputfile,barcodes,1,five_const,three_const,"test","ntest",False)

```


## Dependences
regex 
numpy 
pandas 
biopython 
