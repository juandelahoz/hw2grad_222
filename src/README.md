# sequence aligner.
# to run:

index:
$ python variants_discoverer.py ${reference}

align:
$ python variants_discoverer.py ${reference} ${reads}

sort:
$ sort -k4 -n ${aligned_reads} > ${aligned_reads}.srt

pileup:
$ python variants_discoverer.py ${reference} ${aligned_sorted_reads in the right format!} 
