printf "n\ttime\n"
for n in 40000 100000 200000 300000 1000000 2000000 3000000 4000000; do
  samtools view -h /ssd/research/rsem-stan/rsem_1.transcript.bam | head -$n|samtools view -b > example.bam
  t=$(/bin/time -p ./conn_comps example.bam 2>&1 >/dev/null | grep real | cut -d' ' -f2)
  printf "%d\t%f\n" "$n" "$t"
done
