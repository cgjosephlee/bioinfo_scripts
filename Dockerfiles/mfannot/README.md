Fix privilege problem.

```
docker build . -t mfannot_patched

docker run --rm -u $(id -u):$(id -g) -v /mnt:/mnt -v $(pwd):/project mfannot_patched mfannot -g 4 --sqn --tbl -d sequence.fasta
```
