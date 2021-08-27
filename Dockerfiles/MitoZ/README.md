Source: https://github.com/linzhi2013/MitoZ/blob/master/version_2.4-alpha/Dockerfile

Fix privilege problem, expnad --clade option.

```
docker build . -t mitoz_patched:2.4-alpha

docker run --rm -u $(id -u):$(id -g) -v /mnt:/mnt -v $(pwd):/project mitoz_patched:2.4-alpha /app/MitoZ/MitoZ.py annotate --genetic_code 4 --thread_number 12 --outprefix mitoZ --fastafile input.fasta
```
