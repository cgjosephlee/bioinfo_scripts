#!/usr/bin/env python3

import sys
import os
import gzip
from tqdm import tqdm

# hack tqdm to behave like progressbar
class tqdm(tqdm):
    def update_to_value(self, n):
        if n > self.last_print_n:
            self.n = n
            with self._lock:
                self.display()
            self.last_print_n = self.n

in_file = sys.argv[1]
in_size = os.path.getsize(in_file)
sum_size = 0

print(in_size)

pb = tqdm(total=in_size)
with open(in_file, 'rb') as f:
    # for line in gzip.GzipFile(fileobj=f):
    for line in f:
        print(line)
        # print(f.tell())
        pb.update_to_value(f.tell())
    print(f.tell())
pb.close()

# print(sum_size)
print('finish.')
