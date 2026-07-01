import gzip
import sys
path = sys.argv[1]
with gzip.open(path, "rb") as fh:
    while fh.read(1024 * 1024):
        pass
