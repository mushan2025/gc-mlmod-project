import sys
path = sys.argv[1]
with open(path, "rt", encoding="utf-8", errors="strict") as fh:
    for _ in range(5):
        if fh.readline() == "":
            break
