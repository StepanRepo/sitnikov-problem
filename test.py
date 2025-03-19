#! venv/bin/python3

size = 8


print(f"THREAD {0}: PRINTING")

for i in range(size - 1):
    print(f"THREAD {i+1}: 0.0 %")

source = 7
nlines = size - source

for i in range(10):
    print(f"\033[{nlines}A", end="")
    print(f"\rTHREAD {source}: {i}.0 %")

    if nlines > 1:
        print(f"\033[{nlines-1}B", end="")

print(f"\r here {source} {nlines}          ", end = "")
