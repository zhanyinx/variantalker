# Taken from http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
import os


def count_lines(filename):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b:
                break
            yield b

    with open(filename, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count
