#!/usr/bin/env python

import sys

def describe_module_components(module, include_dunders):
    for e in sorted(dir(module)):
        if include_dunders or not e.startswith("__"):
            print("%-30s %s" % (e, type(module.__getattribute__(e))))

if __name__ == "__main__":
    include_dunders = False
    names = sys.argv[1:]
    if len(names) >= 1 and names[0] == "-a":
        include_dunders = True
        names = names[1:]
    for name in names:
        module = __import__(name)
        describe_module_components(module, include_dunders)
