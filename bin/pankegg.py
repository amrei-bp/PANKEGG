#!/usr/bin/env python
from __future__ import print_function

import PANKEGG
print(dir(PANKEGG))
import PANKEGG.parser
import sys

if __name__ == '__main__':
    sys.exit(PANKEGG.parser.main(arg_src=sys.argv[1:]).run())