#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import sys

try:
    content = json.load(sys.stdin)
    print(json.dumps(content, sort_keys=True, ensure_ascii=False, indent=4))
except json.decoder.JSONDecodeError as e:
    print(e)
    sys.exit(1)

sys.exit(0)
