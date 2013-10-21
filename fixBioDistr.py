#!/usr/bin/env python

import sys, fnmatch, os, re

for root, dirnames, filenames in os.walk('SAP/Bio'):
  for filename in fnmatch.filter(filenames, '*.py'):
      name = os.path.join(root, filename)
      with open(name) as f:
          new_content = re.sub(r'from (Bio[. ])', r'from SAP.\1', f.read())
          new_content = re.sub(r'import (Bio[. ])', r'import SAP.\1', new_content)
          #new_content = content.replace('from Bio', 'from SAP.Bio').replace('import Bio', 'import SAP.Bio')
#          print new_content
          out = open(name, 'w')
          out.write(new_content)
          out.close()


