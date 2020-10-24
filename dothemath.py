#!/usr/bin/env python3
# Script takes a text file as argument. Prints the processed file to stdout.
# Makes use of the simpleeval module, beware of security risks. You might execute
# malicious code...
import sys
import simpleeval
import re

class _SimpleEvalAllowNumbers(simpleeval.SimpleEval):
   """Subclass SimpleEval to allow numbers to be passed to it"""
   def eval(self, expr):
      if isinstance(expr,float):
         tmp = str(expr)
      elif isinstance(expr,int):
         tmp = str(expr)
      else:
         tmp = expr
      return super(_SimpleEvalAllowNumbers, self).eval(tmp)

def main():
    assert len(sys.argv) == 2
    inFile = sys.argv[1]
    symbols = ['+', '-', '*', '/']
    math = _SimpleEvalAllowNumbers()
    with open(inFile, 'r') as f:
        for line in f.readlines():
            if any(s in line for s in symbols):
                line = re.split('(\s)', line)
                for i,el in enumerate(line):
                    if any(s in el for s in symbols):
                        try:
                            line[i] = str(math.eval(line[i]))
                        except:
                            continue
                line = "".join(line)
            print(line, end='')

if __name__ == "__main__":
    main()



