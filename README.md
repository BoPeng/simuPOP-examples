# simuPOP examples

Collection of simuPOP scripts.

## NOTE:

This repository collects many scripts that were written in Python 2.x and uses `simuOpt.Param` for parameter handling. Because newer version of simuPOP (1.1.8+) only supports Python 3.5+, and has deprecated `simuOpt.Param` (in favor of Python [argparse](https://docs.python.org/3/library/argparse.html), **many scripts will not work directly**.

Therefore,
* If you are interested only in the core simulation function, you can copy/paste the function to run it under Python 3 because the **simuPOP core stays compatible**.
* If you would like to run the scripts, please
  * Use `2to3 -w script.py` to convert the script to Python 3.
  * Manually change the parameter passing function or use [misc/simuOptParam2argparse.py](misc/simuOptParam2argparse.py) to assist the conversion process.
  * Submit a pull-request to fix the script.

## Contributions

If you have a simuPOP script that you would like to share with others, please fork this project, add your scripts, and submit a pull request. Please document your script as well as possible.
