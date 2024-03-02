#!/usr/bin/env python3 
# Plot absorption spectrum

# Import local files 
import common.cross as cross
import common.dace as dace
import common.utils as utils

import os, argparse

# Main function
def main(formula:str, source:str, target_p:str, target_t:str, yunits:str):

    source = utils.sourcesafe(source)

    formula_path =  os.path.join(utils.dirs[source], formula.strip())

    match source:
        case "dace":
            close_path = dace.find_bin_close(formula_path, float(target_p), float(target_t))
        case "hitran":
            close_path = ""
        case "exomol":
            close_path = ""

    yunits = yunits.strip().lower()
    match yunits:
        case "cm2g-1":          yunits_int=0
        case "cm2molecule-1":   yunits_int=1
        case "m2kg-1":          yunits_int=2
        case _:
            raise Exception("Invalid units [%s]"%yunits)

    xc = cross.xsec(formula, source, close_path)
    xc.read()
    xc.plot(yunits=yunits_int)


# Run main function
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Plot absorption spectrum')
    parser.add_argument('absorber', type=str, help='Absorber name')
    parser.add_argument('source',   type=str, help='Source database name')
    parser.add_argument('pres',     type=str, help='Target pressure [bar]')
    parser.add_argument('temp',     type=str, help='Target temperature [K]')
    parser.add_argument('--yunits', type=str, default="cm2g-1", help='y-axis units')

    args = parser.parse_args()
    
    
    main(args.absorber,   # absorber
         args.source,     # database
         args.pres,       # target pressure [bar]
         args.temp,       # target temperature [K],
         args.yunits       # y-axis units
         )


# End of file