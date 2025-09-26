"""
main, provides CLI interface for optimizing M_unit and generating table.
"""
import argparse
import sys
from runner import IPoleRunner
from table import create_munit_table

def main():
    #parse all the cli args
    parser = argparse.ArgumentParser(description="optimize m_unit to match a target flux using GRMHD simulation dumps.")
    parser.add_argument('-o', '--output', type=str, default="Munit_table.csv",
                        help="output CSV table file path.")
    parser.add_argument('-f', '--folders', type=str, nargs='+', required=True,
                        help="one or more folders with dump files.")
    parser.add_argument('-i', '--inclinations', type=float, nargs='+', default=[50.0],
                        help="observer inclination angle(s), deg.")
    parser.add_argument('--flux_goal', type=float, default=2.0,
                        help="target flux in Jansky.")
    parser.add_argument('--tolerance', type=float, default=0.05,
                        help="fractional flux diff allowed for convergence.")
    parser.add_argument('--Munit_guess', type=float, default=1e19,
                        help="initial guess for M_unit.")
    parser.add_argument('--format', choices=['Rhigh', 'criticalBeta'], default='Rhigh',
                        help="mode: Rhigh or criticalBeta.")
    parser.add_argument('--Rhigh', type=float, nargs='+', dest='Rhigh_list',
                        help="Rhigh values to try (if format=Rhigh).")
    parser.add_argument('--beta_crit_coefficient', type=float, default=0.5,
                        help="beta crit coeff (for criticalBeta only).")
    parser.add_argument('--beta_crit', type=float, default=1.0,
                        help="crit beta val (for criticalBeta only).")
    parser.add_argument('--dump_range', type=int, nargs=2, metavar=('START', 'END'), default=(1000, 2000),
                        help="dump num range (inclusive).")
    parser.add_argument('--dump_samples', type=int, default=100,
                        help="# of dumps to sample from that range.")
    parser.add_argument('--MBH', type=float, default=4.14e6,
                        help="black hole mass in M_sun.")
    parser.add_argument('--dsource', type=float, default=8.127e3,
                        help="distance to source in parsecs.")
    parser.add_argument('--kappa', type=float, default=3.5,
                        help="kappa param.")
    parser.add_argument('--emission_type', type=int, default=4,
                        help="emission model (4 = thermal synchrotron).")
    parser.add_argument('--freq_Hz', type=float, default=230e9,
                        help="obs frequency in Hz.")
    parser.add_argument('--npixel_max', type=int, default=129,
                        help="image res (pixels per side).")
    parser.add_argument('--npixel_min', type=int, default=17,
                        help="min res for adaptive.")
    parser.add_argument('--fov', type=float, default=200.0,
                        help="field of view in µas.")
    parser.add_argument('--rmax_geo', type=float, default=50.0,
                        help="max geodesic radius.")
    parser.add_argument('--phicam', type=float, default=0.0,
                        help="cam azimuth angle in deg.")
    parser.add_argument('--counterjet', type=int, choices=[0, 1], default=0,
                        help="include counterjet or not.")
    parser.add_argument('--target_nturns', type=int, default=-1,
                        help="# of half orbits for ray tracing.")
    parser.add_argument('--no_flip_retrograde', dest='flip_retrograde', action='store_false',
                        help="don't flip inclination if retrograde.")
    parser.add_argument('--overwrite', action='store_true',
                        help="overwrite entries in table.")
    parser.add_argument('--no_adaptive', dest='adaptive', action='store_false',
                        help="disable adaptive res (npixel_min = npixel_max).")
    parser.add_argument('--ipole_path', type=str, default=None,
                        help="manual path to ipole bin.")
    parser.add_argument('--verbose', action='store_true',
                        help="print ipole command.")
    parser.add_argument('--generate_images', action='store_true',
                        help="turn off '-quench' (actually make images).")

    args = parser.parse_args()

    #force fixed res if no adaptive
    if not args.adaptive:
        args.npixel_min = args.npixel_max

    #set up runner
    try:
        runner = IPoleRunner(ipole_path=args.ipole_path, generate_images=args.generate_images)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    #handle formatting mode + optional Rhigh list
    formatting = args.format
    Rhigh_list = args.Rhigh_list if formatting == 'Rhigh' else None
    if formatting != 'Rhigh' and args.Rhigh_list:
        print("Note: Rhigh values provided but format is 'criticalBeta'; Rhigh_list will be ignored.")

    #run the thing
    try:
        create_munit_table(
            runner,
            output_path=args.output,
            folders=args.folders,
            inclination_list=args.inclinations,
            dump_samples=args.dump_samples,
            dump_range=tuple(args.dump_range),
            MBH=args.MBH,
            dsource=args.dsource,
            kappa=args.kappa,
            emission_type=args.emission_type,
            flux_goal=args.flux_goal,
            fractional_tolerance=args.tolerance,
            formatting=formatting,
            Rhigh_list=Rhigh_list,
            beta_crit_coefficient=args.beta_crit_coefficient,
            beta_crit=args.beta_crit,
            flip_retrograde=args.flip_retrograde,
            overwrite=args.overwrite,
            freq_Hz=args.freq_Hz,
            npixel_max=args.npixel_max,
            npixel_min=args.npixel_min,
            fov=args.fov,
            rmax_geo=args.rmax_geo,
            phicam=args.phicam,
            counterjet=args.counterjet,
            target_nturns=args.target_nturns,
            verbose=args.verbose
        )
    except Exception as e:
        print(f"an error occurred during optimization: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()

"""
python main.py \
  -f /home/aaron/Desktop/devenv/dump_samples/MAD_a-0.5 \
  -i 50 \
  --format Rhigh \
  --Rhigh 20 \
  --dump_range 4000 4300 \
  --dump_samples 3 \
  --overwrite \
  --output Munit_table_dev.csv

"""