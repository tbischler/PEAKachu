import argparse


def frac_float(string_value):
    value = float(string_value)
    if value < 0.0 or value > 1.0:
        raise argparse.ArgumentTypeError("{} not in range [0.0, 1.0]".format(
            value))
    return value
