#! /usr/bin/env python

""" 
Tools to help moving between real and reciprocal space. Calculated values
are stored in a dictionary. Units of length are angstroms (in real space) 
and 1/angstroms (in reciprocal space).
"""

import numpy
import argparse

def vector_length(vector):
    """ Takes a 3x1 array and returns a length
    """
    length =  numpy.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    length = round(length, 7)
    return length

def array_to_list(array):
    formatted_array = numpy.around(array, decimals=7)
    aslist = formatted_array.tolist()
    return aslist


def transform_space(a_in, b_in, c_in, space, scaling=1):
    """ Takes three 3x1 arrays, input space (real/reciprocal)
    and <optional> scaling. Returns three 3x1 arrays
    for transformed space (reciprocal/real) and all vector legnths
    """
    a_in_length, b_in_length, c_in_length = (
                                            vector_length(a_in), 
                                            vector_length(b_in), 
                                            vector_length(c_in))
    volume = numpy.dot(a_in, numpy.cross(b_in, c_in))
    a_out = 2*numpy.pi*numpy.cross(b_in,c_in) / volume
    b_out = 2*numpy.pi*numpy.cross(c_in,a_in) / volume
    c_out = 2*numpy.pi*numpy.cross(a_in,b_in) / volume
  
    # reciprocal lengths can also be calc'd with 2*pi/a etc
    a_out_length, b_out_length, c_out_length = (
                                                  vector_length(a_out), 
                                                  vector_length(b_out), 
                                                  vector_length(c_out))
    if space == 'real':
        real = True
    values = {'a': a_in if real else a_out,
              'b': b_in if real else b_out,
              'c': c_in if real else b_out,
              'a_star': a_out if real else a_in,
              'b_star': b_out if real else b_in,
              'c_star': c_out if real else c_in,
              'a_star_length': a_out_length if real else a_in_length,
              'b_star_length': b_out_length if real else b_in_length,
              'c_star_length': c_out_length if real else c_in_length,
              'a_length': a_in_length if real else a_out_length,
              'b_length': b_in_length if real else b_out_length,
              'c_length': c_in_length if real else c_out_length}
    return values

def calculate_distance(a_star_length,b_star_length,c_star_length, 
                       high_symmetry_point):
    """ Takes high symmetry point in reciprocol space as 
    fractional coordinates. Returns the distance to this high 
    symmetry point 
    """
    vector = (high_symmetry_point[0]*a_star_length,
              high_symmetry_point[1]*b_star_length,
              high_symmetry_point[2]*c_star_length)

    values['distance_to_hsp'] = vector_length(vector)
    return values

def calculate_number_kpoints(sampling_density, distance):
    """ Takes k point density and distance in reciprocal space.
    Returns the number of k points required to sample gamma to 
    the high symmetry point at this density.
    """
    values['number_kpoints'] = (1/sampling_density)*distance
    return values

def calculate_sampling_distance(number_kpoints, distance): 
    """Takes number of k points and distance in reciprocal space.
    Returns the distance between each k-point.
    """
    values['sampling_distance'] = distance / number_kpoints
    return values    


def print_to_screen(values):
    """ prints dictionary keys and values to screen
    """
      
    for key, value in sorted(values.items()):
        if type(value) == numpy.ndarray:
            value = array_to_list(value)
        if type(value) == float:
            value = round(value, 7)
        print (key, "=", value)

def calculate_direction(a, b):
    """ Takes two 1x3 numpy arrays. Returns the direction between
    them.
    """
    # we need to find the smallest none-zero value within a-b
    # return array with invalid entries where values are equal
    direction_masked = np.ma.masked_equal(a - b, 0)
    # fill invalid elements of array with a large number s
    direction_filled = np.ma.filled(direction_masked, 10**6)
    # return absolute values of each element 
    direction_absolute = np.absolute(direction_filled)
    smallest = np.amin(direction_absolute)

    # use the minimum absolute value as a divisor a-b
    direction = a - b / smallest
    return direction

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""A program for moving 
    between real and reciprocal space""")
    parser.add_argument( '--space', '-s', type=str, required=True,
                        help="""<required> space of input vector""", 
                        choices=['real', 'reciprocal'])
    parser.add_argument('--vectors', '-v', type=float, required=True,
                        help="""<required> lattice vectors""",
                        metavar='a1 a2 a3 b1 b2 b3 c1 c2 c3'
                        )
    parser.add_argument('-p', '--high_symmetry_point', type=float, 
                        metavar='0 0.5 0',
                        help="""<optional> position of high symmetry point
                        in fractional coordinates of k-space""",
                        required=False, default=None)
    parser.add_argument('-d', '--sampling_distance', action='store', type=float,
                        help="""<optional> the k point sampling distance
                        in units 1/Angstrom """, metavar='0.01', 
                        required=False, default=None)
    parser.add_argument('-n', '--number_kpoints', type=int, metavar='50',
                       help="""<optional> the number of k points between 
                       gamma and high_symmetry_point """,
                       required=False, default=None)
    parser.add_argument('-m', '--scaling', type=float, metavar='5.235',
                        help="""<optional> scales input lattice vectors 
                        by given 
                        amount (defalt: 1)""", required=False, 
                        default=1)
    args = parser.parse_args()
    
    vectors = numpy.array(args.vectors)
    high_symmetry_point = numpy.array(args.high_symmetry_point)
    
    a_in = vectors[:3]
    b_in = vectors[3:6]
    c_in = vectors[6:9]
    
    values = transform_space(a_in, b_in, c_in, 
                    args.space, scaling=args.scaling)
    
    if args.high_symmetry_point is not None:
        calculate_distance(values['a_star_length'],
                           values['b_star_length'],
                           values['c_star_length'],
                           high_symmetry_point)
        if args.sampling_distance is not None:
            calculate_number_kpoints(args.sampling_density[0], 
                                     values['distance_to_hsp'])
        elif args.number_kpoints is not None:
            calculate_sampling_distance(args.number_kpoints[0], 
                                       values['distance_to_hsp'])
    print_to_screen(values)
