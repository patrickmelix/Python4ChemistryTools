#!/usr/bin/env python3

from ase import io
import numpy as np
from numpy.linalg import svd

def fit_plane(points):
    """Fit a plane to a set of 3D points.
    based on https://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    point, normal = fit_plane(points)

    Given an array, points, of shape (n_points, d)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    points = np.swapaxes(points, 0, 1)
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trailing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return ctr, svd(M)[0][:,-1]

def calc_dists(points, plane_point, plane_normal):
    """Calculate the distances between points and a plane.
    """
    dists = []
    for point in points:
        # get the vector from the plane to the point
        v = point - plane_point
        # project the vector onto the plane normal
        dists.append(np.abs(np.dot(v, plane_normal)))
        
    return dists

def main(file, plane_tag=None, ads_tag=0, dimension='c', direction='+', verbose=True):
    """Calculate the distance of an adsorbate to the above periodic slab.
    
    Parameters:
    -----------
    file: str
        Path to the structure file
        plane_tag: int
        Tag of the plane
        dimension: str
        Dimension of the plane
        
    Returns:
    --------
    float: distance
    """
    # read structure.
    structure = io.read(file)
    nAtoms = len(structure)

    if plane_tag is None:
        plane_tag = max(structure.get_tags())
        if verbose:
            print("Using tag {} for plane".format(plane_tag))


    dim = ['a','b','c'].index(dimension.lower())
    supercell = np.ones(3, dtype=int)
    supercell[dim] = 2
    structure_super = structure*supercell
    second_cell = structure_super[nAtoms:]

    # get the plane
    if direction == '-':
        plane_atoms = structure[[atom.index for atom in structure if atom.tag == plane_tag]]
    else:
        plane_atoms = second_cell[[atom.index for atom in second_cell if atom.tag == plane_tag]]
    plane_point, plane_normal = fit_plane(plane_atoms.positions)
    if verbose:
        print("Plane normal: {}".format(plane_normal))
        print("Plane point: {}".format(plane_point))

    # get adsorbate
    if direction == '-':
        adsorbate = second_cell[[atom.index for atom in second_cell if atom.tag == ads_tag]]
    else:
        adsorbate = structure[[atom.index for atom in structure if atom.tag == ads_tag]]

    # calculate distance
    dists = calc_dists(adsorbate.positions, plane_point, plane_normal)
    min_dist = np.min(dists)
    min_idx = np.argmin(dists)
    if verbose:
        print("Current cell vector length in this direction: {}".format(np.linalg.norm(structure.cell[dim])))
        print("Minimum distance: {}".format(min_dist))
        print("Minimum distance adsorbate atom index: {}".format(min_idx))
        print("Minimum distance adsorbate atom position: {}".format(adsorbate.positions[min_idx]))
    return min_dist

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Get the distance of an adsorbate to the above periodic slab.')
    parser.add_argument('file', type=str, help='ASE compatible structure file')
    parser.add_argument('--plane_tag', type=int, help='Tags of plane-defining atoms, None will use the largest tag', default=None)
    parser.add_argument('--ads_tag', type=int, help='Tags of adsorbate atoms', default=0)
    parser.add_argument('--dimension', type=str, help='Dimension of the vacuum', default='c')
    parser.add_argument('--direction', type=str, help='Direction of the vacuum', default='+')
    args = parser.parse_args()
    main(args.file, plane_tag=args.plane_tag, ads_tag=args.ads_tag, dimension=args.dimension, direction=args.direction, verbose=True)
