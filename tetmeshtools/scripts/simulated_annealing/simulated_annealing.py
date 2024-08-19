import argparse
import sys
import mrcfile
import pathlib
import tetmeshtools.mrclininterpolate as mi
from tetmeshtools.scripts.simulated_annealing.copy_mesh import read_tetgen_file, remove_suffix, write_tetgen_file
from tetmeshtools.scripts.simulated_annealing.compute_densities import compute_node_densities, get_density_at_point
from tetmeshtools.meshtools.tetgenstructs import NodePoint
import tetmeshtools.tetprops as tp
from tetmeshtools.meshtools.trisurface import TriSurface
import numpy as np
import random
import time
import math
import pandas

def get_args():
    description = "Performs simulated annealing given an MRC file and a tetmesh file."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-m",
                        "--mrcfile",
                        type=pathlib.Path,
                        required=True,
                        help="Input MRC file")
    parser.add_argument("-i",
                        "--tetmeshfile",
                        type=pathlib.Path,
                        required=True,
                        help="Input tetmesh file")
    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="Output file name root, (no suffix)")
    parser.add_argument("--density_weight",
                        "-A",
                        type=float,
                        default=1,
                        required=False,
                        help="Sum of Density Weights")
    parser.add_argument("--surface_area_weights",
                        "-B",
                        type=float,
                        default=1,
                        required=False,
                        help="Weights for variance of triangle surface area")
    parser.add_argument("--mutation_probability",
                        type=float,
                        required=True,
                        help="Probability of mutating a node.")
    parser.add_argument("--mutation_multiplier",
                        type=float,
                        required=True,
                        help="Multiplier controls how large mutation changes can be.")
    parser.add_argument("-t",
                        "--temperature",
                        "--initial_temperature",
                        type=float,
                        required=True,
                        help="Controls the initial temperature.")
    parser.add_argument("--temperature_decrement",
                        type=float,
                        required=False,
                        default=1,
                        help="Decrement of temperature.")
    parser.add_argument("--iterations_per_temperature",
                        type=int,
                        required=True,
                        help="Determines how many iterations are done at each temperature change.")
    parser.add_argument("-s",
                        "--seed",
                        type=int,
                        required=False,
                        help="Seed to control randomness.")
    parser.add_argument("--target_density_sum",
                        type=float,
                        required=True,
                        help="Target density sum during the annealing process.")
    parser.add_argument("-k",
                        type=float,
                        required=False,
                        default=1,
                        help="Multiplier against T for checking probability of keeping worse fitness.")
    parser.add_argument("-p",
                        "--progress",
                        action="store_true",
                        help="Print optimization progress (may slow down optimization)")

    return parser.parse_args()

def get_sum_of_density_error(node_densities, target_density):
    """
    This function returns the sum of density error given a list of values.

    Parameters
    ----------
    node_densities : list
        A list of all the node densities.
    target_density : float
        Density to be targeted
    
    Returns
    -------
    float
        a number representing the sum of density error.
    """
    sum_of_density_error = 0
    for density in node_densities:
        sum_of_density_error += abs(target_density - density)
    return sum_of_density_error

def get_surface_area_variance(surface_area_values):
    """
    This function returns the surface area variance of the mesh.

    Parameters
    ----------
    surface_area_values : list
        A list of all face areas.
    
    Returns
    -------
    float
        Surface area variance of the solution.
    """
    return np.var(surface_area_values)

def fitness_function(node_densities, surface_area_values, target_density, A=1, B=1):
    """
    This function returns fitness score of the current solution.

    @TODO Add more parameters in order to optimize tetrahedrons.

    Parameters
    ----------
    node_densities : list
        A list of all node densities.
    surface_area_values : list
        A list of all face areas.
    target_density : float
        Target density of nodes.
    A : float, optional
        Weighting for the sum of node density error. (Default is 1)
    B : float, optional
        Weighting for the surface area variance. (Default is 1)
    
    Returns
    -------
    float
        Fitness score of the current solution.
    """
    sum_of_density_error = get_sum_of_density_error(node_densities, target_density=target_density)
    surface_area_variance = get_surface_area_variance(surface_area_values)
    return (A * sum_of_density_error) + (B * surface_area_variance)
        

def mutate(node, temperature, mutation_multiplier=1):
    """
    This function takes a node and creates a mutated copy of it.

    Parameters
    ----------
    node : NodePoint
        A node to be mutated.
    temperature : float
        Current temperature.
    mutation_multiplier : float, optional
        Scale of mutation. (Default is 1)
    
    Returns
    -------
    NodePoint
        New node that has been mutated.
    """
    
    mutation_max_value = mutation_multiplier * temperature
    new_node = NodePoint(node.index,
        node.x + random.uniform(-mutation_max_value, mutation_max_value),
        node.y + random.uniform(-mutation_max_value, mutation_max_value),
        node.z + random.uniform(-mutation_max_value, mutation_max_value)
    )
    return new_node # Maybe also return how close it was to max value (might help with skipping temperature for faster processing)

def compute_probability_of_keeping(current_fitness, new_fitness, temperature, k=1):
    """
    This function returns the probability of accepting a worse solution
    based on the delta fitness, temperature, and the value of k.

    Parameters
    ----------
    current_fitness : float
        Fitness of the current solution.
    new_fitness : float
        Fitness of the new solution.
    temperature : float
        Current temperature.
    k : float, optional
        Weight to adjust range of accepted probabilities. (Default is 1)
    
    Returns
    -------
    float
        a number representing the probability of acceptance.
    """
    if new_fitness == current_fitness:
        return 1
    return math.exp((current_fitness - new_fitness)/ (temperature * k))

def metropolis_algorithm(current_fitness, new_fitness, temperature, k=1):
    """
    This function returns whether the new solution should be accepted
    using the metropolis algorithm.

    Parameters
    ----------
    current_fitness : float
        Fitness of the current solution.
    new_fitness : float
        Fitness of the new solution.
    temperature : float
        Current temperature.
    k : float, optional
        Weight to adjust range of accepted probabilities. (Default is 1)
    
    Returns
    -------
    tuple of bool
        A tuple containing two booleans:
            - First boolean represents whether the solution was accepted.
            - Second boolean represents whether it is a worse solution that was accepted.
    """

    accept_change = current_fitness > new_fitness
    worse_solution_accepted = False
    if not accept_change:
        keep_probability = compute_probability_of_keeping(current_fitness, new_fitness, temperature, k=k)
        if random.random() <= keep_probability:
            accept_change, worse_solution_accepted = True, True
    return accept_change, worse_solution_accepted



def simulated_annealing(
        mrc_image,
        nodes,
        faces,
        target_density,
        A=1,
        B=1,
        k = 1,
        mutation_probability=1,
        mutation_multiplier = 1,
        iterations_per_temp = 100,
        initial_temperature = 100,
        temperature_decrement = 1,
        print_progress = False,
        seed = None
    ):
    """
    Performs simulated annealing on the nodes & faces of a mesh.

    Parameters
    ----------
    mrc_image : MRCImage
        The MRC image of the protein.
    nodes : dict
        Dictionary of all nodes.
    faces : dict
        Dictionary of all faces.
    target_density : float
        Target density of nodes.
    A : float, optional
        Weighting for the sum of node density error. (Default is 1)
    B : float, optional
        Weighting for the surface area variance. (Default is 1)
    k : float, optional
        Weight to adjust range of accepted probabilities. (Default is 1)
    mutation_probability : float, optional
        Probability of mutating a node. (Default is 1)
    mutation_multiplier : float, optional
        Scale of mutation. (Default is 1)
    iterations_per_temp : int, optional
        Iterations (Steps) per temperature. (Default is 100)
    initial_temperature : float, optional
        Initial temperature for the algorithm. (Default is 100)
    temperature_decrement : float, optional
        Temperature decrement for the algorithm. (Default is 1)
    print_progress : boolean, optional
        Print algorithm progress on the console. (Default is False)
    seed : int?, optional
        Seed for randomization. (Default is None)
    
    Returns
    -------
    tuple of (dict, pandas.DataFrame)
        A tuple containing two elements:
            - dict : A dictionary of all the new mutated nodes.
            - pandas.DataFrame : Dataframe documenting algorithm progress over each temperature.
    """
    
    # Initialize Seed
    random.seed(seed)

    # Compute initial density values
    node_densities = compute_node_densities(mrc_image, nodes)

    # If one of the nodes does not have a density, terminate.
    if None in node_densities.values():
        print("One or more node(s) is not fitted to the MRC Image. The optimization cannot continue.")
        exit()

    # Compute initial surface area values
    triangle_surface_areas = dict()
    tri_surface = TriSurface(nodes, faces)
    for face_index in tri_surface.get_faces().keys():
        face_nodes = tri_surface.get_triangle_nodes(face_index)
        triangle_surface_areas[face_index] = (tp.area_of_triangle(face_nodes))

    surface_area_list = list(triangle_surface_areas.values())

    current_fitness = fitness_function(node_densities.values(), surface_area_list, target_density, A=A, B=B)
    current_node = None

    # Print initial values.
    if print_progress:
        sum_of_density_error = get_sum_of_density_error(node_densities.values(), target_density=target_density)
        surface_area_variance = get_surface_area_variance(surface_area_list)
        print("Initial Raw Density Error and Surface Area Variance:      ", sum_of_density_error, surface_area_variance)
        print("Initial Weighted Density Error and Surface Area Variance: ", (A * sum_of_density_error),(B * surface_area_variance))
        print("Initial Fitness: ", current_fitness)

    # Start logging
    df = pandas.DataFrame({
        "Temperature": [],
        "Avg Surface Area": [],
        "Surface Area Variance": [],
        "Avg Node Density": [],
        "Node Density Variance": [],
        "Node Density Error": []
    })

    # Start
    start_time = time.time()
    try:
        temperature = initial_temperature
        while temperature > 0:
            if print_progress:
                print(f"Current Temperature: {temperature}")
            for i in range(iterations_per_temp): # Maybe change it to not be a value
                successful_mutations = 0
                worse_solutions_accepted = 0
                worse_solutions_total = 0
                old_fitness = current_fitness
                for node in nodes.values():
                    if random.random() <= mutation_probability:
                        current_node = node

                        # mutate node
                        new_node = mutate(node, temperature, mutation_multiplier=mutation_multiplier)

                        # Calculate new density
                        new_density = get_density_at_point(mrc_image, new_node.x, new_node.y, new_node.z)
                        if new_density is None:
                            continue

                        # Save node as replacement for original (Might be swapped back later)
                        nodes[node.index] = new_node
                        
                        old_triangle_areas = dict()
                        for face in faces.values():
                            # Check if a face is associated with a node
                            if new_node.index in {face.vert0, face.vert1, face.vert2}:
                                face_nodes = tri_surface.get_triangle_nodes(face.index)
                                area = tp.area_of_triangle(face_nodes)                                    
                                # Cache old triangle area
                                old_triangle_areas[face.index] = triangle_surface_areas[face.index]
                                # Apply new area
                                triangle_surface_areas[face.index] = area

                                # If inverted, something has gone wrong and the algorithm should terminate.
                                if (np.sign(area) != np.sign(old_triangle_areas[face.index])):
                                    print("inverted.")
                                    exit()

                        # Cache old density
                        old_density = node_densities[node.index]

                        # Apply new density
                        node_densities[node.index] = new_density

                        # Calculate fitness of the new solution
                        new_fitness = fitness_function(node_densities.values(), list(triangle_surface_areas.values()), target_density, A=A, B=B)

                        # Metropolis algorithm to determine acceptance
                        accept_change, accepted_worse_solution = metropolis_algorithm(current_fitness, new_fitness, temperature, k=k)

                        if accept_change:
                            if accepted_worse_solution:
                                worse_solutions_accepted += 1
                                worse_solutions_total += 1
                            current_node = None
                            current_fitness = new_fitness
                            successful_mutations += 1
                        else:
                            worse_solutions_total += 1
                            nodes[node.index] = node
                            # revert cached density and triangle area
                            node_densities[node.index] = old_density
                            for face_index, area in old_triangle_areas.items():
                                triangle_surface_areas[face_index] = area
                if print_progress:
                    if successful_mutations > 0:
                        if worse_solutions_total > 0:
                            percentage_worse_solutions_accepted = round(worse_solutions_accepted/worse_solutions_total*100,1)
                        else:
                            percentage_worse_solutions_accepted = 0
                        print(f"New Best Fitness: {current_fitness} (Improvement: {old_fitness-current_fitness}) at iteration {i} (temperature: {temperature}) (successful mutations: {successful_mutations}) (worse solutions accepted: {worse_solutions_accepted} ({percentage_worse_solutions_accepted}%)) (Time Since Start: {round(time.time() - start_time, 2)} seconds)")
                        sum_of_density_error = get_sum_of_density_error(node_densities.values(), target_density=target_density)
                        surface_area_variance = get_surface_area_variance(list(triangle_surface_areas.values()))
                        print("Raw Density Error and Surface Area Variance:      ", sum_of_density_error, surface_area_variance)
                        print("Weighted Density Error and Surface Area Variance: ", (A * sum_of_density_error),(B * surface_area_variance))

            
            # Log results
            area_values = np.array(list(triangle_surface_areas.values()))
            density_values = np.array(list(node_densities.values()))
            df.loc[len(df.index)] = [
                temperature,
                area_values.mean(),
                np.var(area_values),
                density_values.mean(),
                np.var(density_values),
                get_sum_of_density_error(node_densities.values(), target_density=target_density)
            ]

            # Decrement temperature
            temperature -= temperature_decrement
    except KeyboardInterrupt:
        if current_node is not None:
            # Swap back for safety measures
            nodes[current_node.index] = current_node
        print(f"Early shutdown.")
    except Exception as E:
        print("An error has occured. Early Shutdown.")
        # print(E)
        import traceback
        print(traceback.format_exc())
    if print_progress:
        print(f"Finished in {round(time.time() - start_time, 2)} seconds.")
        sum_of_density_error = get_sum_of_density_error(node_densities.values(), target_density=target_density)
        surface_area_variance = get_surface_area_variance(list(triangle_surface_areas.values()))
        print("Final results:")
        print("Raw Density Error and Surface Area Variance:      ", sum_of_density_error, surface_area_variance)
        print("Weighted Density Error and Surface Area Variance: ", (A * sum_of_density_error),(B * surface_area_variance))
        print(f"Fitness score: {current_fitness}")
    return nodes, df

def command_validation(args):
    """
    Validates arguments provided through the commandline.
    """
    # Verify that the inputs exist
    if not args.mrcfile.exists():
        return f"Error: file {args.mrcfile} does not exist!"
    if not args.tetmeshfile.exists():
        return f"Error: file {args.tetmeshfile} does not exist!"
    if args.mutation_probability < 0:
        return f"Error: Mutation probability must be greater than 0."
    if args.mutation_probability > 1:
        return f"Error: Mutation probability must be less than 1."
    if args.iterations_per_temperature < 1:
        return f"Error: Iterations per temperature must be greater than 0."
    if args.temperature <= 0:
        return f"Error: Initial temperature must be greater than 0."
    return None

def main():
    # Verify arguments
    args = get_args()
    error_message = command_validation(args)
    if error_message is not None:
        print(error_message, file=sys.stderr)
        return
    
    # load tetgen mesh 
    nodes, faces, tets, _ = read_tetgen_file(pathlib.Path(remove_suffix(str(args.tetmeshfile))))
    path = pathlib.Path(str(args.mrcfile))

    with mrcfile.mmap(path, mode='r+') as mrc:
        mrc_image = mi.MRCImage(mrc)

    # Assign parameters
    A = args.density_weight
    B = args.surface_area_weights
    mutation_probability= args.mutation_probability
    mutation_multiplier = args.mutation_multiplier
    target_density=args.target_density_sum
    iterations_per_temp = args.iterations_per_temperature
    initial_temperature = args.temperature
    temperature_decrement = args.temperature_decrement
    seed = args.seed
    k = args.k
    print_progress = args.progress
    
    # Perform mesh optimization
    nodes, df = simulated_annealing(
        mrc_image,
        nodes,
        faces,
        A=A,
        B=B,
        target_density=target_density,
        mutation_multiplier=mutation_multiplier,
        mutation_probability=mutation_probability,
        initial_temperature=initial_temperature,
        iterations_per_temp= iterations_per_temp,
        seed=seed,
        k=k,
        temperature_decrement=temperature_decrement,
        print_progress=print_progress
    )
    # Save evolution file
    df.to_csv(str(args.output)+'_evolution.csv', index=False)
    comment = f"Simulated Annealing with: A={A}, B={B}, k={k}, mutation_probability={mutation_probability}, target_density={target_density} mutation_multiplier={mutation_multiplier}, iterations_per_temp={iterations_per_temp}, initial_temperature={initial_temperature}, temperature_decrement={temperature_decrement}"
    write_tetgen_file(pathlib.Path(str(args.output)), nodes, faces, tets, comment=comment)
    if print_progress:
        print("Saved Model.")

if __name__ == "__main__":
    main()